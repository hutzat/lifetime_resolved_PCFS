%% t3PCFS.m V1.0 @ HENDRIK UTZAT 2017
%  29/10/2017
%%Class definition for analysis of pulsed photon resolved PCFS files as generates
%from the Labview instrument control software by HENDRIK UTZAT 

% compiles the lifetimes of all photon-streams first, then sorts them and
% creates PCFS children for each tau-limits.  



classdef t3PCFS_simul<dynamicprops
    properties
        folder=[];
        buffer_size=[];
        pico_path=[];
        mode=[];
        PCFS_FID=[];
    end
    methods
        function obj=t3PCFS_simul(PCFS_folder,buffer_size)
            %constructor method which is called when a PCFS object is created.
            %save the fundamental arguments as object properties.
            %mode is either t2 or t3
            
            %PCFS_FOLDER is the full path to the folder contaiing the
            %photon stream data and meta data of the PCFS run.
            %Buffer_size is number of photon records to keep in memory -
            %default 1E6.
            %pico_path is the system installation folder of Tom Bischof's
            %photon_intensity_correlate code - ONLY WORKS FOR MAC OS.
            
            obj.folder=PCFS_folder;
            PCFS_folder = strrep(PCFS_folder,'\','/');
            PCFS_FID=strsplit(PCFS_folder,'/');
            obj.PCFS_FID=strrep(strcat(PCFS_FID(end-2),'/',PCFS_FID(end-1)),'_',' ');
            obj.buffer_size=buffer_size;
        
            
            %get list of all files in the PCFS directory.
            files=dir;
            dummy=addprop(obj,'files');
            obj.files=files;
            
            %Read in the metadata of the .pcfslog file and store it as property.
            for i=3:numel(files)
                [path,f_name,ext]=fileparts(files(i).name);
                
                %read in the metadata from the pcfs logfile.
                if strcmp(ext,'.pcfslog')== true;
                    
                    fid = fopen(files(i).name);
                    
                    tline = fgetl(fid);
                    
                    while ischar(tline)
                        
                        semicolon=find(tline == '=')
                        
                        if ~isempty(semicolon) == true;
                            prop_name=tline(1:(semicolon-1));
                            dummy=addprop(obj,(prop_name));
                            obj.(prop_name)=str2double(tline((semicolon+1):end));
                            disp(tline)
                        end
                        tline = fgetl(fid);
                    end
                    
                    fclose(fid);
                    
                end
                
                %read in the position file as array.
                if strcmp(ext,'.pos')== true
                    dummy=addprop(obj,'stage_positions');
                    obj.stage_positions=textread(files(i).name);
                end
                
                %create photons object for the photon stream at each
                %interferometer path length difference.
                if strcmp(ext,'.photons')== true
                    dummy=addprop(obj,(f_name));
                    obj.(f_name)=photons('dummy.stream',obj.buffer_size);
                end
            end
        end
        
        %% get and parse photonstream/get correlation functions
        function obj=get_photons_all(obj)
            %just batch-run picoquant_bin of the photon objects to get the photon arrival data.
            
            files=dir;
            obj.files=files;
            for i=3:numel(obj.files)
                [path,f_name,ext]=fileparts(obj.files(i).name);
                if strcmp(ext,'.stream')== true
                    obj.(f_name).get_photon_records()
                end
            end
        end
        function obj=get_sum_signal_all(obj)
           %%get the auto_correlation of the sum signal of the two
            %%detectors for all photon arrival files.
            files=dir;
            obj.files=files;
            for i=3:numel(obj.files)
                [path,f_name,ext]=fileparts(obj.files(i).name);
                if strcmp(ext,'.photons')== true
                    if strcmp(f_name(1:3),'sum')==0;
                        obj.(f_name).write_photons_to_one_channel(f_name,strcat('sum_signal_',f_name))
                    end
                end
            end
        
        end
        function obj=get_all_lifetimes(obj,resolution)
            files=dir;
            obj.files=files;
            for i=3:numel(obj.files)
                [path,f_name,ext]=fileparts(obj.files(i).name);
                if strcmp(ext,'.photons')== true
                    if strcmp(f_name(1:3),'sum')==true %looking to lifetimes of the sum signals.
                        obj.(f_name(12:end)).lifetime_histo(f_name(12:end),'life',resolution)
                    end
                end
            end
        end
        function obj=fit_all_lifetimes(obj)
            %% call fit_lifetime_biexponential for all original photon-streams of the PCFS run.
            % initiate some arrays to capture the fit parameters.
            
            a=zeros(numel(obj.stage_positions),1);
            b=zeros(numel(obj.stage_positions),1);
            c=zeros(numel(obj.stage_positions),1);
            d=zeros(numel(obj.stage_positions),1);
            
            files=dir;
            obj.files=files;
            for i=3:numel(obj.files);
                [path,f_name,ext]=fileparts(obj.files(i).name);
                if strcmp(ext,'.stream')== true;
                    obj.(f_name).fit_lifetime_biexponential('life');
                    k=-1;
                    r=[];
                    while true
                        k=k+1;
                        if isstrprop(f_name(end-k),'digit')==true;
                            r=[f_name(end-k),r];
                        else
                            correlation_number=str2num(r);
                            r=[];
                            break
                        end
                    end
                    
                    
                    a_temp=obj.(f_name).life.exp2.fit.a;
                    b_temp=-1./obj.(f_name).life.exp2.fit.b;
                    c_temp=obj.(f_name).life.exp2.fit.c;
                    d_temp=-1./obj.(f_name).life.exp2.fit.d;
                    
                    if b_temp>d_temp
                        a(correlation_number+1)=a_temp;
                        b(correlation_number+1)=b_temp;
                        c(correlation_number+1)=c_temp;
                        d(correlation_number+1)=d_temp;
                    else
                        a(correlation_number+1)=c_temp;
                        b(correlation_number+1)=d_temp;
                        c(correlation_number+1)=a_temp;
                        d(correlation_number+1)=b_temp;
                        
                    end
                end
            end
            if ~isprop(obj,'lifetime_fits')
                dummy=addprop(obj,'lifetime_fits');
            end
            
            obj.lifetime_fits=struct('a',a,'b',b,'c',c,'d',d);
            
        end
        function obj=parse_photons(obj,tau_id,tau)
            % sorts the photon-stream at each state position according to photon-arrival time
            % between tau: [lower,upper]
            % we do that for both the sum_signals and the original
            % photon-stream.
            % tau is a cell array containig the fids for all the different
            % sorting parameters.
            files=dir;
            obj.files=files;
            for k=1:numel(tau_id);
            for i=3:numel(obj.files)
                [path,f_name,ext]=fileparts(obj.files(i).name);
                if strcmp(ext,'.photons')== true;
                    if strcmp(f_name(1:3),'tau')==0;
                        if strcmp(f_name(1:3),'sum')==0;
                            obj.(f_name).arrival_time_sorting(f_name,strcat(char(tau_id(k)),f_name),tau(k,:))
                        end
                          if strcmp(f_name(1:3),'sum')==1;
                            obj.(f_name(12:end)).arrival_time_sorting(f_name,strcat(char(tau_id(k)),f_name),tau(k,:))
                        end
                    end
                end
            end
            end
        end
        function obj=get_all_intensity_traces(obj,bin_width)
            files=dir;
            obj.files=files;
            for i=3:numel(obj.files)
                [path,f_name,ext]=fileparts(obj.files(i).name);
                if strcmp(ext,'.photons')== true
                    if strcmp(f_name(1:3),'sum')==true %looking to lifetimes of the sum signals.
                        
                        if isprop(obj,(f_name(12:end)))
                            obj.(f_name(12:end)).get_intensity(f_name,'intensity',bin_width)
                        end
                    end
                end
            end
            
        end
        function obj=get_all_intensity_traces_fid(obj,bin_width,fid)
            files=dir;
            obj.files=files;
            for i=3:numel(obj.files)
                [path,f_name,ext]=fileparts(obj.files(i).name);
                if strcmp(ext,'.photons')== true
                    if strcmp(f_name(1:3),'sum')==true %looking to lifetimes of the sum signals.
                        if strcmp(f_name(11:(10+length(fid))),fid)==true
%                        if isprop(f_name((11+length(fid)):end))
                            disp((f_name((11+length(fid)):end)))
                            obj.(f_name((11+length(fid)):end)).get_intensity(f_name,'intensity',bin_width)
                      %  end
                        end
                        
                    end
                end
            end
        end
        
        
        function obj=plot_all_intensity_traces(obj)
            props = properties(obj);
            figure()
            count=1
            for iprop = 1:length(props)
                if isa(obj.(props{iprop}),'photons')==true
                    subplot(10,10,count)
                    plot(obj.(props{iprop}).intensity.time,obj.(props{iprop}).intensity.trace)
                    hold on
                    count=count+1;
                end
            end
        end
        function obj=get_all_intensity_traces_for_Fourier(obj,fid,bin_width)
            files=dir;
            obj.files=files;
            lstr=length(fid);
            for i=3:numel(obj.files)
                [path,f_name,ext]=fileparts(obj.files(i).name);
                if strcmp(ext,'.photons')== true
                    if strcmp(f_name(1:3),'sum')==false
                        if strcmp(f_name(1:lstr),fid)==true%looking to lifetimes of the sum signals.
                            obj.(f_name((lstr+1):end)).get_intensity(f_name,strcat(fid,'_intensity'),bin_width)
                        end
                    end
                end
            end
            
        end
        function obj=sort_intensity_traces(obj,fid,bin_width,bounds)
            % bounds is [lower,upper;lower,upper;.....] for the different
            % stage positions.
            
            files=dir;
            obj.files=files;
            for i=3:numel(obj.files)
                
                [path,f_name,ext]=fileparts(obj.files(i).name);
                if strcmp(ext,'.photons')== true
                    if strcmp(f_name(1:3),'sum')==false
                        if  isprop(obj,f_name)==true
                            
                            % extract correlation number.
                            k=-1;
                            r=[];
                            while true
                                k=k+1;
                                if isstrprop(f_name(end-k),'digit')==true;
                                    r=[f_name(end-k),r];
                                else
                                    correlation_number=str2num(r);
                                    r=[];
                                    break
                                end
                            end
                            
                            lower_limit=bounds(correlation_number+1,1);
                            upper_limit=bounds(correlation_number+1,2);
                            obj.(f_name).get_intensity(f_name,strcat(fid,f_name),bin_width,lower_limit, upper_limit)
                            
                        end
                    else % looking at sum signal also to be intensity sorted.
                        if  isprop(obj,f_name(12:end))==true
                            
                            % extract correlation number.
                            k=-1;
                            r=[];
                            while true
                                k=k+1;
                                if isstrprop(f_name(end-k),'digit')==true;
                                    r=[f_name(end-k),r];
                                else
                                    correlation_number=str2num(r);
                                    r=[];
                                    break
                                end
                            end
                            
                            lower_limit=bounds(correlation_number+1,1);
                            upper_limit=bounds(correlation_number+1,2);
                            obj.(f_name(12:end)).get_intensity(f_name(12:end),strcat('sum_signal',fid,f_name(12:end)),bin_width,lower_limit, upper_limit)
                        end
                    end
                end
            end
            
        end
        function obj=divide_photons(obj,tau)
            %tau is a vector that contains the boundaries for the
            %photon-compartments for photon-sorting such that for tau=[x,y,z,...]
            % the photon stream is parsed into x<tau<y, y<tau<z ....
            n_children=numel(tau)-1; %number of iterations for photon-parsing.
            for i=1:n_children
                obj.parse_photons(['tau',num2str(i),'_'],[tau(i),tau(i+1)])
            end
            if ~isprop(obj,'n_children')
                dummy=addprop(obj,'n_children');
            end
            obj.n_children=n_children;
            
            
            if ~isprop(obj,'tau_bounds')
                dummy=addprop(obj,'tau_bounds');
            end
            obj.tau_bounds=tau;
            
        end
        function obj=get_all_blinking_statistics(obj,intensity_ID,boundaries)
            
            if ~isprop(obj,'counts_length_blinking')
               dummy=addprop(obj,'counts_length_blinking')
            end
           
            
            %get list of files.
            files=dir;
            obj.files=files;
            
          %  l_str=numel(fid);
            
            counts_g=[];
            counts_b=[];
            counts_d=[];
            for i=3:numel(obj.files) %exclusing . and ..
                [path,f_name,ext]=fileparts(obj.files(i).name);
                if strcmp(ext,'.photons')== true
                    if strcmp(f_name(1:13),'sum_signalall')==1; %looking at only the streams that contain all photons.
                        obj.(f_name(14:end)).get_blinking_analysis(intensity_ID,boundaries)
                        counts_g=[counts_g,obj.(f_name(14:end)).counts_length_intensity.count_grey_length];
                        counts_d=[counts_d,obj.(f_name(14:end)).counts_length_intensity.count_bright_length];
                        counts_b=[counts_b,obj.(f_name(14:end)).counts_length_intensity.count_dark_length];
                        
                    end
                end
            end
            
            obj.counts_length_blinking=struct('grey',counts_g,'dark',counts_d,'bright',counts_b);
        end

        function obj=get_all_FLID(obj,file_out_key,resolution,pulses_per_bin,Sync_Rate)
            files=dir;
            obj.files=files;
            for i=3:numel(obj.files)
                [path,f_name,ext]=fileparts(obj.files(i).name);
                if strcmp(ext,'.photons')== true
                    if strcmp(f_name(1:3),'sum')==true
                        if strcmp(f_name(11:13),'all')==true
                        obj.(f_name(14:end)).FLID(f_name(14:end),file_out_key,resolution,pulses_per_bin,Sync_Rate);
                        end
                    end
                    end
                end
        end
        function obj=create_PCFS_children(obj,time_bounds,lag_precision,Sync_Rate)
            
            n_children=obj.n_children; %number of iterations for photon-parsing.
            multi_dim_auto=[];
            multi_dim_cross=[];
            multi_dim_PCFS_inter=[];
            
            for i=1:n_children
                % create a PCFS child with the sorted photon streams and
                % call the get_correlation functions with the respective file
                % identifiers.
                
                fid=strcat('tau',num2str(i),'_');
                if ~isprop(obj,fid)
                    dummy=addprop(obj,fid);
                end
                
                disp(fid)
                command=strcat('obj.',fid,'=PCFS_simul(','''',obj.folder,'''',',','2E6,','''','t3','''',')');
                eval(command);
                
                % now call the PCFS function that calculates the auto-corr
                % and cross-corr and populates the respective matrices.
                
                obj.(fid).get_intensity_corrs_fid('t3',time_bounds,lag_precision,fid,Sync_Rate);
                
                
                %% write  multidimensional matrices containing the interferograms.
                multi_dim_auto(:,:,i)=obj.(fid).auto_corr_sum_interferogram;
                multi_dim_cross(:,:,i)=obj.(fid).cross_corr_interferogram;
                multi_dim_PCFS_inter(:,:,i)=obj.(fid).PCFS_interferogram;
            end  
            
            if ~isprop(obj,'multi_dim_auto')
                    dummy=addprop(obj,'multi_dim_auto');
                end
                
                if ~isprop(obj,'multi_dim_cross')
                    dummy=addprop(obj,'multi_dim_cross');
                end
                
                if ~isprop(obj,'multi_dim_PCFS_inter')
                    dummy=addprop(obj,'multi_dim_PCFS_inter');
                end
                
                obj.multi_dim_PCFS_inter=multi_dim_PCFS_inter;
                obj.multi_dim_cross=multi_dim_cross;
                obj.multi_dim_auto=multi_dim_auto; 
                
        end
        
        function obj=create_single_PCFS_child(obj,fid,time_bounds,lag_precision,mode,Sync_Rate)
               
                
                if ~isprop(obj,fid)
                    dummy=addprop(obj,fid);
                end
                switch mode
                    case 2
                command=strcat('obj.',fid,'=PCFS_simul(','''',obj.folder,'''',',','2E6,','''','t2','''',')');
                eval(command);
                
                % now call the PCFS function that calculates the auto-corr
                % and cross-corr and populates the respective matrices.
                obj.(fid).get_intensity_corrs_fid('t2',time_bounds,lag_precision,fid,Sync_Rate);
                
                    case 3
                         command=strcat('obj.',fid,'=PCFS_simul(','''',obj.folder,'''',',','2E6,','''','t3','''',')');
                eval(command);
                
                % now call the PCFS function that calculates the auto-corr
                % and cross-corr and populates the respective matrices.
                obj.(fid).get_intensity_corrs_fid('t3',time_bounds,lag_precision,fid,Sync_Rate);
                
                end
        
        end
        function obj=get_intensity_correlations(obj,mode,time_bounds,lag_precision)
            % obtains the cross-correlations and the auto-correlation 
            % of the sum signal at each stage position and returns them 
            % in a matrix obj.cross_corr_interferogram where the first line
            % is the obj.stage_positions. Also returns the auto-correlation
            % function of the sum signal to a matirx of similar structure. 
            
            
            if ~isprop(obj,'time_bounds')
                dummy=addprop(obj,'time_bounds')
            end
            obj.time_bounds=time_bounds;
            
            if ~isprop(obj,'lag_precision')
                dummy=addprop(obj,'lag_precision')
            end
            obj.lag_precision=lag_precision;
            
            
            
            %get list of files.
            files=dir;
            obj.files=files;
            
            for i=3:numel(obj.files) %exclusing . and ..
                [path,f_name,ext]=fileparts(obj.files(i).name);
                if strcmp(ext,'.photons')== true
                    if strcmp(f_name(1:3),'sum')==false %looking to get cross-corrs
                        
                        %%get_cross correlation and do stuff with it.
                       
                         obj.(f_name).photon_corr(f_name,'cross_corr',[0,1], time_bounds, lag_precision, 0, 0);
                         tau=obj.(f_name).cross_corr.lags;
                         l_tau=length(tau);
                        
                        if ~isprop(obj,'cross_corr_interferogram')
                            dummy=addprop(obj,'cross_corr_interferogram')
                        end
                        
                        %%create an array containing to be filled up with the PCFS
                        %%interferogram.
                        if exist('cross_corr_interferogram')==0
                            cross_corr_interferogram=zeros((l_tau+1),(length(obj.stage_positions)));
                        end
                        
                        cross_corr_interferogram(1,:)=obj.stage_positions;
                      
                        %extract the number of the correlation
                        %measurement from the filename.
                        k=-1;
                        r=[];
                        while true
                            k=k+1;
                            if isstrprop(f_name(end-k),'digit')==true;
                                r=[f_name(end-k),r];
                            else
                                correlation_number=str2num(r);
                                r=[];
                                break
                            end
                        end
                        
                        cross_corr_interferogram(2:end,(correlation_number+1))=obj.(f_name).cross_corr.corr_norm;
                        obj.cross_corr_interferogram=cross_corr_interferogram;

                    else %means that we are looking at the sum signal photon stream of both detectors
                        
                        %%get and do stuff with the autocorrelation of the sum signal.
                      
                         obj.(f_name(12:end)).photon_corr(f_name,'auto_corr',[0,0], time_bounds, lag_precision, 0, 0);
                         tau=obj.(f_name(12:end)).cross_corr.lags;
                         l_tau=length(tau);                     
                        
                        
                        if ~isprop(obj,'auto_corr_sum_interferogram')
                            dummy=addprop(obj,'auto_corr_sum_interferogram')
                        end
                        
                        if ~isprop(obj,'tau')
                            dummy=addprop(obj,'tau')
                        end
                       
                        %%create an array containing to be filled up with the PCFS
                        %%interferogram.
                        if exist('auto_corr_sum')==0
                            auto_corr_sum=zeros((l_tau+1),(length(obj.stage_positions)));
                        end
                        
                        %fill in the stage positions in top line.
                        auto_corr_sum(1,:)=obj.stage_positions;
                        
                        %extract the number of the correlation
                        %measurement from the filename.
                        k=-1;
                        r=[];
                        while true
                            k=k+1;
                            if isstrprop(f_name(end-k),'digit')==true;
                                r=[f_name(end-k),r];
                            else
                                correlation_number=str2num(r);
                                r=[];
                                break
                            end
                        end
                        
                        auto_corr_sum(2:end,(correlation_number+1))=obj.((f_name(12:end))).auto_corr.corr_norm;
                        obj.auto_corr_sum_interferogram=auto_corr_sum;
                        obj.tau=obj.((f_name(12:end))).auto_corr.lags;
                        
                        % substract auto-correlation of sum signal from the
                        % cross correlation. 
                        PCFS_interferogram=cross_corr_interferogram-auto_corr_sum;
                        PCFS_interferogram(1,:)=obj.stage_positions;
                        
                        
                        if ~isprop(obj,'PCFS_interferogram')
                            dummy=addprop(obj,'PCFS_interferogram')
                        end
                        obj.PCFS_interferogram=PCFS_interferogram;
                        
                    end
                end
            end
        end
        
        
        %% wrangle and plot PCFS-type data
        function obj=plot_interferogram(obj,tau_select)
          
            tau=obj.tau;
            
            %Plot raw interferogram at different tau.
            legend_str=[];
            figure();
            for i=1:length(tau_select);
                [dummy,index]=min((tau- tau_select(i)).^2);
                plot(obj.stage_positions,obj.cross_corr_interferogram(index+1,:));
                legend_str=[legend_str; {num2str(tau_select(i))}];
                hold on
            end
            xlabel('Stage Position [mm]')
            ylabel('g^{(2)}_{cross}')
            columnlegend(2,legend_str, 'Location', 'NorthWest', 'boxoff')
            title({char(obj.PCFS_FID),' Uncorrected PCFS Interferogram at Different Tau'})
            cm=colormap(jet(length(obj.stage_positions)));
            
            
            %% plot the cross-correlation functions.
            legend_str=[];

            figure()
            subplot(2,1,1)
            for i=1:(length(obj.stage_positions))
                h{i}=semilogx(tau,obj.cross_corr_interferogram(2:end,i),'color',cm(i,:))
                legend_str=[legend_str;{num2str(obj.stage_positions(i))}];
                hold on
            end
            %xlabel('\tau [pulses or ps]')
            ylabel('g^{(2)}_{cross}')
            legend([h{1},h{length(obj.stage_positions)}],{num2str(obj.stage_positions(1)),num2str(obj.stage_positions(end))})
            %columnlegend(3,legend_str, 'Location', 'NorthWest', 'boxoff')
            title({char(obj.PCFS_FID),' Uncorrected Cross-Correlations at Different Stage Positions'})
            xlim([1E7,1E13])
            ylim([0.5,2.5])
            set(gca,'FontSize',12)
            set(gca,'XTick',[]);
            
            
            %% plot the autocorrelation functions.
            legend_str=[];
            %figure('rend','painters','pos',[10 10 1300 300])
            subplot(2,1,2)
            for i=1:(length(obj.stage_positions))
                j{i}=semilogx(tau,obj.auto_corr_sum_interferogram(2:end,i),'Color',cm(i,:))
                %legend_str=[legend_str;{num2str(obj.stage_positions(i-4))}];
                hold on
            end
            xlabel('\tau [ps]')
            ylabel('g^{(2)}_{auto}')
            legend([j{1},j{length(obj.stage_positions)}],{num2str(obj.stage_positions(1)),num2str(obj.stage_positions(end))})
            %columnlegend(3,legend_str, 'Location', 'NorthWest', 'boxoff')
            title({char(obj.PCFS_FID),' Auto-Correlations at Different Stage Positions'})
            xlim([1E7,1E13])
            ylim([0.8,2.5])
            set(gca,'FontSize',12)
            
            
            %% plot a heat map of the total interferogram
            figure()
            colormap hsv
            contourf(obj.stage_positions,tau,obj.cross_corr_interferogram(2:end,1:end),1000,'edgecolor','none')
            title({char(obj.PCFS_FID),' Cross-correlation (Uncorrected)PCFS-interferogram'})
            xlabel('Stage Position [mm]')
            ylabel('\tau [ps]')
            ylim([1E7,1E12])
            %zlim([0.5,1.5])
            colorbar
            set(gca, 'CLim', [0.5, 1.5]);
            set(gca,'yscale','log');
            set(gca,'FontSize',12)
            
        end
        function obj=get_SPCFS_contributions(obj,tau_plateau,tau_single,tau_limits,stage_limits,tau_diffusion,n_slice_avg,fit_gaussian)
            
            % tau_plateau forms the lower and upper bounds of tau values (pulses or ps) in [a,b] form that are used to
            % calculate the ensemble interferogram.
            
            % tau_single are the two tau values [a,b] marking high single
            % molecule contrast in the spectral correlation - used to average
            % the single molecule component.
            
            %tau_limits and stage_limits are vectors of the form [lowerlimit,upperlimit]
            %to define the limits for parsing the interferogram. Data outside
            %these limits will be disregarded for the further analysis.
            
            %%takes the auto-correlation of the sum signal and the
            %%cross-correlations, obtained from .get_intensity_correlations and
            %%separates the degree of anti-correlation at each stage position
            %%into the single molecule and ensemble contribution.
            
            %The ensemble contribution is taken to be the plateau-value at times < the dwell time of the emitter
            %g(cross,single)(tau,delta)=(g(cross,delta)-g(cross,ensembe tau->inf,delta))/(g(auto,tau,delta)-1)
            
            %tau_plateu is a vactor containing the lower and upper bound (in
            %pulses) for taking the plateu value
            
            %correct the observed cross correlation for the auto-correlation of
            %the sum signal. creates a matrix that contains the anti-correlation
            % due to emission coherence.
            
            %Also allows to plot slices of the single emitter interferogram to
            %analyze for spectral diffusion.
            
            
            %%
            %easier to work with matrix representation (wo tau and stage positions values) of the
            %correlations -- redefine later.
            cross_corr=obj.cross_corr_interferogram(2:end,:);
            auto_corr=obj.auto_corr_sum_interferogram(2:end,:);
            PCFS_interferogram=obj.PCFS_interferogram(2:end,:);% substract the auto_corr (FCS) from the cross-corr
            
            %get ensemble g(cross) at late tau from the input
            tau=obj.tau;
            index_ensemble=zeros(2,1);
            for i=1:2;
                [dummy,index_ensemble(i)]=min((tau- tau_plateau(i)).^2);
            end
            
            %initialize and populate a ensemble contribution matrix of the
            %interferogram.
            g_ensemble=zeros(size(cross_corr));
            for i=1:length(obj.stage_positions)
                g_ensemble_average=mean(PCFS_interferogram(index_ensemble(1):index_ensemble(2),i));
                g_ensemble(:,i)=ones(length(cross_corr),1)*g_ensemble_average;
            end
            
            %populate the single emitter contribution to the interferogram
            %i.e. substracting the ensemble contributon from the difference
            %between the autocorrelation of the sum-signal and the cross
            %correlation.
            g_single=[];
            g_single=(PCFS_interferogram-g_ensemble);
            
            
            %normalize and average selected single molecule and ensemble contributions.
            index_single=zeros(2,1);
            for i=1:2;
                [dummy,index_single(i)]=min((tau- tau_single(i)).^2);
            end
            
            %normalize the entire single emitter and ensemble PCFS interferogram matrix
            g_single_norm=[];
            for i=1:length(obj.stage_positions);
                g_single_norm(:,i)=g_single(:,i)/min(g_single(:,i));
            end
            
            g_ensemble_norm=g_ensemble(1,:)/min(g_ensemble(1,:));
            g_single_avg=sum(g_single_norm(index_single(1):index_single(2),:))/max(sum(g_single_norm(index_single(1):index_single(2),:)));
            
            
            %%getting the part of the interferogram selected by tau_limits and stage_limit.
            index_tau=zeros(2);
            index_stage=zeros(2);
            for i=1:2;
                [dummy,index_tau(i)]=min((tau- tau_limits(i)).^2);
                [dummy,index_stage(i)]=min((obj.stage_positions- stage_limits(i)).^2);
            end
            g_single_select=g_single(index_tau(1):index_tau(2),index_stage(1):index_stage(2));
            tau_select=obj.tau(index_tau(1):index_tau(2));
            stage_positions_select=obj.stage_positions(index_stage(1):index_stage(2));
            
            %%normalize the selected part of the single emitter interferogram
            g_single_select_norm=[];
            for i=1:length(tau_select);
                g_single_select_norm(i,:)=g_single_select(i,:)/min(g_single_select(i,:));
            end
            
            %% write new observables as object properties.
               
            if ~isprop(obj,'g_single')
                dummy=addprop(obj,'g_single');
            end
            
            if ~isprop(obj,'g_ensemble')
                dummy=addprop(obj,'g_ensemble');
            end
            
            if ~isprop(obj,'g_ensemble_norm')
                dummy=addprop(obj,'g_ensemble_norm');
            end
            
            if ~isprop(obj,'g_single_norm')
                dummy=addprop(obj,'g_single_norm');
            end
            
            if ~isprop(obj,'g_single_avg')
                dummy=addprop(obj,'g_single_avg');
            end
            
            if ~isprop(obj,'g_single_select_norm')
                dummy=addprop(obj,'g_single_select_norm');
            end
            
            if ~isprop(obj,'tau_select')
                dummy=addprop(obj,'tau_select');
            end
            
            if ~isprop(obj,'stage_positions_select')
                dummy=addprop(obj,'stage_positions_select');
            end
            
            
            
            obj.g_ensemble=g_ensemble;
            obj.g_ensemble_norm=g_ensemble_norm;
            
            obj.g_single=g_single;
            obj.g_single_avg=g_single_avg;
            obj.g_single_norm=g_single_norm;
            
            obj.tau_select=tau_select;
            obj.stage_positions_select=stage_positions_select;
            obj.g_single_select_norm=g_single_select_norm;
            
            %% plot the averaged single and ensemble contributions as a function of stage position
            
            figure
            plot(obj.stage_positions,obj.g_single_avg)
            hold on
            plot(obj.stage_positions,obj.g_ensemble_norm(1,:))
            title(strcat(obj.PCFS_FID, ' Normalized PCFS Interferogram of Single Emitter and Ensemble'))
            xlabel('Stage Position [mm]')
            ylabel('Normalized g^{(2)}_{cross}-g^{(2)}_{auto}')
            legend('Single','Ensemble')
            
            
            %% plot a heat map of the selected normalized single emitter interferogram.
            
            figure()
            contourf(obj.stage_positions_select,obj.tau_select,obj.g_single_select_norm,20,'edgecolor','none')
            title(strcat(obj.PCFS_FID,' Single emitter PCFS-interferogram'))
            xlabel('Stage Position [mm]')
            ylabel('\tau [ps or pulses]')
            ylim([str2num(obj.corr_binwidth),1E12])
            %zlim([0.5,1.5])
            colorbar
            set(gca, 'CLim', [0.9, 1]);
            set(gca,'yscale','log');
            colorbar
            
            %% spectral diffusion analysis
            %plotting normalized single molecule interferograms of selected tau
            spectral_diffusion=[];
            figure();
            for i=1:length(tau_diffusion);
                [dummy,index]=min((obj.tau_select-tau_diffusion(i)).^2);
                interferogram=obj.g_single_select_norm(index,:);%mean(obj.g_single_select_norm((index-n_slice_avg):(index+n_slice_avg),:));%/min(mean(obj.g_single_select_norm((index-n_slice_avg):(index+n_slice_avg),:)));
                %spectral_diffusion=[spectral_diffusion;interferogram];
                % disp(size(spectral_diffusion))
                plot(obj.stage_positions_select,interferogram)%spectral_diffusion(i,:));
                hold on
            end
            xlabel('Stage Position [mm]')
            % legend(strread(num2str(tau_select),'%s'))
            ylabel('g^{(2)}_{cross} - g^{(2)}_{auto}')
            title(strcat(obj.PCFS_FID, ' Single emitter interferogram at different \tau'))
            
            %fit the single emitter interferogram with gaussians
            c_param=zeros(176,3);
            if fit_gaussian==true
                for i=1:153
                    disp(i)
                    f=fit(obj.stage_positions(1:end-3)',obj.g_single_norm(i,1:end-3)','gauss1')
                    c_param(i,1)=f.c1;
                    a=confint(f)
                    c_param(i,2:3)=a(:,3)';
                end
                
                if ~isprop(obj,'sigma_spectral_diffusoin')
                    dummy=addprop(obj,'sigma_spectral_diffusoin');
                end
                
                obj.sigma_spectral_diffusoin=c_param/sqrt(2);
            end
            
            %%plot sigma of the single NC interferogram as a function of tau
            figure()
            semilogx(obj.tau(1:176),obj.sigma_spectral_diffusoin(:,1),'linewidth',1,'color','red')
            hold on
            semilogx(obj.tau(1:176),obj.sigma_spectral_diffusoin(:,2),'--','color','black','linewidth',0.3)
            semilogx(obj.tau(1:176),obj.sigma_spectral_diffusoin(:,3),'--','color','black','linewidth',0.3)
            xlim([1E5,1E10])
            ylim([0,5E-3])
            xlabel('\tau [ps]')
            ylabel('Single NC interferogram width \sigma [\mum]')
            title(strcat(obj.PCFS_FID, 'Spectral Diffusion Analysis'))
            set(gca,'FontSize',16)
            
        end
        function obj=center_Gaussian_interferograms(obj,file_out_key,stage_positions,interferogram)
            %%takes the average single and ensemble interferograms, fits them
            %%to a Gaussian and interpolates the data around the center of the
            %%Gaussian to reduce error when taking the Fourier transform to
            %%obtian the spectral correlation.
            n=length(stage_positions);
            delta=stage_positions(2)-stage_positions(1);
            Gaussian_fit=fit(stage_positions',interferogram','Gauss1');
            
            %interpolation vector centered around the center of the Gaussian.
            xq=linspace((Gaussian_fit.b1-n/2*delta),(Gaussian_fit.b1+n/2*delta),n+1);
            interpolated_interferogram=interp1(stage_positions, interferogram,xq);
            %replace NaN with zeros in the interferogram.
            interpolated_interferogram(isnan(interpolated_interferogram)) = 0;
            
            if ~isprop(obj,(file_out_key))
                dummy=addprop(obj,(file_out_key))
            end
            obj.(file_out_key)=struct('fit',Gaussian_fit,'interp_stage_pos',xq,'interp_interferogram',interpolated_interferogram);
            
        end
        function obj=get_spectral_corr(obj,file_out_key,stage_positions,interferogram, white_fringe_pos)
                %%calculates, returns and plots the spectral correlation of an
                %%interterferogram parsed as two vectors containing the stage_positions
                %%(not path length differences!), the corresponding interferogram
                %%values and the white-fringe position.
                %file_key_out corresponds to the identifier for the spectral
                %correlation vector property of the PCFS instance.
                
                eV2cm=8065.54429;
                cm2eV=1/eV2cm;
                
                N=length(stage_positions);
                path_length_difference=2*(stage_positions-white_fringe_pos)*0.1;%% NOTE: This is where we convert to path length difference space in cm.
                delta=(max(path_length_difference)-min(path_length_difference))/N;
                
                %get reciprocal space (wavenumbers).
                increment=1/delta;
                zeta_eV=linspace(-0.5*increment,0.5*increment,N)*cm2eV*1000; %converted to meV
                
                %%take the FFT of the interferogram to get the spectral correlation.All
                %%that shifting is to shift the zero frequency component to the middle
                %%of the FFT vector. We take the real part of the FFT because the
                %%interferogram is by definition entirely symmetric.
                spectral_corr = real(fftshift(fft(ifftshift(interferogram),N)));
                normalized_spectral_correlatoin=spectral_corr;
                
              %  figure()
            %    plot(zeta_eV,normalized_spectral_correlatoin)
             %   title('Spectral Correlation')
             %   xlabel('\zeta [meV]')
             %   ylabel('intensity [a.u.]')
                
                
                %%perform Gaussian fit for the respective spectral correlation
                f=fit(zeta_eV',normalized_spectral_correlatoin','gauss1')
                
                if ~isprop(obj,(file_out_key))
                    dummy=addprop(obj,(file_out_key))
                end
                obj.(file_out_key)=struct('Gaussian_fit',f,'spectral_corr',normalized_spectral_correlatoin,'zeta',zeta_eV);
                
            end
        function obj=plot_single_vs_ensemble_spectral_corr(obj,white_fringe_pos);
                %% plots the single and ensemble spectral correation together with Gaussian Fits.
                % The single and ensemble interferograms need to be saved as
                % object properties like obj.g_single_avg and
                % obj.g_ensemble_norm by calling .get_SPCFS_contributions.
                
                %call get_spectral_corr for the single and ensemble
                %interferogram.
                obj.center_Gaussian_interferograms('ensemble',obj.stage_positions,obj.g_ensemble_norm);
                obj.center_Gaussian_interferograms('single',obj.stage_positions,obj.g_single_avg);

                
                obj.get_spectral_corr('single',obj.single.interp_stage_pos,obj.single.interp_interferogram, white_fringe_pos);
                obj.get_spectral_corr('ensemble',obj.ensemble.interp_stage_pos,obj.ensemble.interp_interferogram, white_fringe_pos);
                
                %normalize the spectral correlations.
                spectral_corr_single=obj.single.spectral_corr/max(obj.single.spectral_corr);
                spectral_corr_ensemble=obj.ensemble.spectral_corr/max(obj.ensemble.spectral_corr);
                Gaussian_fit_single=fit(obj.single.zeta',spectral_corr_single','Gauss1');
                Gaussian_fit_ensemble=fit(obj.ensemble.zeta',spectral_corr_ensemble','Gauss1');
                
                
                if ~isprop(obj,'Norm_single_spectral_corr')
                    dummy=addprop(obj,'Norm_single_spectral_corr')
                end
                
                
                if ~isprop(obj,'Norm_ensemble_spectral_corr')
                    dummy=addprop(obj,'Norm_ensemble_spectral_corr')
                end
                
                if ~isprop(obj,'Gaussian_fit_single')
                    dummy=addprop(obj,'Gaussian_fit_single')
                end
                
                if ~isprop(obj,'Gaussian_fit_ensemble')
                    dummy=addprop(obj,'Gaussian_fit_ensemble')
                end
                
                obj.Norm_ensemble_spectral_corr=spectral_corr_ensemble;
                obj.Norm_single_spectral_corr=spectral_corr_single;
                
                obj.Gaussian_fit_ensemble=Gaussian_fit_ensemble;
                obj.Gaussian_fit_single=Gaussian_fit_single;
                
                single_conf_int=confint(obj.Gaussian_fit_single);
                ensemble_conf_int=confint(obj.Gaussian_fit_ensemble);
                
                
                %%plottting the comparison of single and ensemble spectral
                %%corr.
                plot(obj.single.zeta,spectral_corr_single,'linewidth',1.5)
                hold on
                plot(obj.ensemble.zeta,spectral_corr_ensemble,'linewidth',1.5)
                plot(linspace(-1000,1000,5000),Gaussian_fit_single(linspace(-1000,1000,5000)),'--');
                plot(linspace(-1000,1000,5000),Gaussian_fit_ensemble(linspace(-1000,1000,5000)),'--');
                legend('Average Single NC','Ensemble')
                ylim([-0.1,1.1])
                xlim([-400,400])
                xlabel('\zeta [meV]')
                ylabel('Norm. Spectral Corr. [a.u.]')
                %title(strcat(obj.PCFS_FID,'   Average Single vs. Ensemble Spectral Corr.'))
                title('Single vs. Ensemble Spectral Correlation')
                text(-380,0.8,{strcat('Single FWHM:',num2str(obj.Gaussian_fit_single.c1/sqrt(2)*2.3548,'%0.2f'),' meV'),strcat(' (',num2str(single_conf_int(1,3)/sqrt(2)*2.3548,'%0.2f'),',',num2str(single_conf_int(2,3)/sqrt(2)*2.3548,'%0.2f'),')')},'fontsize',13)
                text(-380,0.9,{strcat('Ensemble FWHM:',num2str(obj.Gaussian_fit_ensemble.c1/sqrt(2)*2.3548,'%0.2f'),' meV'),strcat(' (',num2str(ensemble_conf_int(1,3)/sqrt(2)*2.3548,'%0.2f'),',',num2str(ensemble_conf_int(2,3)/sqrt(2)*2.3548,'%0.2f'),')')},'fontsize',13)
                set(gca,'FontSize',16)
                
                
        end
            
        %% additional analysis functions
        function obj=plot_Fourier_Spectrum(obj)
                %%Plots the Fourier Interferogram and Spectrum created with the
                %Labview .scan functionality (file extension: .intf).
                files=dir;
                obj.files=files;
                for i=3:numel(obj.files)
                    [path,f_name,ext]=fileparts(obj.files(i).name);
                    if strcmp(ext,'.intf')== true;
                        %read data
                        filename = strcat(f_name,'.intf');
                        delimiter = '\t';
                        startRow = 8;
                        formatSpec = '%f%f%f%f%f%f%[^\n\r]';
                        fileID = fopen(filename,'r');
                        textscan(fileID, '%[^\n\r]', startRow-1, 'WhiteSpace', '', 'ReturnOnError', false);
                        dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'ReturnOnError', false);
                        fclose(fileID);
                        data = [dataArray{1:end-1}];
                        clearvars filename delimiter startRow formatSpec fileID dataArray ans;
                        if ~isprop(obj,strcat(f_name,'_intf'));
                            dummy=addprop(obj,strcat(f_name,'_intf'));
                        end
                        obj.(strcat(f_name,'_intf'))=data;
                        %make a figure and plot the interferogram and the
                        %Fourier transformed spectrum.
                        path_length_difference=2*(data(:,1))*0.1;%% NOTE: This is where we convert to path length difference space in cm.
                        
                        figure()
                        plot(path_length_difference,data(:,3),path_length_difference,data(:,4))%converted in cm.
                        hold on
                        xlabel('Path Length Difference + Const. [mm]')
                        ylabel('Intensity [cps]')
                        
                        %calculate Fourier Transform.
                        wavepacket=(obj.(strcat(f_name,'_intf'))(:,3)-obj.(strcat(f_name,'_intf'))(:,4))-mean((obj.(strcat(f_name,'_intf'))(:,3)-obj.(strcat(f_name,'_intf'))(:,4)));
                        
                        figure()
                        plot(path_length_difference,wavepacket)
                        
                        
                        %
                        % %%getting the FFT of the wavepacket.
                        %
                        eV2cm=8065.54429;%cm^-1/eV
                        cm2eV=1/eV2cm;%eV/cm^-1
                        
                        N=length(path_length_difference);
                        %
                        %
                        steps=(max(path_length_difference)-min(path_length_difference))/N;
                        max_energy=1/steps;
                        energy_eV=linspace(0,max_energy,N)*cm2eV;
                        
                        %
                        % %%take the FFT of the interferogram to get the spectral correlation.All
                        % %%that shifting is to shift the zero frequency component to the middle
                        % %%of the FFT vector. We take the real part of the FFT because the
                        % %%interferogram is by definition entirely symmetric.
                        %
                        spectrum = abs(fft(wavepacket))/max(abs(fft(wavepacket)));
                        figure()
                        plot(energy_eV(1:N/2),spectrum(1:N/2))
                        xlim([0,3])
                        
                        if ~isprop(obj,strcat(f_name,'_Fourier_spectrum'));
                            dummy=addprop(obj,strcat(f_name,'_Fourier_spectrum'));
                        end
                        obj.(strcat(f_name,'_Fourier_spectrum'))=struct('energy',energy_eV(1:N/2),'spectrum',spectrum(1:N/2));
                        
                        
                    end
                end
            end
        function obj=validation_ensemble_spectral_correlation(obj,camera_spectrum)
                %calculates the auto-correlation of an ensemble spectrum taken
                %with the camera and Fourier transforms it to obtein the PCFS
                %interferogram of the camera data. Camera spectrum should be a
                %512x2 matrix containing (:,1) the spectrum and (:,2) the
                %calibration file in eVas aquired with PI's Lightfield.
                
                %returns the energy difference variable zeta in eV and cm,
                %the spectral correlation, and the interferogram.
                
                spectral_corr_cam=xcorr(camera_spectrum(:,1)');
                
                eV=camera_spectrum(:,2);
                nm=1240./eV;
                cm=nm/1E7;
                
                N=(length(spectral_corr_cam)+1)/2;%just the lenght of the original spectrum.
                delta_zeta=camera_spectrum(1,2)-camera_spectrum(2,2);%delta in energy space is equal to delta in zeta space.
                disp(delta_zeta)
                
                zeta=linspace((-N-1),(N-1),(2*N-1))*delta_zeta;
                delta_per_cm=8065.54429*delta_zeta;
                zeta_percm=linspace((-N-1),(N-1),(2*N-1))*delta_per_cm;
                
                figure()
                plot(zeta_percm, spectral_corr_cam)
                title('spectral correlation from camera data')
                xlabel('wavenumbers [cm^{-1}]')
                ylabel('intensity [a.u.]')
                
                figure()
                plot(zeta, spectral_corr_cam)
                title('spectral correlation from camera data')
                xlabel('eV')
                ylabel('intensity [a.u.]')
                
                %get reciprocal space.
                increment=1/delta_per_cm;
                zeta_cm=linspace(-increment/2,increment/2,(2*N-1))*1E4;%converted to \mu m
                
                %%take the FFT of the spectral_correlation to get interferogram of the
                %%camera data.
                interferogram_cam = real(fftshift(fft(ifftshift(spectral_corr_cam),2*N-1))/sqrt(N+1));
                
                figure()
                plot(zeta_cm,interferogram_cam)
                title('PCFS ensemble interferogram from camera data')
                xlabel('\delta [\mu m]')
                ylabel('intensity [a.u.]')
                
                
                %packing the results up as properties.
                if ~isprop(obj,'zeta')
                    dummy=addprop(obj,'zeta')
                end
                if ~isprop(obj,'zeta_cm')
                    dummy=addprop(obj,'zeta_cm')
                end
                if ~isprop(obj,'spectral_corr_cam')
                    dummy=addprop(obj,'spectral_corr_cam')
                end
                if ~isprop(obj,'interferogram_cam')
                    dummy=addprop(obj,'interferogram_cam')
                end
                
                obj.zeta=zeta;
                obj.zeta_cm=zeta_cm;
                obj.spectral_corr_cam=spectral_corr_cam;
                obj.interferogram_cam=interferogram_cam;
                
            end
        function obj=fit_FCS_traces(obj)
                %% fits all auto-correlations  of the sum signal swith a simple FCS equation  and creates an array with the fit_
                %  curves, parameters, and residuals.
                
                %call photons.fit_auto_corr_FCS_trace(obj,file_key_in,p0,afterpulsing_time)
                %for each photonstream.
                files=dir;
                obj.files=files;
                
                FCS_curves=zeros(size(obj.auto_corr_sum_interferogram(2:end,5:end)));
                P=zeros(4,length(obj.stage_positions));
                rel_errors=zeros(size(obj.auto_corr_sum_interferogram(2:end,5:end)));
                
                for i=3:numel(obj.files) %exclusing . and ..
                    [path,f_name,ext]=fileparts(obj.files(i).name);
                    if strcmp(ext,'.intensity_corr')== true
                        if strcmp(f_name(1:3),'sum')==true
                            %we are looking at an auto-correlation function of
                            %the sum signal.
                            obj.(f_name(12:end)).fit_auto_corr_FCS_trace((f_name));
                            
                            %extract the correlation number from the filename.
                            k=-1;
                            r=[];
                            while true
                                k=k+1;
                                if isstrprop(f_name(end-k),'digit')==true;
                                    r=[f_name(end-k),r];
                                else
                                    correlation_number=str2num(r);
                                    r=[];
                                    break
                                end
                            end
                            FCS_curves(:,correlation_number+1)=obj.(f_name(12:end)).FCS_fit.fit_curve;
                            P(:,correlation_number+1)=obj.(f_name(12:end)).FCS_fit.P';
                            residuals(1:length(obj.(f_name(12:end)).FCS_fit.residuals),correlation_number+1)=obj.(f_name(12:end)).FCS_fit.residuals';
                            tau_for_fit=obj.(f_name(12:end)).FCS_fit.tau_for_fit;
                        end
                    end
                end
                
                if ~isprop(obj,'FCS_analysis')
                    dummy=addprop(obj,'FCS_analysis');
                end
                
                rel_errors( ~any(rel_errors,2), : ) = [];  %delete empty rows
                wrapped_data=struct('FCS_curves',FCS_curves,'residuals',residuals,'Params',P,'tau_for_fit',tau_for_fit);
                obj.FCS_analysis=wrapped_data;
                
                figure()
                subplot(5,1,1)
                plot(obj.stage_positions,sum((obj.FCS_analysis.residuals).^2,1))
                ylabel('SRS')
                xlabel('Stage position')
                title('Sum of residuals squared of FCS fit to g^{(2)}_{auto} of sum-signal ')
                
                subplot(5,1,2)
                plot(obj.stage_positions,obj.FCS_analysis.Params(2,:))
                ylabel('Dwell Time')
                xlabel('Stage position')
                title('Dwell Time of FCS fit to g^{(2)}_{auto} of sum-signal ')
                
                subplot(5,1,3)
                plot(obj.stage_positions,obj.FCS_analysis.Params(4,:))
                ylabel('g^{(2)}(\tau=0)')
                xlabel('Stage position')
                title('g^{(2)}(\tau=0) of FCS fit to g^{(2)}_{auto} of sum-signal ')
                
                subplot(5,1,4)
                plot(obj.stage_positions,obj.FCS_analysis.Params(3,:))
                ylabel('\alpha')
                xlabel('Stage position')
                title('Aspect ratio of convocal volume ')
                
                subplot(5,1,5)
                plot(obj.stage_positions,obj.FCS_analysis.Params(1,:))
                ylabel('g^{(2)}_{(\tau= long)}')
                xlabel('Stage position')
                title('g^{(2)}_{(\tau= long)}')
                
                
        end
        function obj=get_Fourier_spectrum_from_stream(obj,fid)
             
                Fourier=[];
                int=zeros(numel(obj.stage_positions),1);
                pos=obj.stage_positions;
                 
                files=dir;
                obj.files=files;
                
                lstr=length(fid);
                
                for i=3:numel(obj.files) %exclusing . and ..
                    [path,f_name,ext]=fileparts(obj.files(i).name);
                    if strcmp(ext,'.photons')== true
                        if strcmp(f_name(1:3),'sum')==false %looking to get cross-corrs
                             if strcmp(f_name(1:lstr),fid)==true
                            %if isprop(obj,f_name)==true
                                k=-1;
                                r=[];
                                while true
                                    k=k+1;
                                    if isstrprop(f_name(end-k),'digit')==true;
                                        r=[f_name(end-k),r];
                                    else
                                        correlation_number=str2num(r);
                                        r=[];
                                        break
                                    end
                                end
                                
                                int_id=strcat(fid,'_intensity');
                                intensity=obj.(f_name((lstr+1):end)).(int_id).trace;
                                
                                NA_elements=[];
                                
                                if isempty(intensity)==true
                                    NA_elements=[NA_elements,(correlation_number+1)]
                                else
                                    Fourier(correlation_number+1)=(sum(intensity(:,1))-sum(intensity(:,2)))./(sum(intensity(:,1))+sum(intensity(:,2)));
                                    int(correlation_number+1)=(sum(intensity(:,1))+sum(intensity(:,2)));
                                end
                             end
                        end
                    end
                end
            
             % deleting empty elements. 
             stage_pos=obj.stage_positions;
             stage_pos(NA_elements)=[];
             Fourier(NA_elements)=[];
             int(NA_elements)=[];
                
             Fourier_ID=strcat(fid,'_Fourier');
             if ~isprop(obj,Fourier_ID)
                    dummy=addprop(obj,Fourier_ID)
             end
             obj.(Fourier_ID)=struct('Fourier',Fourier,'stage_pos',stage_pos,'intensity',int);
        end
            
        %% depreciated functions
        function obj=bin_2_int_all(obj)
                %batch convert all binary photon files to .photons_int to be
                %useable with Tom's correlation code (photon_intensity_correlate)
                files=dir;
                obj.files=files;
                for i=3:numel(obj.files)
                    [path,f_name,ext]=fileparts(obj.files(i).name);
                    if strcmp(ext,'.photons')== true
                        if strcmp(f_name(1:3),'sum')==0;
                            obj.(f_name).bin_2_int(obj.(f_name).fname,obj.(f_name).fname)
                        end
                        if strcmp(f_name(1:3),'sum')==1;
                            obj.(f_name(12:end)).bin_2_int(f_name,f_name)
                        end
                    end
                end
                files=dir;
                obj.files=files;
                
            end
        function obj=intensity_correlate_pipeline_mac(obj,mode,bin_width,time_scale)
                %basically a systems call to pipe the picoquant output into
                %photon_intensity_correlate. This avoids creating the .photons and .photons_int file
                %to save storage space, but does not allow to perform further
                %photon stream analysis.
                
                
                n_channels=3;
                for i=3:numel(obj.files)
                    [path,f_name,ext]=fileparts(obj.files(i).name);
                    %loop over all .stream files in obj.files and pipe picoquant output
                    %into photon_intensity correlate.
                    
                    
                    
                    if strcmp(ext,'.stream')== true
                        
                        command=['picoquant --file-in ', [path,f_name,'.stream'], '| photons --copy-to-channel 2 --mode ', mode, '| photon_intensity_correlate --file-out ', [path,f_name,'.intensity_corr'], ' -m ', mode, ' -g 2 -w ', bin_width,' --time-scale ',time_scale,' -c ', num2str(n_channels)]
                        disp(command)
                        system(command)
                        
                        %get the correlation function of the sum signal
                        
                        
                        
                        
                        
                    end
                end
            end
            
    end
end
    
