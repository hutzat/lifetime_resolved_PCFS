%% MATLAB Implementation for analysis of photon stream data from the Picoqunt GmbH
% Hydraharp.

%% photons.m V5.0 @ HENDRIK UTZAT, KATIE SHULENBERGER, BORIS SPOKOYNY, TIMOTHY SINCLAIR, 2017
% 10/29/2017

% This is an object-oritented implementation of a comprehensive
% analysis library for time-tagged single photon arrival time data as
% recorded with the Picoquant GmbH Hydraharp picosecond photon timing
% card. 

% After reading the metadata with the constructor methods (.photons), the function
% .get_photon_records() parses the raw hoton-arrival time data recorded with an 
% absolute experimental clock (T2 data, .ht2) or a relative
% experimental clock (T3 data, .ht3) into the true photon-arrival time
% records. These photon-records are written to a .photons
% file in uint64 representation. 

% Subseqent analysis of the photon-stream reads the true photon arrival 
% events from the .photons files. Various photon-sorting functions exist. 
% As such the photon-stream can be sorted according to intensity, 
% detection channel, arrival after the sync pulse etc. Detector after-pulsing
% can be removed. Each call to a sorting/cleaning function creates a new 
% .photons file that contains only the photon-records meeting the specified 
% requirements. 

% The .photons files can be read by various analysis functions to obtain:
% intensity traces, fluorescence lifetimes, linearly time-spaced
% coincidence counts (g(2),g(3)...), logarithmically time-spaced
% correlations functions, or photon-number-resolved lifetimes. A function
% to fit FCS decay of emitters diffusing in solution exists as well. 



classdef photons<dynamicprops
    %%%% creating sub-class of dynamicprops which allows
    %to dynamically add propterties to an instance of that class
    
    %Initialize some properties of photons class.
    properties
        file_path=[];pathstr=[];fname=[];fext=[];                   %identifier for the photon-stream file.
        header_size=[];                                             %size of header in bytes
        header=[];                                                  %all header info from the .ht2 or .ht3 files.
        in_ch0=[];in_ch1=[];in_ch2=[];in_ch3=[];                    %average channel counts
        buffer_size=1E6;                                            %default buffer size for reading in photon records in binary
    end
    
    
    methods
        %% constructor method
        function obj=photons(file_path,buffer_size)
            %constructor mehtod called when creating an instance of the
            %photon class. It reads the header info from the input photon-stream file.
            
            %file_path: path to the photon-stream file
            %buffer_size: number of photon records to read to memory at
            %once. default is 1E6
            
            %set buffer_size for reading in photons to buffer_size if
            %sepcified.
            if nargin<2
                obj.buffer_size=buffer_size;
            end
            
            %extract photon-stream file path info.
            [pathstr,name,ext] = fileparts(file_path);
            obj.fname=name;
            obj.pathstr=pathstr;
            obj.fext=ext;
            obj.file_path=file_path;
            
            %read in the header file info.
            obj.get_photon_stream_header();
            
            disp('-------------------------------')
            disp('Instance of photon class created')
            disp('-------------------------------')
            
        end
        
        %% parse the photon stream and meta data.
        function obj=get_photon_stream_header(obj)
            % read the ASCI and binary header of the photon stream file
            % code adapted from Picoquant GmBH, Germany.
            
            fid=fopen(obj.file_path);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            % ASCII file header
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            Ident = char(fread(fid, 16, 'char'));
            fprintf(1,'               Ident: %s\n', Ident);
            
            FormatVersion = deblank(char(fread(fid, 6, 'char')'));
            fprintf(1,'      Format version: %s\n', FormatVersion);
            
            if not(strcmp(FormatVersion,'2.0'))
                fprintf(1,'\n\n      Warning: This program is for version 2.0 only. Aborted.');
                STOP;
            end;
            
            CreatorName = char(fread(fid, 18, 'char'));
            fprintf(1,'        Creator name: %s\n', CreatorName);
            
            CreatorVersion = char(fread(fid, 12, 'char'));
            fprintf(1,'     Creator version: %s\n', CreatorVersion);
            
            FileTime = char(fread(fid, 18, 'char'));
            fprintf(1,'    Time of creation: %s\n', FileTime);
            
            CRLF = char(fread(fid, 2, 'char'));
            
            Comment = char(fread(fid, 256, 'char'));
            fprintf(1,'             Comment: %s\n', Comment);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            % Binary file header
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % The binary file header information is indentical to that in HHD files.
            % Note that some items are not meaningful in the time tagging modes
            % therefore we do not output them.
            
            NumberOfCurves = fread(fid, 1, 'int32');
            
            BitsPerRecord = fread(fid, 1, 'int32');
            fprintf(1,'       Bits / Record: %d\n', BitsPerRecord);
            
            ActiveCurve = fread(fid, 1, 'int32');
            
            MeasurementMode = fread(fid, 1, 'int32');
            fprintf(1,'    Measurement Mode: %d\n', MeasurementMode);
            
            SubMode = fread(fid, 1, 'int32');
            fprintf(1,'            Sub-Mode: %d\n', SubMode);
            
            Binning = fread(fid, 1, 'int32');
            fprintf(1,'             Binning: %d\n', Binning);
            
            Resolution = fread(fid, 1, 'double');
            fprintf(1,'          Resolution: %f ps\n', Resolution);
            
            Offset = fread(fid, 1, 'int32');
            fprintf(1,'              Offset: %d\n', Offset);
            
            Tacq = fread(fid, 1, 'int32');
            fprintf(1,'    Acquisition Time: %d ms \n', Tacq);
            
            StopAt = fread(fid, 1, 'uint32');
            StopOnOvfl = fread(fid, 1, 'int32');
            Restart = fread(fid, 1, 'int32');
            DispLinLog = fread(fid, 1, 'int32');
            DispTimeAxisFrom = fread(fid, 1, 'int32');
            DispTimeAxisTo = fread(fid, 1, 'int32');
            DispCountAxisFrom = fread(fid, 1, 'int32');
            DispCountAxisTo = fread(fid, 1, 'int32');
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            for i = 1:8
                DispCurveMapTo(i) = fread(fid, 1, 'int32');
                DispCurveShow(i) = fread(fid, 1, 'int32');
            end;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            for i = 1:3
                ParamStart(i) = fread(fid, 1, 'float');
                ParamStep(i) = fread(fid, 1, 'float');
                ParamEnd(i) = fread(fid, 1, 'float');
            end;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            RepeatMode = fread(fid, 1, 'int32');
            RepeatsPerCurve = fread(fid, 1, 'int32');
            Repaobjime = fread(fid, 1, 'int32');
            RepeatWaiobjime = fread(fid, 1, 'int32');
            ScriptName = char(fread(fid, 20, 'char'));
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            %          Hardware information header
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            fprintf(1,'-------------------------------------\n');
            
            HardwareIdent = char(fread(fid, 16, 'char'));
            fprintf(1,' Hardware Identifier: %s\n', HardwareIdent);
            
            HardwarePartNo = char(fread(fid, 8, 'char'));
            fprintf(1,'Hardware Part Number: %s\n', HardwarePartNo);
            
            HardwareSerial = fread(fid, 1, 'int32');
            fprintf(1,'    HW Serial Number: %d\n', HardwareSerial);
            
            nModulesPresent = fread(fid, 1, 'int32');
            fprintf(1,'     Modules present: %d\n', nModulesPresent);
            
            for i=1:10
                ModelCode(i) = fread(fid, 1, 'int32');
                VersionCode(i) = fread(fid, 1, 'int32');
            end;
            for i=1:nModulesPresent
                fprintf(1,'      ModuleInfo[%02d]: %08x %08x\n', i-1, ModelCode(i), VersionCode(i));
            end;
            
            BaseResolution = fread(fid, 1, 'double');
            fprintf(1,'      BaseResolution: %f\n', BaseResolution);
            
            InputsEnabled = fread(fid, 1, 'ubit64');
            fprintf(1,'      Inputs Enabled: %x\n', InputsEnabled); %actually a bitfield
            
            InpChansPresent  = fread(fid, 1, 'int32');
            fprintf(1,' Input Chan. Present: %d\n', InpChansPresent);
            
            RefClockSource  = fread(fid, 1, 'int32');
            fprintf(1,'      RefClockSource: %d\n', RefClockSource);
            
            ExtDevices  = fread(fid, 1, 'int32');
            fprintf(1,'    External Devices: %x\n', ExtDevices); %actually a bitfield
            
            MarkerSeobjings  = fread(fid, 1, 'int32');
            fprintf(1,'     Marker Seobjings: %x\n', MarkerSeobjings); %actually a bitfield
            
            SyncDivider = fread(fid, 1, 'int32');
            fprintf(1,'        Sync divider: %d \n', SyncDivider);
            
            SyncCFDLevel = fread(fid, 1, 'int32');
            fprintf(1,'      Sync CFD Level: %d mV\n', SyncCFDLevel);
            
            SyncCFDZeroCross = fread(fid, 1, 'int32');
            fprintf(1,'  Sync CFD ZeroCross: %d mV\n', SyncCFDZeroCross);
            
            SyncOffset = fread(fid, 1, 'int32');
            fprintf(1,'         Sync Offset: %d\n', SyncOffset);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            %          Channels' information header
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            for i=1:InpChansPresent
                InputModuleIndex(i) = fread(fid, 1, 'int32');
                InputCFDLevel(i) = fread(fid, 1, 'int32');
                InputCFDZeroCross(i) = fread(fid, 1, 'int32');
                InputOffset(i) = fread(fid, 1, 'int32');
                
                fprintf(1,'\n-------------------------------------\n');
                fprintf(1,'Input Channel No. %d\n', i-1);
                fprintf(1,'-------------------------------------\n');
                fprintf(1,'  Input Module Index: %d\n', InputModuleIndex(i));
                fprintf(1,'     Input CFD Level: %d mV\n', InputCFDLevel(i));
                fprintf(1,' Input CFD ZeroCross: %d mV\n', InputCFDZeroCross(i));
                fprintf(1,'        Input Offset: %d\n', InputOffset(i));
            end;
            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            %                Time tagging mode specific header
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            fprintf(1,'\n-------------------------------------\n');
            for i=1:InpChansPresent
                InputRate(i) = fread(fid, 1, 'int32');
                fprintf(1,'     Input Rate [%02d]: %d\n', i-1, InputRate(i));
            end;
            
            fprintf(1,'\-------------------------------------\n');
            
            SyncRate = fread(fid, 1, 'int32');
            fprintf(1,'           Sync Rate: %d Hz\n', SyncRate);
            
            StopAfter = fread(fid, 1, 'int32');
            fprintf(1,'          Stop After: %d ms \n', StopAfter);
            
            StopReason = fread(fid, 1, 'int32');
            fprintf(1,'         Stop Reason: %d\n', StopReason);
            
            ImgHdrSize = fread(fid, 1, 'int32');
            fprintf(1,' Imaging Header Size: %d bytes\n', ImgHdrSize);
            
            nRecords = fread(fid, 1, 'uint64');
            fprintf(1,'   Number of Records: %d\n', nRecords);
            
            % Special header for imaging. How many of the following ImgHdr array elements
            % are actually present in the file is indicated by ImgHdrSize above.
            % Storage must be allocated dynamically if ImgHdrSize other than 0 is found.
            
            ImgHdr = fread(fid, ImgHdrSize, 'int32');  % You have to properly interpret ImgHdr if you want to generate an image
            
            % The header section end after ImgHdr. Following in the file are only event records.
            % How many of them actually are in the file is indicated by nRecords in above.
            
            obj.header_size=ftell(fid) %header size is current byte read in .ht3 file
            fclose(fid)
            
            %wrapping the header info into a struct array for clarity.
            
            header=struct();
            header.('Ident')=Ident;
            header.('FormatVersion')=FormatVersion;
            header.('CreatorVersion')=CreatorVersion;
            header.('Comment')=Comment;
            header.('BitsPerRecord')=BitsPerRecord;
            header.('FileTime')=FileTime;
            header.('CRLF')=CRLF;
            header.('NumberOfCurves')=NumberOfCurves;
            header.('MeasurementMode')=MeasurementMode;
            header.('SubMode')=SubMode;
            header.('Binning')=Binning;
            header.('Resolution')=Resolution;
            header.('Offset')=Offset;
            header.('Tacq')=Tacq;
            header.('StopAt')=StopAt;
            header.('StopOnOvfl')=StopOnOvfl;
            header.('Restart')=Restart;
            header.('DispLinLog')=DispLinLog;
            header.('DispTimeAxisFrom')=DispTimeAxisFrom;
            header.('DispTimeAxisTo')=DispTimeAxisTo;
            header.('DispCountAxisFrom')=DispCountAxisFrom;
            header.('HardwareIdent')=HardwareIdent;
            header.('HardwarePartNo')=HardwarePartNo;
            header.('HardwareSerial')=HardwareSerial;
            header.('nModulesPresent')=nModulesPresent;
            header.('BaseResolution')=BaseResolution;
            header.('InputsEnabled')=InputsEnabled;
            header.('InpChansPresent')=InpChansPresent;
            header.('ExtDevices')=ExtDevices;
            header.('RefClockSource')=RefClockSource;
            header.('SyncDivider')=SyncDivider;
            header.('SyncDivider')=SyncDivider;
            header.('SyncCFDLevel')=SyncCFDLevel;
            header.('SyncCFDZeroCross')=SyncCFDZeroCross;
            header.('SyncOffset')=SyncOffset;
            header.('SyncDivider')=SyncDivider;
            header.('SyncDivider')=SyncDivider;
            header.('SyncDivider')=SyncDivider;
            header.('SyncRate')=SyncRate;
            header.('nRecords')=nRecords;
           
            
            
            %Channels information header.
            for i=1:InpChansPresent
                header.(strcat('channel_',num2str(i)))=struct('InputModuleIndex',InputModuleIndex(i),'InputCFDLevel',InputCFDLevel(i),'InputCFDZeroCross',InputCFDZeroCross(i),'InputOffset',InputOffset(i));
            end
            
            obj.header=header;
            
            
        end
        function get_photon_records(obj)
            %reads the photon stream records and writes 'true' photon
            %records to a binary file. See comment below:
            
            %COMMENT: It is important to understand the difference between
            %'records' and 'photons' according to the convention used by
            %Picoquat GmBH
            %Each record in the .ht3 or .ht2 file has 32 bit which encode a
            %set of paramters:
            
            % .ht2:  channel,absolute time, special character
            % .ht3   channel,number_of_sync_pulse, time after the
            % syncpulse, specical character
            
            %The special cahracter is important to distinguish whether a record is a true
            %photon. If the special bit is 0, the 32 bit record
            %corresponds to a true photon arrival event. If the record is
            %1, the record is a special overflow record and not a true photon arrival event
            
            %Use this overflow (OFL) character, we can encode time (or
            %syncpulses) to infinity, although we only have a finite
            %bit-depth. Specifically, to recover the true arrival time
            %of a photon, we add the number of overflows*the
            %resolution*2^(number of bits) to the actual recorded arrival
            %time.
            
            % In this function, we parse the binary photon-stream data,
            % count the overflow characters, and write a binary file
            % containing the true photon arrival times.
            
            buffer_size=obj.buffer_size;
            header_size=obj.header_size;
            
            outfile = strcat(obj.pathstr,obj.fname,'.photons'); %create an output file that contains the photon records.
            fpout = fopen(outfile,'W');
            fid=fopen(obj.file_path);                           %open the binary photon stream file
            fseek(fid,header_size,'bof');                       %skip over the header to the photon data
            
            
            switch obj.header.MeasurementMode
                
                % reading in t3 data
                case 3
                    
                    syncperiod = 1E9/obj.header.SyncRate;                     %syncperiod in nanoseconds
                    OverflowCorrection = 0;
                    T3WRAPAROUND=1024;                                  %if overflow occured, the true n_sync is n_sync+1024
                    true_photons=zeros(3,length(buffer_size));          %initialize an array to store true_photon records.
                    
                    while 1 %while true
                        
                        batch=fread(fid,buffer_size,'ubit32');          %reading in a multiple of 32 bit registers
                        lbatch=length(batch);
                        
                        k=0;                                            %true photon counting variable
                        for i=1:lbatch;                                 %looping over all records in batch
                            
                            %read and decode the 32 bit register of the ith record
                            nsync = bitand(batch(i),1023);                  %the lowest 10 bits of the ith photon
                            dtime = bitand(bitshift(batch(i),-10),32767);   %the next 15 bits
                            channel = bitand(bitshift(batch(i),-25),63);    %the next 6 bits:%0-4
                            special = bitand(bitshift(batch(i),-31),1);     %the last bit:% MSB - for overflow handling
                            
                            if special == 0                                 %this means a true 'photon' arrival event.
                                true_nSync = OverflowCorrection + nsync;
                                %one nsync time unit equals to "syncperiod" which can be calculated from "SyncRate"
                                time =dtime*obj.header.Resolution;
                                k=k+1;                                      %counting the real photons that we see.
                                true_photons(:,k)=[channel;true_nSync;time];%writing the true photon to an array.
                            else                                            %this means we have a special record; the 'record' is not a 'photon'
                                if channel == 63                            %overflow of nsync occured
                                    if(nsync==0)                            %if nsync is zero it is an old style single oferflow
                                        OverflowCorrection = OverflowCorrection + T3WRAPAROUND;
                                    else                                    %otherwise nsync indicates the number of overflows - THIS IS NEW IN FORMAT V2.0
                                        OverflowCorrection = OverflowCorrection + T3WRAPAROUND*nsync;
                                    end;
                                end;
                            end;
                        end;
                        fwrite(fpout,true_photons(:,1:k),'uint64');         %writing the true photons to the output file in binary.
                        
                        %break the while loop when we have reached the end of the
                        %.ht3 file.
                        if lbatch <buffer_size;
                            break
                        end
                        
                    end
                    fclose(fid);
                    fclose(fpout);
                    
                    %% for t2 data
                case 2
                    
                    cnt_OFL=0; cnt_MAR=0;  cnt_SYN=0;                   %just counters
                    OverflowCorrection = 0;
                    T2WRAPAROUND=33554432;                              % = 2^25  IMPORTANT! THIS IS NEW IN FORMAT V2.0
                    
                    true_photons=zeros(2,length(buffer_size));          %initialize an array to store true_photon records.
                    
                    fid=fopen(obj.file_path);%opdn the binary .ht2 file
                    
                    fseek(fid,header_size,'bof');%skip over the header to the photon data
                    while 1 %while true
                        
                        batch=fread(fid,buffer_size,'ubit32');%reading in a multiple of 32 bit registers
                        lbatch=length(batch);
                        
                        
                        k=0;%true photon counting variable
                        for i=1:lbatch;%looping over all records in batch
                            
                            %read and decode the 32 bit register of the ith record
                            dtime = bitand(batch(i),33554431);   % the last 25 bits:
                            channel = bitand(bitshift(batch(i),-25),63);   % the next 6 bits:
                            special = bitand(bitshift(batch(i),-31),1);   % the last bit:
                            truetime = OverflowCorrection + dtime;
                            
                            if special == 0   % this means a true 'photon' arrival event.
                                k=k+1;%counting the real photons that we see.
                                true_photons(:,k)=[channel;truetime]; %writing the true photon to a binary array.
                                
                            else    % this means we have a special record; the 'record' is not a 'photon'
                                
                                if channel == 63  % overflow of dtime occured
                                    if(dtime==0) % if dtime is zero it is an old style single oferflow
                                        OverflowCorrection = OverflowCorrection + T2WRAPAROUND;
                                        cnt_OFL=cnt_OFL+1;
                                    else         % otherwise dtime indicates the number of overflows - THIS IS NEW IN FORMAT V2.0
                                        OverflowCorrection = OverflowCorrection + T2WRAPAROUND*dtime;
                                        cnt_OFL=cnt_OFL+dtime;
                                    end;
                                end;
                            end;
                        end
                        fwrite(fpout,true_photons(:,1:k),'uint64');% writing the true photons to the output file in binary.
                        
                        %break the while loop when we have reached the end of the
                        %.ht2 file.
                        if lbatch <buffer_size;
                            break
                        end
                        
                    end
                    fclose(fid);
                    fclose(fpout);
            end
            
            disp('----------------')
            disp(sprintf('%s %s%s', 'Photon records written to',obj.fname,'.photons'))
            disp('----------------')
            
        end
        
        %% photon stream data manipulation/sorting/cleaning.
        function obj=photons_2_channels(obj,file_in_key,file_out_key,nchannels)
            %% creates four .photons output files containing the photon arrival data of each channel.
            %file_in_key is the filename of the .photons file without
            %extension. file_out_key is the filename (-_ch0 ....3) of the new . photons
            %files.n_channels is the number of channels used (2 for PCFS).
            
            buffer_size=obj.buffer_size;
            fid=fopen(strcat(obj.pathstr,file_in_key,'.photons'),'r');
            mode=obj.header.MeasurementMode;
            
            switch mode
                case 3
                    
            if nargin<4
                nchannels=4;
            end
            
            
            %create output files for each channel
            for i=1:nchannels;
                fid_out(i)=fopen(strcat(obj.pathstr,file_out_key,'_ch_',num2str(i),'.photons'),'w');
            end
            
            while 1 %~while true
                batch=fread(fid,[3,buffer_size],'uint64')';%;'%reading in an array of the binary to process.
                lbatch=length(batch);
                
                for k=1:lbatch;
                    fwrite(fid_out(uint64(batch(k,1))+1),batch(k,:)','ubit64');
                end
                if lbatch <buffer_size; %stop reading in batches when we have reached end of .ht3 file.
                    break;
                end
                
            end
            fclose(fid);
            for i=1:nchannels;
                fclose(fid_out(i));
            end
                case 2
               
            if nargin<4
                nchannels=4;
            end
            
            
            %create output files for each channel
            for i=1:nchannels;
                fid_out(i)=fopen(strcat(obj.pathstr,file_out_key,'_ch_',num2str(i),'.photons'),'w');
            end
            
            while 1 %~while true
                batch=fread(fid,[2,buffer_size],'uint64')';%;'%reading in an array of the binary to process.
                lbatch=length(batch);
                
                for k=1:lbatch;
                    fwrite(fid_out(uint64(batch(k,1))+1),batch(k,:)','ubit64');
                end
                
                if lbatch <buffer_size; %stop reading in batches when we have reached end of .ht3 file.
                    break;
                end
                
            end
            
            fclose(fid);
            for i=1:nchannels;
                fclose(fid_out(i));
            end
            
            
            
            end
        end
        function obj=bin_2_decimal(obj,file_in_key,file_out_key)
            %converts any picoquant_bin output file into a decimal (human readable), comma separated
            %file. This file could be parsed to Thomas Bischof's
            %photon-stream analysis package.
            %(https://github.com/tsbischof/photon_correlation)
            % the human-readable photon-arrival time data takes up a lot of
            % disk-space. Storing photon-arrival time data in this
            % way is discouraged. 
            
            %file_in_key: filename of the .photons file without extension.
            %file_out_key: filename of the created .photons_dec
            
            buffer_size = obj.buffer_size;
            
            fid = fopen(strcat(obj.pathstr,file_in_key,'.photons'),'r');
            fid_out = fopen(strcat(obj.pathstr,file_out_key,'.photons_dec'),'w');
            
            switch obj.header.MeasurementMode
                
                case 3
                    while 1 %~while true
                        batch = fread(fid,[3,buffer_size],'uint64')';%reading in an array of the binary to process.
                        fprintf(fid_out,'%2u,%9u,%6u\n',batch');
                        lbatch = length(batch);
                        
                        if lbatch < buffer_size; %stop reading in batches when we have reached end of file.
                            break
                        end
                    end
                    fclose(fid)
                    fclose(fid_out)
                    
                    
                case 2
                    
                    while 1 %~while true
                        batch = fread(fid,[2,buffer_size],'uint64')';%reading in an array of the binary to process.
                        fprintf(fid_out,'%2u,%9u\n',batch');
                        lbatch = length(batch);
                        
                        if lbatch < obj.buffer_size; %stop reading in batches when we have reached end of file.
                            break
                        end
                    end
                    fclose(fid)
                    fclose(fid_out)
            end
        end
        function obj=write_photons_to_one_channel(obj,file_in_key,file_out_key)
            %writing the photons detected on different channels to one
            %channel. Useful to get the autocorrelation of the sum
            %signal for PCFS analysis.
            %open the binary photon file and change the channel to a
            %new number (0) if the photon was detected at one of the
            %channels given in channels_to_combine.
            
            buffer_size=obj.buffer_size;
            fid=fopen(strcat(obj.pathstr,file_in_key,'.photons'),'r');
            fid_out=fopen(strcat(obj.pathstr,file_out_key,'.photons'),'w');
            
            
            switch obj.header.MeasurementMode
                
                case 3
                    while 1
                        batch=fread(fid,[3,buffer_size],'uint64')'; %reading in an array of the binary to process.
                        lbatch=length(batch);
                        batch(:,1)=zeros(lbatch,1);
                        fwrite(fid_out,batch','ubit64');
                        if lbatch <buffer_size; %stop reading in batches when we have reached end of .ht3 file.
                            break;
                        end
                    end
                    
                case 2
                    
                    while 1
                        batch=fread(fid,[2,buffer_size],'uint64')';%reading in an array of the binary to process.
                        lbatch=length(batch);
                        batch(:,1)=zeros(lbatch,1);
                        fwrite(fid_out,batch','ubit64');
                        if lbatch <buffer_size; %stop reading in batches when we have reached end of file.
                            break;
                        end
                    end
                    
            end
            
            fclose(fid);
            fclose(fid_out);
            
        end
        function obj=arrival_time_sorting(obj,file_in_key,file_key_out,tau_window)
            % sort photon data according to photon arrival time.
            % for t2 data the time in ps is the absolute arrival time of
            % the photons
            % for t3 data the time is relative to the sync pulse.
            % A new .photons file is written containing
            % the photons detected withing tau_window
            
            %tau_window: [lower_tau, upper_tau] in ps.
            
            switch obj.header.MeasurementMode
                case 2
                    
                    %counters
                    bin_count=0;
                    run_count=[0 0 0];
                    buffer_size=obj.buffer_size;
                    
                    fid=fopen(strcat(obj.pathstr,file_in_key,'.photons'),'r');
                    fid_out=fopen(strcat(obj.pathstr,file_key_out,'.photons'),'w');
                    
                    while 1
                        batch=fread(fid,[2,buffer_size],'uint64')';%reading in an array of the binary to process.
                        lbatch=length(batch);
                        for i=1:lbatch; %looping over all photons in batch.
                            
                            if batch(i,2)<= tau_window(2) & batch(i,2) > tau_window(1);
                                fwrite(fid_out,batch(i,:),'ubit64');
                            end
                            
                            
                        end
                        if lbatch <buffer_size; %stop reading in batches when we have reached end of .ht2 file.
                            break
                        end
                    end
                    
                case 3
                    
                    %counters
                    bin_count=0;
                    run_count=[0 0 0 0];
                    buffer_size=obj.buffer_size;
                    
                    fid=fopen(strcat(obj.pathstr,file_in_key,'.photons'),'r');
                    fid_out=fopen(strcat(obj.pathstr,file_key_out,'.photons'),'w');
                    
                    while 1
                        batch=fread(fid,[3,buffer_size],'uint64')';%reading in an array of the binary to process.
                        lbatch=length(batch);
                        for i=1:lbatch; %looping over all photons in batch.
                            
                            if batch(i,3)<= tau_window(2) & batch(i,3) > tau_window(1);
                                fwrite(fid_out,batch(i,:),'ubit64');
                            end
                            
                            
                        end
                        if lbatch <buffer_size; %stop reading in batches when we have reached end of .ht3 file.
                            break
                        end
                    end
                    fclose(fid);
                    fclose(fid_out);
            end
        end
        function obj=prep_change_point(obj,file_in_key,file_out_key)
            % Need two column vectors to input into the change point code provided
            % by the Schlau-Cohen group. We need what they call micro and macro
            % time. Microtime is the lifetime of the particle (for each photon).
            % Macrotime is the time since the experiment started when the photon
            % arrives. For t3 data, we need to convert pulse number and rep rate
            % and time after pulse to a macro time. This method preps the matlab
            % object that is needed to run the change point program.

            % Also, it is important to note that you cannot run this on the whole
            % photon stream if it is long. 


            if obj.mode==3;
                fid=fopen(strcat(obj.pathstr,file_in_key,'.photons'),'r');
                buffer_size = obj.buffer_size;
                tt1 = [];
                kin1 = [];

               while 1 %~while true
                        batch=fread(fid,[3,buffer_size],'uint64')';%reading in an array of the binary photons to process.
                        lbatch=length(batch);

                        microtime_temp = batch(:,3);
                        macrotime_temp = batch(:,2)/(obj.sync_rate)*10^15+batch(:,3);

                        tt1 = [tt1; macrotime_temp/(10^12)];
                        kin1 = [kin1; microtime_temp/(10^3)];

                        if lbatch <buffer_size; %stop reading in batches when we have reached end of .ht3 file.
                            break
                        end
               end

               save(file_out_key,'tt1','kin1')
               fclose(fid);
            end
        end

        %% photon stream analysis (intensity, lifetimes, g(n),....)
        function obj=get_intensity(obj,file_in_key,file_out_key,bin_width,lower_limit,upper_limit)
            %compiles and stores the intensity trace as a property of the photons class.
            %file_in_key: filename of the .photons file without ending.
            %file_out_key: name of the .photons property containing the intenisty trace.
            %bin_width: width of the bin for intensity compilation - ps for t2 data; number of pulses for t3 data.
            
            %If upper, lower_limit and file_out_key are given, the photons
            %emitted in bins  with intensities in these bounds are written
            %to a file_out_key.photons file. Important, the upper and lower
            %limit refer to the total counts per bin (all detectors).
            disp('GET_INTENSITY CALLED')
            if nargin<5
                photon_sorting=false;
            else
                photon_sorting=true;
            end
            
            buffer_size=obj.buffer_size;
            
            switch obj.header.MeasurementMode
                case 3
                    
                    
                    if photon_sorting==true
                        
                        %counters.
                        bin_count=0;
                        run_count=[0 0 0 0];%counting variable for all four channels.
                        counts=[];
                        
                        fid=fopen(strcat(obj.pathstr,file_in_key,'.photons'),'r');
                        fid_out=fopen(strcat(obj.pathstr,file_out_key,'.photons'),'w');
                        
                        while 1 %~while true
                            batch=fread(fid,[3,buffer_size],'uint64')';     %reading in an array of the binary photons to process.
                            lbatch=length(batch);
                            indices=[];                                     %vector storing the indices of photons in a particular bin.
                            for i=1:lbatch;
                                while 1
                                    bin=[bin_count,bin_count+1]*bin_width;  %define a new bin in units of pulses
                                    if batch(i,2)<= bin(2) & batch(i,2) > bin(1);
                                        run_count(batch(i,1)+1)=run_count(batch(i,1)+1)+1;% add the photon
                                        indices=[indices;i];%store photons indeces that are in the bin.
                                        break
                                    else %if the photon is not in the bin in units of pulses
                                        bin_count=bin_count+1; %move to the next bin
                                        counts(bin_count,:)=run_count; %add the counts of the previous bin to the counts variable
                                        if sum(run_count)<=upper_limit & sum(run_count)>=lower_limit & ~isempty(indices)
                                            buffer=batch(indices(1):indices(end),:);%buffer the photons from the last bin
                                            fwrite(fid_out,buffer','ubit64'); %buffer needs to be transposed to match dim of records variable in picoquant_bin.
                                        end
                                        run_count=[0 0 0 0];%reset the counters for the next bin.
                                        indices=[];
                                    end
                                    
                                    
                                end
                            end
                            if lbatch <buffer_size; %stop reading in batches when we have reached end of .ht3 file.
                                break
                            end
                        end
                        fclose(fid)
                        fclose(fid_out)
                        
                        %store the counts as property of the photons object
                        if ~isprop(obj,file_out_key)
                            dummy=addprop(obj,file_out_key);
                        end
                        time=linspace(0,length(counts),length(counts))*bin_width;
                        obj.(file_out_key)=struct('trace',counts,'time',time,'bin_width',bin_width,'upper_limit','n.a','lower_limit','n.a');
                    else
                        fid=fopen(strcat(obj.pathstr,file_in_key,'.photons'),'r');
                        buffer_size=obj.buffer_size
                        
                        bin_count=0
                        run_count=[0 0 0 0]
                        counts=[]
                        
                        while 1
                            batch=fread(fid,[3,buffer_size],'uint64')';%reading in an array of the binary to process.
                            lbatch=length(batch);
                            for i=1:lbatch;
                                while 1 % ~ while true
                                    window=[bin_count,bin_count+1]*bin_width; %define a new window in units of pulses
                                    if batch(i,2)<= window(2) & batch(i,2) > window(1);
                                        run_count(batch(i,1)+1)=run_count(batch(i,1)+1)+1;% add the photon
                                        break
                                    else %if the photon is not in the window in units of pulses
                                        bin_count=bin_count+1; %move to the next bin
                                        counts(bin_count,:)=run_count; %add the counts of the previous bin to the counts variable
                                        run_count=[0 0 0 0];%reset the running count to zero; reiterate until the ith photon is fit into a window
                                    end
                                end
                            end
                            if lbatch <buffer_size; %stop reading in batches when we have reached end of .ht3 file.
                                break
                            end
                        end
                        fclose(fid)
                        if ~isprop(obj,file_out_key)
                            dummy=addprop(obj,file_out_key);
                        end
                        time=linspace(0,length(counts),length(counts))*bin_width;
                        obj.(file_out_key)=struct('trace',counts,'time',time,'bin_width',bin_width,'upper_limit','n.a','lower_limit','n.a');
                    end
                    
                    
                    
                    %%dealing with t2 data
                case 2
                    
                    %counters.
                    bin_count=0;
                    run_count=[0 0 0 0];%counting variable for all four channels.
                    counts=[];
                    
                    fid=fopen(strcat(obj.pathstr,file_in_key,'.photons'),'r');
                    fid_out=fopen(strcat(obj.pathstr,file_out_key,'.photons'),'w');
                    if photon_sorting==true
                        
                        while 1 %~while true
                            batch=fread(fid,[2,buffer_size],'uint64')';%reading in an array of the binary photons to process.
                            lbatch=length(batch);
                            indices=[]; %vector storing the indices of photons in a particular bin.
                            for i=1:lbatch;
                                while 1 % ~ while true
                                    bin=[bin_count,bin_count+1]*bin_width; %define a new bin in units of ps
                                    if batch(i,2)<= bin(2) & batch(i,2) > bin(1);
                                        run_count(batch(i,1)+1)=run_count(batch(i,1)+1)+1;% add the photon
                                        indices=[indices;i];%store photons indeces that are in the bin.
                                        break
                                    else %if the photon is not in the bin in units of ps
                                        bin_count=bin_count+1; %move to the next bin
                                        counts(bin_count,:)=run_count; %add the counts of the previous bin to the counts variable
                                        if sum(run_count)<=upper_limit & sum(run_count)>=lower_limit & ~isempty(indices)
                                            buffer=batch(indices(1):indices(end),:);%buffer the photons from the last bin
                                            fwrite(fid_out,buffer','ubit64'); %buffer needs to be transposed to match dim of records variable in picoquant_bin.
                                        end
                                        run_count=[0 0 0 0];%reset the counters for the next bin.
                                        indices=[];
                                    end
                                    
                                    
                                end
                            end
                            if lbatch <buffer_size; %stop reading in batches when we have reached end of .ht2 file.
                                break
                            end
                        end
                        fclose(fid)
                        fclose(fid_out)
                        
                        %store the counts as property of the photons object
                        if ~isprop(obj,file_out_key)
                            dummy=addprop(obj,file_out_key);
                        end
                        time=linspace(0,length(counts),length(counts))*bin_width;
                        obj.(file_out_key)=struct('trace',counts,'time',time,'bin_width',bin_width,'upper_limit','n.a','lower_limit','n.a');
                        
                    else
                        fid=fopen(strcat(obj.pathstr,file_in_key,'.photons'),'r');
                        buffer_size=obj.buffer_size;
                        
                        bin_count=0;
                        run_count=[0 0 0 0];
                        counts=[];
                        
                        while 1 %~while true
                            batch=fread(fid,[2,buffer_size],'uint64')';%reading in an array of the binary to process.
                            lbatch=length(batch);
                            for i=1:lbatch;
                                while 1 % ~ while true
                                    window=[bin_count,bin_count+1]*bin_width; %define a new window in units of ps
                                    if batch(i,2)<= window(2) & batch(i,2) > window(1);
                                        run_count(batch(i,1)+1)=run_count(batch(i,1)+1)+1;% add the photon
                                        break
                                    else %if the photon is not in the window in units of ps
                                        bin_count=bin_count+1; %move to the next bin
                                        counts(bin_count,:)=run_count; %add the counts of the previous bin to the counts variable
                                        run_count=[0 0 0 0];%reset the running count to zero; reiterate until the ith photon is fit into a window
                                    end
                                end
                            end
                            if lbatch <buffer_size; %stop reading in batches when we have reached end of .ht3 file.
                                break
                            end
                        end
                        fclose(fid)
                        if ~isprop(obj,file_out_key)
                            dummy=addprop(obj,file_out_key);
                        end
                        time=linspace(0,length(counts),length(counts))*bin_width;
                        obj.(file_out_key)=struct('trace',counts,'time',time,'bin_width',bin_width,'upper_limit','n.a','lower_limit','n.a');
                    end
            end
        end
        function obj=get_intensity_close_gap(obj,file_in_key,file_out_key,bin_width)
            %takes a photon stream where all photons have been written to
            %the same channel. Function returns a new photons file that has
            %the photon records without all the intensity bins that are
            %zero. This allows for facile correlation-analysis of the
            %blinking dynamics as performed by Rabouw et al.
            
            
            %counters.
            bin_count=0;
            run_count=[0 0 0 0];%counting variable for all four channels.
            counts=[];
            buffer_size=1E7;
            fid=fopen(strcat(obj.pathstr,file_in_key,'.photons'),'r');
            fid_out=fopen(strcat(obj.pathstr,file_out_key,'.photons'),'w');
            N_empty=0;
            
            while 1 %~while true
                batch=fread(fid,[3,buffer_size],'uint64')';     %reading in an array of the binary photons to process.
                lbatch=length(batch);
                indices=[];                                     %vector storing the indices of photons in a particular bin.
                time_difference=[0];
                for i=1:lbatch;
                    while 1
                        bin=[bin_count,bin_count+1]*bin_width;  %define a new bin in units of pulses
                        if batch(i,2)<= bin(2) & batch(i,2) > bin(1);
                            run_count(batch(i,1)+1)=run_count(batch(i,1)+1)+1;% add the photon
                            indices=[indices;i];%store photons indeces that are in the bin.
                            break
                            
                        else %if the photon is not in the bin in units of pulses
                            bin_count=bin_count+1; %move to the next bin
                            counts(bin_count,:)=run_count; %add the counts of the previous bin to the counts variable
                            if sum(run_count)==0
                                N_empty=N_empty+1; % count the number of bins that are empty
                            end
                            if length(indices)>1
                                buffer=batch(indices(1):indices(end),:);%buffer the photons from the last bin
                                buffer(:,2)=buffer(:,2)-N_empty*bin_width; % reduce the time of the last bin's photons by the time occupied by the empty bins.
                                fwrite(fid_out,buffer','ubit64'); %buffer needs to be transposed to match dim of records variable in picoquant_bin.
                            end
                            run_count=[0 0 0 0];%reset the counters for the next bin.
                            indices=[];
                        end
                    end
                end
                if lbatch <buffer_size; %stop reading in batches when we have reached end of .ht3 file.
                    break
                end
            end
            fclose(fid)
            fclose(fid_out)
            
            %store the counts as property of the photons object
            if ~isprop(obj,file_out_key)
                dummy=addprop(obj,file_out_key);
            end
            time=linspace(0,length(counts),length(counts))*bin_width;
            obj.(file_out_key)=struct('trace',counts,'time',time,'bin_width',bin_width,'upper_limit','n.a','lower_limit','n.a');
        end

        function rmv_AP(obj,file_in_key,file_out_key)
            %removes afterpulsing from t3 data. This is accomplished by
            %removing any second photons detected after the same sync pulse
            %for each detection channel.
            
            buffer_size=obj.buffer_size; %number of photons processed per batch
            
            switch obj.header.MeasurementMode 
                case 2
                    
                    disp('--------------')
                    disp('Only for t3 data')
                    disp('--------------')
                    
                case 3
                    
                    fid=fopen(strcat(obj.pathstr,file_in_key,'.photons'),'r'); % open the photon stream. this should not be sorted by channel
                    fid_out=fopen(strcat(obj.pathstr,file_out_key,'.photons'),'w'); %open the file which we write the new photons to
                    
                    while 1
                        batch=fread(fid,[3,buffer_size],'uint64')';%reading in an array of the binary to process.
                        lbatch=length(batch);
                        
                        i=1;
                        while i<=lbatch-3; %looping over the photons.
                            clear channels
                            channels=batch(i,1);
                            
                            if batch(i,2)==batch(i+1,2) %% is the next photon from the same excitation pulse?
                                channels=[channels,batch(i+1)];
                                photons_after_same_pulse=2;
                                if batch(i,2)==batch(i+2,2) % and the next one?
                                    channels=[channels,batch(i+2)];
                                    photons_after_same_pulse=photons_after_same_pulse+1;
                                    if batch(i,2)==batch(i+3,2) % and the one after that?
                                        channels=[channels,batch(i+3)];
                                        photons_after_same_pulse=photons_after_same_pulse+1;
                                    end
                                end
                            else %photon is detected after the next excitation pulse - move on.
                                photons_after_same_pulse=1;
                            end
                            
                            
                            %% write photons detected only on different channels, but after the same excitation pulse to file.
                            used_channels=[];
                            for k=1:photons_after_same_pulse %looping over the number of photons after the same pulse
                                if ismember(batch(i+k-1,1),used_channels)==0
                                    fwrite(fid_out,batch(i+k-1,:),'ubit64');
                                    used_channels=[used_channels,batch(i+k-1,1)];
                                end
                            end
                            
                            i=i+photons_after_same_pulse; %jump over the photons that were recorded after the same excitation pulse.
                            
                        end
                        if lbatch <buffer_size; %stop reading in batches when we have reached end of .ht3 file.
                            break
                        end
                    end
                    fclose(fid)
                    fclose(fid_out)
            end
        end
        function obj=lifetime_histo(obj,file_in_key,file_out_key,resolution)
            
                % histograms the lifetime of a .ht3 file with a given resolution.  
                % the given resolution should be a multiple of the original resolution used
                % for the measurement. For instance, if the measurement
                % rsolution was 64 ps, then the resolution to form the
                % histogram of the photon-arrival records could be 128,
                % 256, 384 or 512 ps ...
                
            switch obj.header.MeasurementMode
                case 2
                    disp('---------------------------------------------')
                    disp('This function is not available for t2 data.')
                    disp('---------------------------------------------')
                    
                case 3
                    
                    %check if given resultion is a multiple of the
                    %original resolution
                    foo=resolution/obj.header.Resolution;
                    
                    if floor(foo)~=foo
                        disp('----------------------')
                        disp('The given resolution must be a multiple of the original resolution')
                        disp('Check original resolution as obj.header.Resolution')
                        disp('----------------------')
                        return 
                    end
                   
                    buffer_size=obj.buffer_size;
                    
                    %counters for the bin and for the four possible
                    %channels.
                    bin_count=0;
                    run_count=[0 0 0 0];
                    
                    %initializations
                    rep_time = 1/obj.header.SyncRate*10^12; % in ps. 
                    n_bins = floor(rep_time/resolution);
                    
                    time = [1/2:n_bins-1/2]*resolution;
                    histo = zeros(size(time));
                                        
                    fid=fopen(strcat(obj.pathstr,file_in_key,'.photons'),'r');
                    
                    %we just compile a histogram of all photons without temporal binning.
                    while 1
                        batch=fread(fid,[3,buffer_size],'uint64')'; %reading in an array of the binary to process.
                        lbatch=length(batch);
                        
                        histo=histo+hist(batch(:,3),time);
                        
                        if lbatch <buffer_size;        % stop reading when reached the end of the file.
                            break
                        end
                    end
                    
                    fclose(fid)
             
                    %store returns as properties of the photons object.
                    if ~isprop(obj,file_out_key)
                        dummy=addprop(obj,file_out_key)
                    end
                    obj.(file_out_key)=struct('time',time,'lifetime',histo);
            end
        end
        function obj=FLID(obj,file_in_key,file_out_key,resolution,pulses_per_bin,Sync_Rate)
            
                % histograms the lifetime of a .ht3 file with a given resolution for 
                % chunks of the photon-stream with a length specified in
                % pulses_per_bin. The histogram of photon-arrival times and
                % the intensity for each chunk of the photon-stream is
                % returned to allow fluorescence lifetime intensity
                % distribution analysis. 

                % the given resolution should be a multiple of the original resolution used
                % for the measurement. For instance, if the measurement
                % rsolution was 64 ps, then the resolution to form the
                % histogram of the photon-arrival records could be 128,
                % 256, 384 or 512 ps ...
            
            switch obj.header.MeasurementMode
                case 2
                    disp('---------------------------------------------')
                    disp('This function is not available for t2 data.')
                    disp('---------------------------------------------')
                    
                case 3
                    
                    %check if given resultion is a multiple of the
                    %original resolution
                    foo=resolution/obj.header.Resolution;
                    
                    if floor(foo)~=foo
                        disp('----------------------')
                        disp('The given resolution must be a multiple of the original resolution')
                        disp('Check original resolution as obj.header.Resolution')
                        disp('----------------------')
                        return
                    end
                    
                    buffer_size=obj.buffer_size;
                    
                    %counters for the bin and for the four possible
                    %channels.
                    bin_count=0;
                    run_count=[0 0 0 0];
                    
                    %initializations
                    rep_time = 1/Sync_Rate*1E12; % in ps.
                    n_bins = floor(rep_time/resolution);
                    
                    time = [1/2:n_bins-1/2]*resolution;
                    histo = zeros(length(time),1E4);
                    
                    fid=fopen(strcat(obj.pathstr,file_in_key,'.photons'),'r');
                    
                    
                    while 1
                        batch=fread(fid,[3,buffer_size],'uint64')';%reading in an array of the binary to process.
                        
                        lbatch=length(batch);
                        bin_count=0;
                        photons=[];
                        
                        
                        
                        for i=1:lbatch; %looping over all photons in batch.
                            while 1 % ~ while true
                                
                                bin=[bin_count,bin_count+1]*pulses_per_bin;             %define a new window in units of pulses
                                
                                if batch(i,2)<= bin(2) & batch(i,2) > bin(1);
                                    run_count(batch(i,1)+1)=run_count(batch(i,1)+1)+1;  % add the photon to intensity counter.
                                    photons=[photons;batch(i,:)];
                                    break
                                else                                                    %if the photon is not in the window in units of pulses
                                    bin_count=bin_count+1;                              %move to the next bin
                                    counts(bin_count,:)=run_count;                      %add the counts of the previous bin to the counts variable
                                    histo(:,bin_count)=hist(photons(:,3),time)';
                                end
                                run_count=[0 0 0 0];%reset counter.
                                photons=[]; %reset photons buffer.
                            end
                        end
                        
                        
                        if lbatch <buffer_size; %stop reading in batches when we have reached end of .ht3 file.
                            break
                        end
                    end
                    
                    histo=histo(:,1:bin_count); % only return non-zero elements of the histogram matrix. 
                    
                    %store returns as properties of the photons object.
                    if ~isprop(obj,file_out_key)
                        dummy=addprop(obj,file_out_key)
                    end
                    obj.(file_out_key)=struct('time',time,'lifetime',histo,'intensity',counts);
            end
        end
        function obj=get_g2(obj,file_in_key,file_out_key,n_pulses,time_bins,rrate)
            
            switch obj.header.MeasurementMode
                case 2
                    disp('---------------------------')
                    disp('Only available for T3 data')
                    disp('---------------------------')
                    
                case 3
            %rrate = obj.sync_rate;
            %time_bins = 1/obj.sync_rate;
            time_bins = time_bins;
            
            buffer_size=obj.buffer_size; %number of photons processed per batch
            fid=fopen(strcat(obj.pathstr,file_in_key,'.photons'),'r'); % open the phooton stream. this should not be sorted by channel
           
            u = 0;
            %g2 = zeros(2*n_pulses+1,channels,channels,time_bins*2+1);
            
            while 1 %~while true
                batch=fread(fid,[3,buffer_size],'uint64')';%reading in an array of the binary to process.
                u = u+1; %counter for batches
                %batch_overflow(4:6,:) = batch(1:3,:);
                
                lbatch=length(batch);
                clear g2_temp
                
                ii = 0;
                
                while 1
                    ii = ii +1;
                    batch_shift_down = batch;
                    batch_shift_down(1:ii,:) = [];
                    %size(batch_shift_down)
                    batch_shift_up = batch;
                    batch_shift_up(end-ii+1:end,:) = [];
                    %size(batch_shift_up)
                    batch_diff = batch_shift_down-batch_shift_up;
                    pulse_sorted = batch_diff;
                    batch_diff_extended = [batch_diff,batch_shift_down(:,3),batch_shift_up(:,3)];
                    Logic1 = batch_diff(:,2) > n_pulses; %If the pulse separation is greater than the max
                    Logic2 = batch_diff(:,2) < -n_pulses; %if the pulse separation is less than the neg of the max
                    
                    Logic = Logic1 | Logic2; %Add up all the conditions (could easily add more). This just cuts off pulse separation
                    pulse_sorted(Logic,:) = []; %Remove those rows
                    pulse_sorted_AP = pulse_sorted;
                    Logic3 = pulse_sorted(:,1) == 0; %If they arrive at the same detector
                    %                     Logic4 = pulse_sorted(:,2) == 0; %Same pulse
                    Logic_AP = Logic3; %& Logic4;
                    pulse_sorted_AP(Logic_AP,:) = []; %Removes After Pulsing
                    
                    %pulse_sorted is now a matrix containing difference in
                    %detector number, difference in time, and pulse
                    %difference. Can now turn pulse difference and time
                    %difference into total time separation and histogram
                    %and plot.
                    pulse_sorted_time = pulse_sorted_AP(:,3)+1/rrate*10^12*(pulse_sorted_AP(:,2));
                    g2_temp(:,ii) = hist(pulse_sorted_time,time_bins*(2*n_pulses));
                    if isempty(pulse_sorted_AP)
                        break;
                    end
                    
                    
                    
                end
                %size(g2_temp)
                g2_batch(u,:) = sum(g2_temp,2);
                
                
                
                
                if lbatch <buffer_size; %stop reading in batches when we have reached end of .ht3 file.
                    break;
                end
            end
            g2 = sum(g2_batch,1);
            size(g2)
            trep = 1/rrate;
            time = linspace(0,2*trep,time_bins*2);
            size(time)
            figure;
            plot(time(2:end),g2(2:end),'b')
            hold on
            plot(-time,fliplr(g2'),'b')
            center = 2*trapz(g2(1:end/4));
            side = trapz(g2(end/4:end));
            ratio = center/side
            noise = g2(7*end/8:8*end/8);
            average_noise = mean(noise);
            g2_noise_sub = g2-average_noise;
            figure;
            plot(time,g2_noise_sub)
            center_noise = 2*trapz(g2_noise_sub(1:end/4));
            side_noise = trapz(g2_noise_sub(end/4:end));
            ratio_noise = center_noise/side_noise
            if ~isprop(obj,file_out_key)
                dummy=addprop(obj,file_out_key)
            end
            obj.(file_out_key)=struct('time',time,'g2',g2)
            end
            
        end
        function obj=get_g2_long_time(obj,file_in_key,file_out_key,n_pulses,time_bins,rrate)
            
             switch obj.header.MeasurementMode
                case 2
                    disp('---------------------------')
                    disp('Only available for T3 data')
                    disp('---------------------------')
                    
                case 3
            
            
            
            %rrate = obj.sync_rate;
            %time_bins = 1/obj.sync_rate;
            time_bins = time_bins;
            
            buffer_size=obj.buffer_size; %number of photons processed per batch
            fid=fopen(strcat(obj.pathstr,file_in_key,'.photons'),'r'); % open the phooton stream. this should not be sorted by channel
            u = 0;
            %g2 = zeros(2*n_pulses+1,channels,channels,time_bins*2+1);
            
            while 1 %~while true
                batch=fread(fid,[3,buffer_size],'uint64')';%reading in an array of the binary to process.
                u = u+1; %counter for batches
                %batch_overflow(4:6,:) = batch(1:3,:);
                
                lbatch=length(batch);
                clear g2_temp
                
                ii = 30000;
                
                while 1
                    clear ind
                    %                     disp('Start Loop')
                    ii = ii +1;
                    batch_shift_down = batch;
                    batch_shift_down(1:ii,:) = [];
                    %size(batch_shift_down)
                    batch_shift_up = batch;
                    batch_shift_up(end-ii+1:end,:) = [];
                    %size(batch_shift_up)
                    batch_diff = batch_shift_down-batch_shift_up;
                    pulse_sorted = batch_diff;
                    batch_diff_extended = [batch_diff,batch_shift_down(:,3),batch_shift_up(:,3)];
                    Logic1 = batch_diff(:,2) ~= n_pulses; %If the pulse separation is greater than the max
                    
                    
                    Logic = Logic1; %Add up all the conditions (could easily add more). This just cuts off pulse separation
                    pulse_sorted(Logic,:) = []; %Remove those rows
                    pulse_sorted_AP = pulse_sorted;
                    Logic3 = pulse_sorted(:,1) == 0; %If they arrive at the same detector
                    %                     Logic4 = pulse_sorted(:,2) == 0; %Same pulse
                    Logic_AP = Logic3; %& Logic4;
                    pulse_sorted_AP(Logic_AP,:) = []; %Removes After Pulsing
                    
                    %pulse_sorted is now a matrix containing difference in
                    %detector number, difference in time, and pulse
                    %difference. Can now turn pulse difference and time
                    %difference into total time separation and histogram
                    %and plot.
                    pulse_sorted_time = pulse_sorted_AP(:,3)+1/rrate*10^12*(pulse_sorted_AP(:,2));
                    g2_temp(:,ii) = hist(pulse_sorted_time,time_bins*2);
                    ind = find(batch_diff(:,2)<=n_pulses,1);
                    disp(ind)
                    disp(length(batch_diff))
                    %                     pause;
                    if isempty(ind) || length(batch_diff) == 1
                        disp('loop')
                        break;
                    end
                    %                     disp('Pass If')
                    
                    
                end
                %size(g2_temp)
                g2_batch(u,:) = sum(g2_temp,2);
                
                
                
                
                if lbatch <buffer_size; %stop reading in batches when we have reached end of .ht3 file.
                    break;
                end
            end
            
            g2 = sum(g2_batch,1);
            size(g2)
            trep = 1/rrate;
            time = linspace(0,trep*2,time_bins*2);
            size(time)
            figure;
            plot(time(2:end),g2(2:end),'b')
            hold on
            plot(-time,fliplr(g2'),'b')
            noise = g2(7*end/8:8*end/8);
            average_noise = mean(noise);
            g2_noise_sub = g2-average_noise;
            figure;
            plot(time,g2_noise_sub)
            
            if ~isprop(obj,file_out_key)
                dummy=addprop(obj,file_out_key)
            end
            obj.(file_out_key)=struct('time',time,'g2',g2)
            
             end
        end
        function obj=get_g3(obj,file_in_key,file_out_key)
            
             switch obj.header.MeasurementMode
                case 2
                    disp('---------------------------')
                    disp('Only available for T3 data')
                    disp('---------------------------')
                    
                case 3
            
            
            %Run remove After Pulsing first
            %This code outputs the pulse integrated information needed to
            %calculate a triexciton quantum yield. It could be made to
            %create an entire 2D map, but we chose that this was simpler
            %and clearer.
            
            buffer_size=obj.buffer_size; %number of photons processed per batch
            fid=fopen(strcat(obj.pathstr,file_in_key,'.photons'),'r'); % open the photon stream. this should not be sorted by channel
            %originate counting variables
            u = 0;
            count = 0;
            g3c = 0;
            g3s = 0;
            
            while 1 %~while true
                batch=fread(fid,[3,buffer_size],'uint64')';%reading in an array of the binary to process.
                %                 size(batch)
                lbatch = length(batch);
                u = u+1; %counter for batches
                %                 nbins = batch(end,2)-batch(1,2)+1
                %                 hist_pulse = histogram(batch(:,2),nbins); %histogram the number of photons per pulse
                %                 mult_pulses = hist_pulse.Values; %extract the values from the histogram
                %                 index_mult = find(mult_pulses==3); %find which have more than 1 photon
                %Create vectors to manipulate
                batch_shift_1 = batch;
                batch_shift_2 = batch;
                %Shift batches so as to compare side by side photons across
                %the whole batch simultaneously
                batch_shift_1(1,:) = [];
                batch_shift_2(1:2,:) = [];
                
                %Logic statements are used to sort for conditions we want
                Logic1 = batch(1:end-2,2) ~= batch_shift_1(1:end-1,2); %same pulse photon 1 and 2
                Logic2 = batch(1:end-2,2) ~= batch_shift_2(:,2); %same pulse photon 2 and 3
                
                Logic1_1p = batch(1:end-2,2) ~= batch_shift_1(1:end-1,2);%same pulse photon 1 and 2
                Logic2_1p = batch(1:end-2,2) ~= batch_shift_2(:,2)-1;%subsequent pulse photon 2 and 3
                
                Logic_det1 = batch(1:end-2,1) ~= batch_shift_1(1:end-1,1);%different detectors
                Logic_det2 = batch(1:end-2,1) ~= batch_shift_2(1:end,1);%different detectors
                
                Logic_det = Logic_det1 | Logic_det2; %combine conditions
                
                Logic = Logic1 | Logic2;
                Logic_1p = Logic1_1p | Logic2_1p | Logic_det;
                
                sort_g3c = batch(1:end-2,:);
                sort_g3c(Logic,:) = []; %remove photons that do not meet the conditions
                sort_g3s = batch(1:end-2,:);
                sort_g3s(Logic_1p,:) = [];
                
                %Add up counts
                count_center = length(sort_g3c);
                count_side = length(sort_g3s);
                
                %Add temporary count to total count
                g3c = g3c+count_center
                g3s = g3s+count_side
                
                
                if lbatch <buffer_size; %stop reading in batches when we have reached end of .ht3 file.
                    break;
                end
            end
            
            if ~isprop(obj,file_out_key)
                dummy=addprop(obj,file_out_key);
            end
            center_counts = g3c;
            side_counts = g3s;
            TXQY = center_counts/side_counts;
            error_side = sqrt(side_counts);
            error_center = sqrt(center_counts);
            error = TXQY*(error_side/side_counts+error_center/center_counts);
            
            
            obj.(file_out_key)=struct('center',center_counts,'side',side_counts,'TX_QY',TXQY,'error',error);
            fclose(fid);
             end
        end
        function obj=get_PNRL2(obj,file_in_key,file_out_key)
            
             switch obj.header.MeasurementMode
                case 2
                    disp('---------------------------')
                    disp('Only available for T3 data')
                    disp('---------------------------')
                    
                case 3
            
            
            test_counter = 0;
            
            buffer_size=obj.buffer_size; %number of photons processed per batch
            fid=fopen(strcat(obj.pathstr,file_in_key,'.photons'),'r'); % open the phooton stream. this should not be sorted by channel
            u = 0;
            %g2 = zeros(2*n_pulses+1,channels,channels,time_bins*2+1);
            
            %             true_nbins=2^15/bin_divisor;%the true number of bins is the systems bin depth divided by an integer defined in bin_divisor to avoid artifacts.
            %             bin_vector=linspace(0.5,true_nbins,(true_nbins-0.5))*obj.resolution*bin_divisor;% center positions for bins.
            %             max_time = 1/obj.sync_rate*(10^15);
            %             time=linspace(0,max_time,2^15);
            rep_rate = 1/obj.sync_rate*10^15;
            bins = rep_rate/(obj.resolution);
            bin_vector = linspace(obj.resolution/2,rep_rate-obj.resolution/2,bins);
            %true_nbins=2^15;%/bin_divisor;%the true number of bins is the systems bin depth divided by an integer defined in bin_divisor to avoid artifacts.
            %bin_vector=linspace(0.5,true_nbins,(true_nbins-0.5))*obj.resolution;% center positions for bins.
            time=bin_vector;%
            size(time)
            figure;
            plot(time)
            PNRL21 = zeros(size(time));
            PNRL22 = zeros(size(time));
            PNWT22 = zeros(size(time));
            
            while 1 %~while true
                batch=fread(fid,[3,buffer_size],'uint64')';%reading in an array of the binary to process.
                u = u+1; %counter for batches
                %batch_overflow(4:6,:) = batch(1:3,:);
                
                lbatch=length(batch);
                clear g2_temp
                
                ii = 0;
                photons_temp = [];
                photons_ext_temp = [];
                
                while 1
                    clear pulse_sorted_AP
                    clear batch_diff
                    ii = ii +1;
                    batch_shift_down = batch;
                    batch_shift_down(1:ii,:) = [];
                    %size(batch_shift_down)
                    batch_shift_up = batch;
                    batch_shift_up(end-ii+1:end,:) = [];
                    %size(batch_shift_up)
                    batch_diff = batch_shift_down-batch_shift_up;
                    pulse_sorted = batch_diff;
                    batch_diff_extended = [batch_diff,batch_shift_down(:,2:3),batch_shift_up(:,3)];
                    batch_diff_extended_t1 = batch_diff_extended;
                    %                     size(batch_diff_extended)
                    Logic1 = batch_diff(:,2) ~= 0; %If the photons don't arrive on the same pulse
                    
                    Logic = Logic1; %Add up all the conditions (could easily add more). This just cuts off pulse separation
                    pulse_sorted(Logic,:) = []; %Remove those rows so we only have multiple photons from the same pulse
                    batch_diff_extended_t1(Logic,:) = [];
                    pulse_sorted_AP = pulse_sorted;
                    Logic3 = pulse_sorted(:,1) == 0; %If they arrive at the same detector
                    %                     Logic4 = pulse_sorted(:,2) == 0; %Same pulse
                    Logic_AP = Logic3; %& Logic4;
                    pulse_sorted_AP(Logic_AP,:) = []; %Removes After Pulsing
                    batch_diff_extended_t1(Logic_AP,:) = [];
                    
                    %pulse_sorted is now a matrix containing difference in
                    %detector number, difference in time, and pulse
                    %difference. Can now turn pulse difference and time
                    %difference into total time separation and histogram
                    %and plot.
                    photons_ext_temp = [photons_ext_temp;batch_diff_extended_t1];
                    size(photons_ext_temp);
                    photons_temp = [photons_temp;pulse_sorted_AP];
                    if isempty(pulse_sorted_AP)
                        break;
                    end
                    
                end
                unique_pulse = unique(photons_temp);%determine which pulses had multiple photons
                unique_pulse_ext = unique(photons_ext_temp(:,4));
                test_counter = test_counter+ size(unique_pulse_ext);
                [count_pulses dummy] = size(unique_pulse_ext); %how many unique pulses
                PNRL21_batch = [];
                PNRL22_batch = [];
                PNWT22_batch = [];
                for jj = 1:count_pulses
                    clear num_ph_temp
                    PNRL21_temp = [];
                    PNRL22_temp = [];
                    PNWT22_temp = [];
                    num_ph_temp = nnz(photons_ext_temp(:,4) == unique_pulse_ext(jj));
                    if num_ph_temp == 1
                        index_temp = find(photons_ext_temp(:,4) == unique_pulse_ext(jj),1);
                        PNRL21_temp = photons_ext_temp(index_temp,6);
                        PNRL22_temp = photons_ext_temp(index_temp,5);
                        PNWT22_temp = photons_ext_temp(index_temp,3);
                    end
                    PNRL21_batch = [PNRL21_batch;PNRL21_temp];
                    PNRL22_batch = [PNRL22_batch;PNRL22_temp];
                    PNWT22_batch = [PNWT22_batch;PNWT22_temp];
                    clear PNRL21_temp
                    clear PNRL22_temp
                    clear PNWT22_temp
                end
                PNRL21_temp_hist = hist(PNRL21_batch,time);
                PNRL22_temp_hist = hist(PNRL22_batch,time);
                PNWT22_temp_hist = hist(PNWT22_batch,time);
                
                
                PNRL21 = PNRL21+PNRL21_temp_hist;
                PNRL22 = PNRL22+PNRL22_temp_hist;
                PNWT22 = PNWT22+PNWT22_temp_hist;
                sum(PNRL21)
                
                
                
                if lbatch <buffer_size; %stop reading in batches when we have reached end of .ht3 file.
                    break;
                end
            end
            
            PNRL21 = PNRL21(1:bins/1000);
            PNRL22 = PNRL22(1:bins/1000);
            PNWT22 = PNWT22(1:bins/1000);
            time = time(1:bins/1000);
            
            if ~isprop(obj,file_out_key)
                dummy=addprop(obj,file_out_key)
            end
            obj.(file_out_key)=struct('time',time,'PNRL21',PNRL21,'PNRL22',PNRL22,'PNWT22',PNWT22)
            figure;
            semilogy(time/1000,PNRL21)
            hold on
            semilogy(time/1000,PNRL22)
            semilogy(time/1000,PNWT22)
            
             end
        end
        function obj=get_PNRL3_test(obj,file_in_key,file_out_key)
            
             switch obj.header.MeasurementMode
                case 2
                    disp('---------------------------')
                    disp('Only available for T3 data')
                    disp('---------------------------')
                    
                case 3
            
            %Run remove After Pulsing first
            
            buffer_size=obj.buffer_size; %number of photons processed per batch
            fid=fopen(strcat(obj.pathstr,file_in_key,'.photons'),'r'); % open the phooton stream. this should not be sorted by channel
            u = 0;
            count = 0;
            PNRL301_string = [];
            PNRL302_string = [];
            PNRL303_string = [];
            
            while 1 %~while true
                batch=fread(fid,[3,buffer_size],'uint64')';%reading in an array of the binary to process.
                %                 size(batch)
                lbatch = length(batch);
                u = u+1; %counter for batches
                %                 nbins = batch(end,2)-batch(1,2)+1
                %                 hist_pulse = histogram(batch(:,2),nbins); %histogram the number of photons per pulse
                %                 mult_pulses = hist_pulse.Values; %extract the values from the histogram
                %                 index_mult = find(mult_pulses==3); %find which have more than 1 photon
                batch_shift_1 = batch;
                batch_shift_2 = batch;
                %                 size(batch_shift_1)
                %                 size(batch_shift_2)
                batch_shift_1(1,:) = [];
                batch_shift_2(1:2,:) = [];
                Logic1 = batch(1:end-2,2) ~= batch_shift_1(1:end-1,2);
                Logic2 = batch(1:end-2,2) ~= batch_shift_2(:,2);
                
                Logic = Logic1 | Logic2;
                %size(Logic)
                %size(batch)
                %size(batch_shift_1)
                %size(batch_shift_2)
                sort_PNRL301 = batch(1:end-2,:);
                sort_PNRL301(Logic,:) = [];
                sort_PNRL302 = batch_shift_1(1:end-1,:);
                sort_PNRL302(Logic,:) = [];
                sort_PNRL303 = batch_shift_2;
                sort_PNRL303(Logic,:) = [];
                
                PNRL301_string = [PNRL301_string; sort_PNRL301];
                PNRL302_string = [PNRL302_string; sort_PNRL302];
                PNRL303_string = [PNRL303_string; sort_PNRL303];
                
                %                 [a,b] = size(index_mult); %determine number of pulses that exist with three photons
                %
                %                 for ii = 1:b
                %                     count = count+1
                %                     pulse_index = find(batch(:,2)==index_mult,1);
                %                     PNRL301_string(count,:) = batch(pulse_index,3);
                %                     PNRL302_string(count,:) = batch(pulse_index+1,3);
                %                     PNRL303_string(count,:) = batch(pulse_index+2,3);
                %                     clear pulse_index
                %                 end
                
                
                if lbatch <buffer_size; %stop reading in batches when we have reached end of .ht3 file.
                    break;
                end
            end
            PNRL312_string = PNRL302_string-PNRL301_string;
            PNRL323_string = PNRL303_string-PNRL302_string;
            disp(size(PNRL301_string))
            disp(size(PNRL302_string))
            disp(size(PNRL303_string))
            disp(size(PNRL312_string))
            disp(size(PNRL323_string))
            rep_rate = 1/(obj.sync_rate+1)*10^12
            bins = rep_rate/(obj.resolution)
            time = linspace(obj.resolution/2,rep_rate-obj.resolution/2,bins);
            PNRL301 = hist(PNRL301_string(:,3),time);
            PNRL302 = hist(PNRL302_string(:,3),time);
            PNRL303 = hist(PNRL303_string(:,3),time);
            PNRL312 = hist(PNRL312_string(:,3),time);
            PNRL323 = hist(PNRL323_string(:,3),time);
            
            if ~isprop(obj,file_out_key)
                dummy=addprop(obj,file_out_key);
            end
            figure;
            semilogy(time/1000,PNRL301)
            hold on
            semilogy(time/1000,PNRL312)
            semilogy(time/1000,PNRL323)
            hold off
            
            obj.(file_out_key)=struct('PNRL301',PNRL301,'PNRL302',PNRL302,'PNRL303',PNRL303,'time',time/1000,'PNRL312',PNRL312,'PNRL323',PNRL323);
            fclose(fid);
             end
        end
        function obj = photon_corr(obj, file_in_key,file_out_key,channels, time_bounds, lag_precision, lag_offset, use_dll,Sync_Rate)
            %% Contributed by Boris Spokoyny. Allow to correlate the photon-stream on a log timescale.
            % file_in key: file ID of the photons-file to correlate
            % file_out_key: ID for the struct array added to the photons
            % object containing the correlation.
            % channels: Hydraharp channels to be correlated. e.g. [0,1] for
            % cross-correlation of channels 0 and 1. 
            % time-bound: upper and lower limit for the correlation. In ps
            % for T2, in pulses for T3. 
            %lag_precision: Number of correlation points between
            %time-spacings of log(2).
            %lag_offset: offset in ps or pulses between the channels. 
            %use_dll for using a C++ code for the correlation algorithm for
            %improved efficiency. 
            
            dot_photons_handle = fopen(strcat(obj.pathstr,file_in_key,'.photons'),'r');
            
            switch obj.header.MeasurementMode
                case 2 %T2 mode
                    data = fread(dot_photons_handle, [2, obj.header.nRecords], 'uint64');
                case 3 %T3 mode
                    data = fread(dot_photons_handle, [3, obj.header.nRecords], 'uint64');
            end
            fclose(dot_photons_handle);
            
            %split inot channels
            ch0 = data(2, data(1, :) == channels(1)); %ch0 syncs
            ch1 = data(2, data(1, :) == channels(2)); % ch1 syncs
            
            start_time = time_bounds(1);
            stop_time = time_bounds(2);
            if use_dll
                [lags, corr_norm] = cross_corr_dll(ch0,ch1, start_time, stop_time,...
                    lag_precision, lag_offset);
                if obj.header.MeasurementMode == 3
                    sync_period = 1E12/Sync_Rate;      % in picoseconds
                    lags = lags.*sync_period;
                    [~, start_bin] = min(abs(lags-start_time.*sync_period));
                    [~, stop_bin] = min(abs(lags-stop_time.*sync_period));
                    lags = lags(start_bin:stop_bin);
                    corr_norm = corr_norm(start_bin:stop_bin);
                end
            else
                [lags, corr_norm] = cross_corr(ch0,ch1, start_time, stop_time, ...
                    lag_precision, lag_offset);
                if obj.header.MeasurementMode == 3
                    sync_period = 1E12/Sync_Rate;      % in picoseconds
                    lags = lags.*sync_period;
                end
            end
            
            if ~isprop(obj,file_out_key)
                dummy=addprop(obj,file_out_key)
            end
            obj.(file_out_key)=struct('lags',lags,'corr_norm',corr_norm);
            
            %helper functions.
            function [lags, corr_norm] = cross_corr(ch1,ch2, start_time, stop_time, coarseness, offset_lag)
                %this algorithm computes the cross-correlation between ch1 and ch2
                %variables.
                %-- For T2 data, ch1 and ch2 correspond to absolute photon arrival
                %   times.
                %-- For T3 data, ch1 and ch2 should correspond to the photon arrival
                %   sync pulse number.
                %-- start_time and stop_time for T2 data should be in time units of the
                %-- Photon arrival times, and for T3 data should be in units of sync pulses.
                %-- The correlation lags are log(2) spaced with coarseness # of lags per
                %   cascade, i.e. if start_time = 1; stop_time = 50; coarseness = 4; the
                %   lag bins will be [ 1,2,3,4 ; 6,8,10,12 ; 16,20,24,28; 36,44 ] if
                %   coarseness = 2, the lag bins become [ 1,2 ; 4,6 ; 10,14 ; 22,30 ; 46 ]
                %-- The algorithm works by taking each photon arrival time and counting the
                %   number of photons that are lag bins away. For example say our lag bin
                %   edges are [1, 2, 4, 6, 10, 14]
                
                %   Time Slot: 1   2   3   4   5   6   7   8   9   10  11  12  13  14  15
                %   Photon?    1   1   0   0   1   1   0   1   1   0   1   0   0   0   1
                %   1st Photon ^
                %   Lag bins   |   |       |       |               |               |
                %   # photons    1     0       2           2                1
                %   2nd Photon     ^
                %   Lag bins       |   |       |       |               |               |
                %   # photons        0     1       1           3               1
                %   3rd Photon                 ^
                %   Lag bins                   |   |       |       |               |               |
                %   # photons                    1     2       1           1               1
                
                %etc..
                % the cross-correlation is the sum of all the photons for each bin, i.e
                % for the three runs above we get [2, 3, 4, 6, 3].
                
                % generates lags bins with log2 spacing
                % Here stop time, refers to a lag time in picoseconds for t2 data or
                % maximum number of pulses for t3 data.
                [lag_bin_edges, lags, division_factor] = generate_log2_lags(stop_time, coarseness);
                
                %find the bin region
                [~, start_bin] = min(abs(lag_bin_edges-start_time));
                [~, stop_bin] = min(abs(lag_bin_edges-stop_time));
                lag_bin_edges = lag_bin_edges(start_bin:stop_bin); %bins
                lags = lags((start_bin+1):stop_bin); % centers of bins
                division_factor = division_factor((start_bin+1):stop_bin); %normalization factor
                
                %%counters etc for normalization
                ch1_min = inf; ch2_min = inf; %minimum time tag
                ch1_count = 0; ch2_count = 0; %photon tally in each channel
                num_ch1 = numel(ch1); num_ch2 = numel(ch2); %number of photons in each channel
                ch1_count = ch1_count+num_ch1; ch2_count = ch2_count+num_ch2;
                ch1_min = min(ch1_min, min(ch1)); ch2_min = min(ch2_min, min(ch2));
                
                %%Correlation
                tic
                fprintf('Correlating Data...\n');
                corr = photons_in_bins(ch1, ch2, lag_bin_edges, offset_lag);
                
                %%Normalization
                ch1_max = max(ch1); ch2_max = max(ch2);
                tag_range = max(ch1_max, ch2_max)- min(ch1_min, ch2_min); %range of tags in the entire dataset
                
                corr_div = corr./division_factor;
                corr_norm = 2.*corr_div.*tag_range./(tag_range-lags)...
                    .*ch1_max./(ch1_count.*ch2_count);
                
                %corr_norm = corr_norm(:, 1:(stop_bin));
                %lags = lags(start_bin:(stop_bin));
                
                fprintf('Done\n')
                toc
                
            end
            
            function [lag_bin_edges, lags, division_factor] = generate_log2_lags(t_end, coarseness)
                %%Generates log2 spaced photon bins
                cascade_end = floor(log2(t_end)); % cascades are collections of lags with equal bin spacing 2^cascade
                nper_cascade = coarseness; % number of equal
                
                n_edges = cascade_end*nper_cascade;
                lag_bin_edges = zeros(1, n_edges); %lag bins
                
                for j = 1:n_edges
                    if j == 1
                        lag_bin_edges(j) = 1;
                    else
                        lag_bin_edges(j) = lag_bin_edges(j-1)+2^(floor((j-1)./nper_cascade));
                    end
                end
                lags = diff(lag_bin_edges)./2+lag_bin_edges(1:end-1); %centers of the lag bins
                division_factor = kron(2.^(1:cascade_end), ones(1, nper_cascade)); %for normalization
            end
            
            
            function acf = photons_in_bins(ch1, ch2, lag_bin_edges, offset_lag)
                %%Counts the number of photons in the photon stream bins
                %%according to a prescription from Ted Laurence: Fast flexible
                %%algirthm for calculating photon correlation, Optics Letters,
                %%31,829,2006
                num_ch1 = numel(ch1);
                n_edges = numel(lag_bin_edges);
                low_inds = ones(1, n_edges-1);
                low_inds(1) = 2;
                max_inds = ones(1, n_edges-1);
                acf = zeros(1, n_edges-1);
                
                for phot_ind = 1:num_ch1
                    bin_edges = ch1(phot_ind)+lag_bin_edges+offset_lag;
                    
                    for k = 1:n_edges-1
                        while low_inds(k) <= numel(ch2) && ch2(low_inds(k)) < bin_edges(k)
                            low_inds(k) = low_inds(k)+1;
                        end
                        
                        while max_inds(k) <= numel(ch2) && ch2(max_inds(k)) <= bin_edges(k+1)
                            max_inds(k) = max_inds(k)+1;
                        end
                        
                        low_inds(k+1) = max_inds(k);
                        acf(k) = acf(k)+(max_inds(k)-low_inds(k));
                    end
                end
            end
            
        end
        function obj = photon_corr_manual_sync(obj, file_in_key,file_out_key,channels, time_bounds, lag_precision, lag_offset, use_dll, Sync_Rate)
            %% Contributed by Boris Spokoyny. Allow to correlate the photon-stream on a log timescale.
            % file_in key: file ID of the photons-file to correlate
            % file_out_key: ID for the struct array added to the photons
            % object containing the correlation.
            % channels: Hydraharp channels to be correlated. e.g. [0,1] for
            % cross-correlation of channels 0 and 1. 
            % time-bound: upper and lower limit for the correlation. In ps
            % for T2, in pulses for T3. 
            %lag_precision: Number of correlation points between
            %time-spacings of log(2).
            %lag_offset: offset in ps or pulses between the channels. 
            %use_dll for using a C++ code for the correlation algorithm for
            %improved efficiency. 
            
            dot_photons_handle = fopen(strcat(obj.pathstr,file_in_key,'.photons'),'r');
            
            switch obj.header.MeasurementMode
                case 2 %T2 mode
                    data = fread(dot_photons_handle, [2, obj.header.nRecords], 'uint64');
                case 3 %T3 mode
                    data = fread(dot_photons_handle, [3, obj.header.nRecords], 'uint64');
            end
            fclose(dot_photons_handle);
            
            %split inot channels
            ch0 = data(2, data(1, :) == channels(1)); %ch0 syncs
            ch1 = data(2, data(1, :) == channels(2)); % ch1 syncs
            
            start_time = time_bounds(1);
            stop_time = time_bounds(2);
            if use_dll
                [lags, corr_norm] = cross_corr_dll(ch0,ch1, start_time, stop_time,...
                    lag_precision, lag_offset);
                if obj.header.MeasurementMode == 3
                    sync_period = 1E12/Sync_Rate;      % in picoseconds
                    lags = lags.*sync_period;
                    [~, start_bin] = min(abs(lags-start_time.*sync_period));
                    [~, stop_bin] = min(abs(lags-stop_time.*sync_period));
                    lags = lags(start_bin:stop_bin);
                    corr_norm = corr_norm(start_bin:stop_bin);
                end
            else
                [lags, corr_norm] = cross_corr(ch0,ch1, start_time, stop_time, ...
                    lag_precision, lag_offset);
                if obj.header.MeasurementMode == 3
                    sync_period = 1E12/Sync_Rate;      % in picoseconds
                    lags = lags.*sync_period;
                end
            end
            
            if ~isprop(obj,file_out_key)
                dummy=addprop(obj,file_out_key)
            end
            obj.(file_out_key)=struct('lags',lags,'corr_norm',corr_norm);
            
            %helper functions.
            function [lags, corr_norm] = cross_corr(ch1,ch2, start_time, stop_time, coarseness, offset_lag)
                %this algorithm computes the cross-correlation between ch1 and ch2
                %variables.
                %-- For T2 data, ch1 and ch2 correspond to absolute photon arrival
                %   times.
                %-- For T3 data, ch1 and ch2 should correspond to the photon arrival
                %   sync pulse number.
                %-- start_time and stop_time for T2 data should be in time units of the
                %-- Photon arrival times, and for T3 data should be in units of sync pulses.
                %-- The correlation lags are log(2) spaced with coarseness # of lags per
                %   cascade, i.e. if start_time = 1; stop_time = 50; coarseness = 4; the
                %   lag bins will be [ 1,2,3,4 ; 6,8,10,12 ; 16,20,24,28; 36,44 ] if
                %   coarseness = 2, the lag bins become [ 1,2 ; 4,6 ; 10,14 ; 22,30 ; 46 ]
                %-- The algorithm works by taking each photon arrival time and counting the
                %   number of photons that are lag bins away. For example say our lag bin
                %   edges are [1, 2, 4, 6, 10, 14]
                
                %   Time Slot: 1   2   3   4   5   6   7   8   9   10  11  12  13  14  15
                %   Photon?    1   1   0   0   1   1   0   1   1   0   1   0   0   0   1
                %   1st Photon ^
                %   Lag bins   |   |       |       |               |               |
                %   # photons    1     0       2           2                1
                %   2nd Photon     ^
                %   Lag bins       |   |       |       |               |               |
                %   # photons        0     1       1           3               1
                %   3rd Photon                 ^
                %   Lag bins                   |   |       |       |               |               |
                %   # photons                    1     2       1           1               1
                
                %etc..
                % the cross-correlation is the sum of all the photons for each bin, i.e
                % for the three runs above we get [2, 3, 4, 6, 3].
                
                % generates lags bins with log2 spacing
                % Here stop time, refers to a lag time in picoseconds for t2 data or
                % maximum number of pulses for t3 data.
                [lag_bin_edges, lags, division_factor] = generate_log2_lags(stop_time, coarseness);
                
                %find the bin region
                [~, start_bin] = min(abs(lag_bin_edges-start_time));
                [~, stop_bin] = min(abs(lag_bin_edges-stop_time));
                lag_bin_edges = lag_bin_edges(start_bin:stop_bin); %bins
                lags = lags((start_bin+1):stop_bin); % centers of bins
                division_factor = division_factor((start_bin+1):stop_bin); %normalization factor
                
                %%counters etc for normalization
                ch1_min = inf; ch2_min = inf; %minimum time tag
                ch1_count = 0; ch2_count = 0; %photon tally in each channel
                num_ch1 = numel(ch1); num_ch2 = numel(ch2); %number of photons in each channel
                ch1_count = ch1_count+num_ch1; ch2_count = ch2_count+num_ch2;
                ch1_min = min(ch1_min, min(ch1)); ch2_min = min(ch2_min, min(ch2));
                
                %%Correlation
                tic
                fprintf('Correlating Data...\n');
                corr = photons_in_bins(ch1, ch2, lag_bin_edges, offset_lag);
                
                %%Normalization
                ch1_max = max(ch1); ch2_max = max(ch2);
                tag_range = max(ch1_max, ch2_max)- min(ch1_min, ch2_min); %range of tags in the entire dataset
                
                corr_div = corr./division_factor;
                corr_norm = 2.*corr_div.*tag_range./(tag_range-lags)...
                    .*ch1_max./(ch1_count.*ch2_count);
                
                %corr_norm = corr_norm(:, 1:(stop_bin));
                %lags = lags(start_bin:(stop_bin));
                
                fprintf('Done\n')
                toc
                
            end
            
            function [lag_bin_edges, lags, division_factor] = generate_log2_lags(t_end, coarseness)
                %%Generates log2 spaced photon bins
                cascade_end = floor(log2(t_end)); % cascades are collections of lags with equal bin spacing 2^cascade
                nper_cascade = coarseness; % number of equal
                
                n_edges = cascade_end*nper_cascade;
                lag_bin_edges = zeros(1, n_edges); %lag bins
                
                for j = 1:n_edges
                    if j == 1
                        lag_bin_edges(j) = 1;
                    else
                        lag_bin_edges(j) = lag_bin_edges(j-1)+2^(floor((j-1)./nper_cascade));
                    end
                end
                lags = diff(lag_bin_edges)./2+lag_bin_edges(1:end-1); %centers of the lag bins
                division_factor = kron(2.^(1:cascade_end), ones(1, nper_cascade)); %for normalization
            end
            
            
            function acf = photons_in_bins(ch1, ch2, lag_bin_edges, offset_lag)
                %%Counts the number of photons in the photon stream bins
                %%according to a prescription from Ted Laurence: Fast flexible
                %%algirthm for calculating photon correlation, Optics Letters,
                %%31,829,2006
                num_ch1 = numel(ch1);
                n_edges = numel(lag_bin_edges);
                low_inds = ones(1, n_edges-1);
                low_inds(1) = 2;
                max_inds = ones(1, n_edges-1);
                acf = zeros(1, n_edges-1);
                
                for phot_ind = 1:num_ch1
                    bin_edges = ch1(phot_ind)+lag_bin_edges+offset_lag;
                    
                    for k = 1:n_edges-1
                        while low_inds(k) <= numel(ch2) && ch2(low_inds(k)) < bin_edges(k)
                            low_inds(k) = low_inds(k)+1;
                        end
                        
                        while max_inds(k) <= numel(ch2) && ch2(max_inds(k)) <= bin_edges(k+1)
                            max_inds(k) = max_inds(k)+1;
                        end
                        
                        low_inds(k+1) = max_inds(k);
                        acf(k) = acf(k)+(max_inds(k)-low_inds(k));
                    end
                end
            end
            
        end
        function obj = get_blinking_analysis(obj,intensity_ID,boundaries)
            % boundaries are up to N=3 right now, so basically 2 or three
            % states. Also intensity compiling needs to have been run
            % before anything can be done with this function.
            if length(boundaries)==3
                
                k=length(obj.(intensity_ID).trace);
                
                count_bright=[];
                count_grey=[];
                count_dark=[];
                
                for i=1:k
                    cpb=obj.(intensity_ID).trace(i,1);
                    if (boundaries(2)<cpb) && (cpb<boundaries(3))
                        count_bright=[count_bright,i];
                    end
                    if (boundaries(1)<cpb) && (cpb<boundaries(2))
                        count_grey=[count_grey,i];
                    end
                    if cpb<boundaries(1)
                        count_dark=[count_dark,i];
                    end
                end
                
                %% now convert this to a vector of the length of the different durations.
                
                count_bright_length=[];
                running_count=1;
                m=(length(count_bright)-1);
                
                for i=1:m
                    a=count_bright(i)+1;
                    b=count_bright(i+1);
                    if a==b
                        running_count=running_count+1;
                    end
                    if a~=b || i==m
                        count_bright_length=[count_bright_length,running_count];
                        running_count=1;
                    end
                end
                
                count_grey_length=[];
                running_count=1;
                m=(length(count_grey)-1);
                
                for i=1:m
                    a=count_grey(i)+1;
                    b=count_grey(i+1);
                    if a==b
                        running_count=running_count+1;
                    end
                    if a~=b || i==m
                        count_grey_length=[count_grey_length,running_count];
                        running_count=1;
                    end
                end
                
                
                count_dark_length=[];
                running_count=1;
                m=(length(count_dark)-1);
                
                for i=1:m
                    a=count_dark(i)+1;
                    b=count_dark(i+1);
                    if a==b
                        running_count=running_count+1;
                    end
                    if a~=b || i==m
                        count_dark_length=[count_dark_length,running_count];
                        running_count=1;
                    end
                end
                
                
                if ~isprop(obj,'counts_length_intensity')
                    dummy=addprop(obj,'counts_length_intensity')
                end
                obj.counts_length_intensity=struct('count_dark_length',count_dark_length,'count_grey_length',count_grey_length,'count_bright_length',count_bright_length);
                
            end
            
            if length(boundaries)==2
               k=length(obj.(intensity_ID).trace);
                
                count_bright=[];
                count_dark=[];
                
                for i=1:k
                    cpb=obj.(intensity_ID).trace(i,1);
                    if (boundaries(1)<cpb) && (cpb<boundaries(2))
                        count_bright=[count_bright,i];
                    end
                    if cpb<boundaries(1)
                        count_dark=[count_dark,i];
                    end
                end
                
                %% now convert this to a vector of the length of the different durations.
                
                count_bright_length=[];
                running_count=1;
                m=(length(count_bright)-1);
                
                for i=1:m
                    a=count_bright(i)+1;
                    b=count_bright(i+1);
                    if a==b
                        running_count=running_count+1;
                    end
                    if a~=b || i==m
                        count_bright_length=[count_bright_length,running_count];
                        running_count=1;
                    end
                end
             
                
                count_dark_length=[];
                running_count=1;
                m=(length(count_dark)-1);
                
                for i=1:m
                    a=count_grey(i)+1;
                    b=count_grey(i+1);
                    if a==b
                        running_count=running_count+1;
                    end
                    if a~=b || i==m
                        count_dark_length=[count_dark_length,running_count];
                        running_count=1;
                    end
                end
                
                if ~isprop(obj,'counts_length_intensity')
                    dummy=addprop(obj,'counts_length_intensity')
                end
                obj.counts_length_intensity=struct('count_dark_length',count_dark_length,'count_bright_length',count_bright_length);
            end
        end
        
        %% high level function for analyzing derived properties of the photon stream (fitting FCS traces etc.)
        function obj=fit_auto_corr_FCS_trace(obj,file_key_in,p0,afterpulsing_time)
            
            if nargin<4
                afterpulsing_time=2E6;%given in ps.
            end
            
            %load the auto-correlation function defined by file_key_in
            data=csvread([obj.pathstr,file_key_in,'.intensity_corr']);
            tau=(data(:,3)+data(:,4))./2;
            
            index=[];
            [foo,index(1)]=min((tau-afterpulsing_time).^2); %only treat data after the detector deadtime.
            index(2)=min(find(data(:,5)==0)); %getting index of the first zero value in the correlation function.
            
            tau_select=tau(index(1):index(2)-1);%discarding the last point of tau.
            autocorrelation=data(index(1):index(2)-1,5);%discarding early (afterpulsing) and late (g(2)==0) data.
            
            % Define exponential function
            fh = @(x,p) p(1) + p(4).*(1./((1+(x./p(2))).*(1+(p(3).^(-2)).*x./p(2)).^0.5));
            
            % define the error function. this is the function to
            % minimize: you can also use norm or whatever:
            errfh = @(p,x,y) sum((y(:)-fh(x(:),p)).^2);
            
            % an initial guess of the exponential parameters
            if nargin<3
                p0 = [min(autocorrelation) 3E9 1 max(autocorrelation)];
            end
            
            % search for solution
            P = fminsearch(errfh,p0,[],tau_select,autocorrelation);
            residuals=(autocorrelation-fh(tau_select,P));
            fit_curve=fh(tau,P);%extrapolating the curve for all tau.
            tau_limits_for_fit=[min(tau) max(tau)];
            
            %writing the fit results as object properties.
            
            if ~isprop(obj,'FCS_fit')
                dummy=addprop(obj,'FCS_fit');
            end
            
            %wrapping results in a structure.
            obj.FCS_fit=struct('P',P,'residuals', residuals,'fit_curve',fit_curve,'tau_limits_for_fit',tau_limits_for_fit,'tau_for_fit',tau_select);
            
        end
        
        %% depreciated functions.
        function obj=int_file_subset(obj,file_in_key,file_out_key,percent_start,percent_end)
            %% Creates .photons_int files containing only a subset of all photon records
            % it is required that the .bin_2_int() method has already been
            % run so you have a record of the earliest and latest photon.
            % Input is the start and stop times by percent of total time of
            % experiment.
            % file_in_key is the filename of the .photons file without
            % extension. file_out_key is the filename of the created
            % .photons_int, except the end will be appended with the
            % percent of the total time of experiment that the file's
            % photons represent.
            
            time_start = obj.earliest_photon + (obj.latest_photon - obj.earliest_photon) * percent_start / 100;
            time_end = obj.earliest_photon + (obj.latest_photon - obj.earliest_photon) * percent_end / 100;
            buffer_size = obj.buffer_size;
            fid_in = fopen(strcat(obj.pathstr,file_in_key,'.photons'),'r');
            fid_out = fopen(strcat(obj.pathstr,file_out_key,'_',num2str(percent_start),'_',num2str(percent_end),'.photons_int'),'w');
            if obj.mode == 3
                while 1 %~while true
                    batch = fread(fid_in,[3,buffer_size],'uint64')';%reading in an array of the binary to process.
                    lbatch = length(batch);
                    for k = 1:lbatch
                        if (batch(k,3) > time_start) && (batch(k,3) < time_end)
                            fprintf(fid_out,'%2u,%9u,%6u\n',batch(k,:)');
                        end % if
                    end % for
                    if lbatch < buffer_size; %stop reading in batches when we have reached end of file.
                        break
                    end % if
                end % while
            end % if
            if obj.mode == 2
                while 1 %~while true
                    batch = fread(fid_in,[2,buffer_size],'uint64')';%reading in an array of the binary to process.
                    lbatch = length(batch);
                    for k = 1:lbatch
                        if (batch(k,2) > time_start) && (batch(k,2) < time_end)
                            fprintf(fid_out,'%2u,%9u\n',batch(k,:)');
                        end % if
                    end % for
                    if lbatch < obj.buffer_size; %stop reading in batches when we have reached end of file.
                        break
                    end % if
                end % while
                
                fclose(fid_in);
                fclose(fid_out);
            end % if
        end % function
        function obj=get_intensity_correlation(obj,file_key_in,file_key_out,mode,bin_width,time_scale,n_channels)
            %get the intensity correlation function from Toms code.The
            %bin_width is in pulses for t3 data and in ps for t2 data.
            
            %THIS ONLY WORKS FOR MAC OS! WINDOWS OR LINUX USERS HAVE TO
            %FIND A DIFFERENT WAY TO INTERFACE WITH TOMS CODE.

            command = ['photon_intensity_correlate.exe --file-in ',' ',obj.pathstr,file_key_in,'.photons_int', ...
                ' --file-out ',' ',obj.pathstr,file_key_out,'.intensity_corr ',' -m ', mode,' -g 2 -w ', bin_width,...
               ' --time-scale ',time_scale,' -c ', num2str(n_channels)];
            system(command)%calls the system command to get intensity_correlation_functions - only on MAC OS
        end
        function obj=get_photon_gn_mac(obj,file_key_in,file_key_out,channels,order,time,pulse)
            command=['photon_gn --file-in ',[obj.pathstr,file_key_in,'.photons_int'], ' --file-out ',[obj.pathstr,file_key_out,'.ht3 '],' --mode ', t3,' --order ', order, '--pulse', pulse, ' --time ', time,' --channels ', num2str(channels)]
            

            system(command)%calls the system command to get intensity_correlation_functions - only on MAC OS    

        end
        function obj=parse_g2(obj,input_folder,file_key_out,rrate,channels);
            %input_folder is the folder that Thomas's code creates without
            %any extension
            %rrate is repetition rate in Hz
            %Parses the g2 and intensity data structures (originally generated by Tom's
            %correlation code) and converts time units to reasonable values.  In this
            %case, we're parsing the results from a t3 run with time and pulse
            %difference.  However, we're going to assume that the run is
            %pulse-integrated.  That is, each intensity bin is the integral over an
            %entire pulse.  In this case, we clip out the time columns and convert from
            %pulse-units to time units

            %Takes folder from Thomas's code and picks out the g2 file
            file_name_corr = {input_folder};
            file_name_corr = strjoin(file_name_corr,'\');
            file_name_corr = {file_name_corr,'g2.run'};
            file_name_corr = strjoin(file_name_corr,'.');
            file_name_corr = {file_name_corr,'g2'};
            file_name_corr = strjoin(file_name_corr,'\');
            g2 = dlmread(file_name_corr);
            
            %Takes folder from Thomas's code and picks out the intensity
            file_name_corr = {input_folder};
            file_name_corr = strjoin(file_name_corr,'\');
            file_name_corr = {file_name_corr,'g2.run'};
            file_name_corr = strjoin(file_name_corr,'.');
            file_name_corr = {file_name_corr,'intensity'};
            file_name_corr = strjoin(file_name_corr,'\');
            intensity = dlmread(file_name_corr);
            
            %clip out time columns
%             g2 = g2(:,[1:4,7]);

            %timeconversion
            corr_t_units = 1000/rrate; %converts from pulse separation to milliseconds
            intensity_t_units = 1/(rrate); %converts from pulse number to seconds

            %PARSE g2
            %Find out how long the data is.  Each g2 structure contains the 00/01/10/11
            %correlation data.
            l = size(g2,1)/channels^2;

            %Extract the correlation data
            for ii=1:channels^2
                g2_split(:,ii)=g2((ii-1)*l+1:ii*l,5);
                tau_split(:,ii)=0.5*(g2((ii-1)*l+1:ii*l,3)+g2((ii-1)*l+1:ii*l,4))*corr_t_units;
              
            end 
                    
            %Extract a centered-bin time axis for each correlation measurement (they
            %should all be the same, but who knows what really happened, I guess...)
            
            %PARSE intensity
            intensity_t = 0.5*(intensity(:,1)+intensity(:,2))*intensity_t_units;
            
            for ii=1:channels
                intensity_split(:,ii)=intensity(:,ii+2)./((intensity(:,2)-intensity(:,1))*intensity_t_units);             
            end 
            
            intensity_total = sum(intensity_split,2);

            %normalize correlation functions; this will be based on all of the time
            %units being the same.
            totaltime = (intensity(end,1)-intensity(1,1))*corr_t_units; 
            for ii=3:channels+2
                AverageI(ii-2) = sum(intensity(:,ii))/totaltime;
            end
            
            ScalingFactor = prod(AverageI)*(g2(1,4)-g2(1,3))*corr_t_units*totaltime;
            

            count = 0;
            for ii = 1:channels
                for jj = 1:channels
                    count = count+1;
                    g2_norm(:,count) = g2_split(:,count)./(AverageI(ii)*AverageI(jj)*(g2(1:l,4)-g2(1:l,3))*corr_t_units*totaltime);
                end
            end

            %Averaging cross correlations to produce symmetric result
            taucross = tau_split(:,2);
            count = 0;
            g2cross = zeros(l,1);
            for ii = 1:channels
                for jj = 1:channels
                    count = count + 1;
                    if ii == jj
                        continue
                    end
                    g2cross(:) = g2cross(:)+g2_norm(:,count);
                end
            end
            g2cross = g2cross/(channels*(channels-1));


            %plot results for sanity check
            figure
            hold all
            for ii = 1:channels^2
                plot(tau_split(:,ii),g2_split(:,ii))
            end
            title('unnormalized correlation functions')
            
            
            figure
            hold all
            plot(intensity_t,intensity_split)
            title('intensity traces')

            figure
            hold all
            for ii = 1:channels^2
                plot(tau_split(:,ii),g2_norm(:,ii))
            end
            title('normalized correlation functions')
        
           
            %Script that performs the solution-g2 analysis.  It assumes you 
            %have run parse_g2 to populate the workspace.  Furthermore, it requires
            %that the data be a pulse-wise correlation function with symmetric,
            %linearly spaced bins.

            %Parameters
            PCH = 0;
            FCS_Allowoffset=0;
            Plot = 'True';

            %Fit the cross correlation to the FCS governing equations to
            %verify normalization and determine average occupation.
            if FCS_Allowoffset
                OffsetUB = 1.01;
                OffsetLB = .99;
            else
                OffsetUB = 1.00001;
                OffsetLB = .99999;
            end
            %FCS FIT OPTIONS
            fiteqn = 'a1./(1+b1*abs(x))+c1';
            coeffs = {'a1','b1','c1'};
            start = [0.3 0.8 1];
            upper = [10 1000 OffsetUB];
            lower = [0 .001 OffsetLB];
            %weights = 1./taucross((l+3)/2:l);
            weights = 1./taucross((l+11)/2:l);
            fo_ = fitoptions('method','NonlinearLeastSquares','Lower',lower,'Upper',upper,'MaxFunEvals',60000,'MaxIter',4000,'TolFun',1e-008,'TolX',1e-008,'Weights',weights);
            set(fo_,'Startpoint',start);
            ft_ = fittype(fiteqn,'dependent',{'y'},'independent',{'x'},'coefficients',coeffs);  
            %DO FCS FIT
            %FCSfitobject = fit(taucross((l+3)/2:l),g2cross((l+3)/2:l),ft_,fo_);
            FCSfitobject = fit(taucross((l+11)/2:l),g2cross((l+11)/2:l),ft_,fo_);
            FCSfit = FCSfitobject(taucross);
            FCSDiffusionTime = 1/FCSfitobject.b1; %characteristic translational diffusion time in msec
            FCSAvgOccupation = 1/FCSfitobject.a1;
            FCSNormalizationOffset = FCSfitobject.c1-1;
            FCSResidual = mean((g2cross((l+3)/2:l)-FCSfit((l+3)/2:l)).^2);
            %PHENOMENOLOGICAL FIT
            Phenomfitobject = fit(taucross((l+3)/2:(l+3)/2+99),g2cross((l+3)/2:(l+3)/2+99),'poly5');
            Phenomfit = Phenomfitobject(abs(taucross((l+1)/2-100:(l+1)/2+100)));
            PhenomAvgOccupation = 1/(Phenomfit(101)-1);

            if Plot
                %Plot fit
                figure
                subplot(3,3,1:6)
                semilogx(taucross((l+3)/2:l),g2cross((l+3)/2:l));
                hold all
                semilogx(taucross((l+3)/2:l),FCSfit((l+3)/2:l));
                hold off
                xlabel('tau (ms)')
                ylabel('FCS cross-correlation')
                title(strcat('AvgOccupation: ',num2str(FCSAvgOccupation),' Diffusion Time (ms): ',num2str(FCSDiffusionTime),'FCS Offset: ',num2str(FCSNormalizationOffset)))
                subplot(3,3,7:9)
                semilogx(taucross((l+3)/2:l),g2cross((l+3)/2:l)-FCSfit((l+3)/2:l));
                xlabel('tau (ms)')
                ylabel('Residual')

                %plot center results
                figure
                hold all
                plot(taucross((l+1)/2-100:(l+1)/2+100)*1e3,g2cross((l+1)/2-100:(l+1)/2+100));
                plot(taucross((l+1)/2-100:(l+1)/2+100)*1e3,FCSfit((l+1)/2-100:(l+1)/2+100));
                plot(taucross((l+1)/2-100:(l+1)/2+100)*1e3,Phenomfit);
                hold off
                xlabel('tau (us)')
                ylabel('Peak-integrated cross-correlation')
                title(strcat('Phenomenological AvgOccupation: ',num2str(PhenomAvgOccupation)))
            end

            %Calculate quantum yield ratios and stddev under two models (Uses scaling
            %factor variable now added to Parse_g2)
            CenterValue_Phenom = g2cross((l+1)/2);
            CenterStdDev_Phenom = sqrt(CenterValue_Phenom*ScalingFactor)/ScalingFactor;
            SideValue_Phenom = 1+1/PhenomAvgOccupation;
            SideStdDev_Phenom = sqrt(SideValue_Phenom*ScalingFactor)/ScalingFactor;

            CenterValue_FCS = g2cross((l+1)/2);
            CenterStdDev_FCS = sqrt(CenterValue_FCS*ScalingFactor)/ScalingFactor;
            SideValue_FCS = 1+1/FCSAvgOccupation;
            SideStdDev_FCS = sqrt(SideValue_FCS*ScalingFactor)/ScalingFactor;

            PhenomRatio = (CenterValue_Phenom-1)/(SideValue_Phenom-1);
            PhenomStdDev = PhenomRatio*sqrt((CenterStdDev_Phenom/(CenterValue_Phenom-1))^2 + (SideStdDev_Phenom/(SideValue_Phenom-1))^2);
            FCSRatio = (CenterValue_FCS-1)/(SideValue_FCS-1);
            FCSStdDev = FCSRatio*sqrt((CenterStdDev_FCS/(CenterValue_FCS-1))^2 + (SideStdDev_FCS/(SideValue_FCS-1))^2);

            if Plot
            %Report
            disp(strcat('Based on the FCS fit, the quantum yield ratio is: ',num2str(FCSRatio), '.  StdDev is: ',num2str(FCSStdDev)));
            disp(strcat('Based on the phenomenological fit, the quantum yield ratio is: ',num2str(PhenomRatio), '.  StdDev is: ',num2str(PhenomStdDev)));
            end

            if PCH
                %Check for aggregation by examining the intensity profile of the focal
                %volume.
                [nelements,centers] = hist(intensity_01, linspace(min(intensity_01),max(intensity_01),100));

                %Fit with Poissonian model
                fiteqn = 'c1*a1.^(b1*x)*exp(-a1)./gamma(b1*x+1)';
                coeffs = {'a1','b1','c1'};
                start = [10 .001 2000];
                fo_ = fitoptions('method','NonlinearLeastSquares','MaxFunEvals',60000,'MaxIter',4000,'TolFun',1e-008,'TolX',1e-008);
                set(fo_,'Startpoint',start);
                ft_ = fittype(fiteqn,'dependent',{'y'},'independent',{'x'},'coefficients',coeffs);
                %DO FCS FIT
                Intfitobject = fit(centers',nelements',ft_,fo_);
                figure
                hold on
                plot(centers,Intfitobject(centers))
                plot(centers,nelements)
            end
        end
        function obj=get_photon_gn_windows(obj,Data_Path,file_key_in,file_key_out,channels,order,time,pulse)
            % All input values need to be strings, not numbers. Data_Path 
            % is the folder where your data is stored. file key in
            % must be an ht3 file. File key out will determine the folder
            % name of the output file. Channels is the number of detectors
            % used (usually 2 or 4). Order is the correlation order. For
            % two photon events this is 2. Time is actually three numbers
            % separated by commas (no spaces). The first and last should be
            % opposites to create a symmetric function. The determine the
            % time away from the pulse which you extend looking for multi
            % photon events in picoseconds. The center number is the number
            % of time bins which that is divided into and must be an odd
            % number. (ie -400500,801,400500 -- goes +/- 400.5 ns from a 
            % pulse and divides that area into 801 bins). Pulse is similar
            % in structure to time. You need the pulse separation which you
            % extend in either direction. The way the code is written, it
            % looks for the pulse before the number you use, so you must
            % extend an extra 0.5 beyond the pulse you want to go to. To
            % just get center and side peaks you would use '-1.5,3,1.5'
            % which takes three pulses the minus tau, tau = 0 and plus tau
            % peaks. The center number should always be the sum of the
            % absolute value of the plus and minus values.
            command = ['"C:\cygwin64\bin\bash" -c "cd ' Data_Path '; C:/cygwin64/bin/picoquant --file-in ' file_key_in '| C:/cygwin64/bin/photon_gn --mode t3 --channels ' channels ...
                ' --order ' order ' --time ' time ' --pulse ' pulse ' --file-out ' file_key_out '"'];
            system(command)
        end         
    end
end
