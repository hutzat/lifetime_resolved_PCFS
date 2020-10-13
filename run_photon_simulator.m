%% Simulate the lifetime-resolved PCFS run for a Static Doublet
% This script takes the input parameters for the numerical PCFS simulation and
% creates photon-streams that contains the color, temporal distribution,
% and coherence properties of photons. It does that by looping over the
% interferometer stage positions and calling Photon_Simulator_two_states.m.

% Photon_Analyzer.m parses the .stream files of the raw photon emission
% events into .photon files containing the photon detection events. The
% necessary parameters are given by the PCFS experiments.

%PCFS_simul.m and t3PCFS_simul are PCFS-analysis classes that perform the
%PCFS-type analysis for the simulated data.

clearvars;

%% set parameters for the simulation
avg_intensity=1E5; % total photon counts per second.
simulation_time=15; % photon-stream simulation time at each stage position in sec.
bin_width=1; % set to 1 by default. The bin-width for a pulsed experiment is in units of pulses.
Sync_rate=2E7; % laser repetition rate.
tau1=100; % photoluminescence lifetime of state 1 in ps.
tau2=1000; % photoluminescence lifetime of state 2 in ps.
color1=500; % color of state 1 in wavelength [nm]
color2=500.0301% color of state 2 in wavelength [nm]
coherence1=60; % coherence time in ps.
coherence2=30; % coherence time in ps.
file_key_out='static_doublet'; % file_ID for the simulation.
buffer_size=10E6; % number of photons processed in memory. 
stage_pos=[-50E-3:1E-3:50E-3]; % center stage positions of the interferometer.

% rename the .pcfslog file if file_key_out has been changed by user.
newfile=strcat(file_key_out,'.pcfslog');
if strcmp(file_key_out,'static_doublet')==0
    movefile('static_doublet.pcfslog',newfile);
end

%% Create the .pos file.
% the .pos file contains the center stage positions of the run and is later
% parsed by the PCFS class. 
%The pcfs_log file is in the folder and does not

pos_name=strcat(file_key_out,'.pos')
csvwrite(pos_name,stage_pos')

%% loop over the stage position vector and create photon-stream (.stream_sim files) for each of the stage positions.
% this routine calls the photon simulator for each stage position and
% writes the binary output to file.

disp('SIMULATING PHOTON EMISSION ... THIS MAY TAKE UP TO HOURS...')
for i=0:(length(stage_pos)-1)
    disp(strcat('CURRENT STAGE POSITION:  ',num2str(i)))
    
    file_ID=strcat(file_key_out,num2str(i));
    Photon_Simulator_two_states(avg_intensity,simulation_time, bin_width,Sync_rate, tau1,tau2,color1,color2,coherence1,coherence2, file_ID);
    
end

%% Set the PCFS experimental parameters here.
dither_amplitude=1000; % dither amplitude in nm. Should be approx. a multiple of the emission color.
dither_time=15; % dither period in seconds.

% get a file list.
files=dir;


for i=3:numel(files)
    [path,f_name,ext]=fileparts(files(i).name);
    
    %read in the metadata from the pcfs logfile.
    if strcmp(ext,'.stream_sim')== true;
        
        foo=f_name((length(file_key_out)+1):end);
        k=str2num(foo);
        % Photon_Analyzer takes the .stream files, simulates the PCFS
        % experiment and writes the results to .photon files.
        Photon_Analyzer(f_name,f_name,stage_pos(k+1), dither_amplitude,dither_time,Sync_rate)
    end
end

%now that the .stream files have been parsed to .photon files, we delete
%the .stream files.
delete *.stream_sim
 
%% We're left with .photon files in the same form than retreived from experiments siliar to those in 
% Utzat et al. Science 2019.

%% We use the PCFS class to perform standard PCFS analysis PCFS analysis.
static_doublet_PCFS=PCFS_simul('/Users/hendrikutzat/Google Drive/2D_PCFS/static_doublet',1E6,3) % change path to the .photon files.
static_doublet_PCFS.get_sum_signal_all()
static_doublet_PCFS.get_intensity_correlations(3,[1,1E8],4,Sync_rate)
static_doublet_PCFS.get_blinking_corrected_PCFS_interferogram()
static_doublet_PCFS.get_spectral_correlation_simul('correlations',[0,0,0.01,0])
static_doublet_PCFS.static_doublet12.get_intensity('static_doublet12','int',1E6)

%% Surface plot of the spectral correlation.
figure()
surf(static_doublet_PCFS.correlations.zeta,static_doublet_PCFS.tau(1:70),static_doublet_PCFS.correlations.spectral_corr(1:70,:))
set(gca,'XScale','lin')
set(gca,'YScale','log')
xlim([-0.2,0.2])
caxis([0,1])
colormap hot
xlabel('\zeta [meV]')
ylabel('\tau ps[')
zlabel('p(\zeta) [a.u.]')
set(gca,'fontsize',18)
pbaspect([1,2,1])
title('Static Doublet - \tau_{coherence} = 40ps - \Delta = 50\mu eV')
colorbar

%% Contour plot of the spectral correlation.

figure()
contourf(static_doublet_PCFS.correlations.zeta,static_doublet_PCFS.tau(1:70),static_doublet_PCFS.correlations.spectral_corr(1:70,:))
set(gca,'XScale','lin')
set(gca,'YScale','log')
xlim([-0.2,0.2])
zlim([0,1])
colormap hot
pbaspect([1,2,1])
xlabel('\zeta [meV]')
ylabel('\tau [ps]')
zlabel('p(\zeta) [a.u.]')
set(gca,'fontsize',18)
pbaspect([1,2,1])
title('Static Doublet - \tau_{coherence} = 40ps - \Delta = 50\mu eV')
colorbar


figure()
contourf(static_doublet_PCFS.correlations.zeta,static_doublet_PCFS.tau(1:70),static_doublet_PCFS.correlations.spectral_corr(1:70,:))
set(gca,'XScale','lin')
set(gca,'YScale','log')
xlim([-0.2,0.2])
zlim([0,1])
colormap hot
pbaspect([1,2,1])
xlabel('\zeta [meV]')
ylabel('\tau [ps]')
zlabel('p(\zeta) [a.u.]')
set(gca,'fontsize',18)
pbaspect([1,2,1])
title('Static Doublet - \tau_{coherence} = 40ps - \Delta = 50\mu eV')
colorbar

%% tau plot turning around the x and y axis. plot 1 - '1D PCFS'
figure()
contourf(static_doublet_PCFS.tau(1:80),static_doublet_PCFS.correlations.zeta,static_doublet_PCFS.correlations.spectral_corr(1:80,:)')
set(gca,'YScale','lin')
set(gca,'XScale','log')
ylim([-0.2,0.2])
zlim([0,1])
colormap hot
pbaspect([1,2,1])
ylabel('\zeta [meV]')
xlabel('\tau [ps]')
zlabel('p(\zeta) [a.u.]')
set(gca,'fontsize',18)
pbaspect([2,1,1])
title('Static Doublet - 1D PCFS')
colorbar


%% tau plot turning around the x and y axis.
figure()
contourf(static_doublet_PCFS.tau(5:70),static_doublet_PCFS.PCFS_interferogram(1,:)*1000,static_doublet_PCFS.PCFS_interferogram(6:71,:)')
set(gca,'YScale','lin')
set(gca,'XScale','log')
%xlim([-0.2,0.2])
%zlim([0,1])
colormap hot
pbaspect([1,2,1])
ylabel('Stage Position [mm]')
xlabel('\tau [ps]')
zlabel('PCFS Interferogram [a.u.]')
set(gca,'fontsize',18)
pbaspect([1,2,1])
title('Static Doublet - \tau_{coherence} = 40ps - \Delta = 50\mu eV')
colorbar

figure()
plot(static_doublet_PCFS.PCFS_interferogram(1,:)*1000,static_doublet_PCFS.PCFS_interferogram(30,:),'-o','linewidth',3)
xlabel('Interferometer Stage Position [mm]')
ylabel('PCFS Interferogram')
set(gca,'fontsize',14)
ylim([-0.6,0.1])
title('PCFS Interferogram - \tau_{PCFS} = 1\mu s')

figure()
plot(static_doublet_PCFS.correlations.zeta,static_doublet_PCFS.correlations.spectral_corr(10,:),'-o','linewidth',3)
ylabel('p(\zeta) [a.u.]')
xlabel('\zeta [meV]')
set(gca,'fontsize',14)
ylim([0,1.1])
title('PCFS Spectral Correlation - \tau_{PCFS} = 1\mu s')
xlim([-0.4,0.4])

figure()
semilogx(static_doublet_PCFS.tau(1:70),2*sqrt(static_doublet_PCFS.correlations.fit_params(1:70,3)),'-o','linewidth',3)
ylabel('Lorentzian Fit FWHM [meV]')
xlabel('\tau [ps]')
set(gca,'fontsize',14)
title('Static Doublet FWHM')


%% create instance of t3PCFS.
static_doublet_t3_PCFS=t3PCFS_simul('/Users/hendrikutzat/Google Drive/2D_PCFS/static_doublet/',1E6)

% define the relative microtime windows along the T-direction.
tau_one=[1,100];
tau_two=[101,200];
tau_three=[201,300]
tau_four=[301,400]
tau_five=[401,2000];
tau_six=[2001,7100];
tau_all=[tau_one;tau_two;tau_three;tau_four;tau_five;tau_six];
tau_id_all={'one','two','three','four','five','six'}

static_doublet_t3_PCFS.parse_photons(tau_id_all,tau_all)

%create new PCFS children for each of the sorted streams.
static_doublet_t3_PCFS.create_single_PCFS_child('one',[1,10E7],4,3,Sync_rate) 
static_doublet_t3_PCFS.create_single_PCFS_child('two',[1,10E7],4,3,Sync_rate) 
static_doublet_t3_PCFS.create_single_PCFS_child('three',[1,10E7],4,3,Sync_rate) 
static_doublet_t3_PCFS.create_single_PCFS_child('four',[1,10E7],4,3,Sync_rate) 
static_doublet_t3_PCFS.create_single_PCFS_child('five',[1,10E7],4,3,Sync_rate) 
static_doublet_t3_PCFS.create_single_PCFS_child('six',[1,10E7],4,3,Sync_rate) 



%% Start of the lifetime-resolved PCFS analysis.

% look at the lifetime at one random stage position.
static_doublet_t3_PCFS.static_doublet32.get_intensity('static_doublet32','int',2E6)
static_doublet_t3_PCFS.static_doublet32.lifetime_histo('static_doublet32','lifetime',16)
gg=fit(static_doublet_t3_PCFS.static_doublet32.lifetime.time(2:6E3)',static_doublet_t3_PCFS.static_doublet32.lifetime.lifetime(2:6E3)','exp2','Start',[1E5,-1/100,1E5,-1/1000])

% plot the lifetime.
figure()
semilogy(static_doublet_t3_PCFS.static_doublet32.lifetime.time,static_doublet_t3_PCFS.static_doublet32.lifetime.lifetime,'o')%'linewidth',3)
hold on
semilogy(static_doublet_t3_PCFS.static_doublet32.lifetime.time,gg.a.*exp(gg.b.*static_doublet_t3_PCFS.static_doublet32.lifetime.time))
semilogy(static_doublet_t3_PCFS.static_doublet32.lifetime.time,gg.c.*exp(gg.d.*static_doublet_t3_PCFS.static_doublet32.lifetime.time))
set(gca,'fontsize',24)
xlim([-100,8000])
ylim([10,1E6])
xlabel('T [ps]')
ylabel('counts')

figure()
semilogy(static_doublet_t3_PCFS.static_doublet32.lifetime.lifetime(2:1E4),static_doublet_t3_PCFS.static_doublet32.lifetime.time(2:1E4),'-o')%'linewidth',3)
hold on
semilogy(gg.a.*exp(gg.b.*static_doublet_t3_PCFS.static_doublet32.lifetime.time(2:1E4)),static_doublet_t3_PCFS.static_doublet32.lifetime.time(2:1E4))
semilogy(gg.c.*exp(gg.d.*static_doublet_t3_PCFS.static_doublet32.lifetime.time(2:1E4)),static_doublet_t3_PCFS.static_doublet32.lifetime.time(2:1E4))
semilogy(gg(static_doublet_t3_PCFS.static_doublet32.lifetime.time(2:1E4)),static_doublet_t3_PCFS.static_doublet32.lifetime.time(2:1E4))
set(gca,'fontsize',24)
xlim([0,2E4])
ylim([1,8E3])
xlabel('counts')
ylabel('T [ps]')
legend('PL lifetime','bright','dark','all states')

%% calculate the ratio between the different states that are emissive.
time=[1:1:8000];
%ratios between the two states
stateA=gg.c.*exp(gg.d.*time)
stateB=gg.a.*exp(gg.b.*time)

probA=(stateA./(stateA+stateB)).^2;
probB=(stateB./(stateA+stateB)).^2;
probAB=(stateA./(stateA+stateB)).*(stateB./(stateA+stateB));

figure()
semilogx(time,probA,'linewidth',3)
hold on
semilogx(time,probB,'linewidth',3)
semilogx(time,2*probAB,'linewidth',3)
xlim([1,8000])
%set ( gca, 'ydir', 'reverse' )
legend('p_{AA}','p_{BB}','2p_{AB}')
set(gca,'fontsize',24)
ylabel('p(T)')
xlabel('T [ps]')

%% extract the relative weights for the 6 T-resolved PCFS slices.

%contains the boundaries for the tau-values.
tau_all;

AreaA=[];
AreaB=[];

for i=1:6
AreaA(i)=sum(stateA(tau_all(i,1):tau_all(i,2)))
AreaB(i)=sum(stateB(tau_all(i,1):tau_all(i,2)))
end

pAA=(AreaA./(AreaA+AreaB)).^2;
pBB=(AreaB./(AreaA+AreaB)).^2;
pA=(AreaA./(AreaA+AreaB));
pB=(AreaB./(AreaA+AreaB));
probabilities=[pA',pB'];

% basically sA(omega) and sB(omega) are unknown. Gamma1, Gamma2 and Omega
% are fitting parameters.


figure()
plot(static_doublet_t3_PCFS.one.PCFS_interferogram(50,:)/max(static_doublet_t3_PCFS.one.PCFS_interferogram(50,:)))
hold on
%plot(static_doublet_t3_PCFS.two.PCFS_interferogram(50,:)/max(static_doublet_t3_PCFS.two.PCFS_interferogram(50,:)))
%plot(static_doublet_t3_PCFS.three.PCFS_interferogram(50,:)/max(static_doublet_t3_PCFS.three.PCFS_interferogram(50,:)))
%plot(static_doublet_t3_PCFS.four.PCFS_interferogram(50,:)/max(static_doublet_t3_PCFS.four.PCFS_interferogram(50,:)))
%plot(static_doublet_t3_PCFS.five.PCFS_interferogram(50,:)/max(static_doublet_t3_PCFS.five.PCFS_interferogram(50,:)))
plot(static_doublet_t3_PCFS.six.PCFS_interferogram(50,:)/max(static_doublet_t3_PCFS.six.PCFS_interferogram(50,:)))
legend('time interval one','time interval three','time intervalfour','time intervalfive','time intervalsix')
xlabel('stage pos index')
ylabel('PCFS interferogram')

%% plot the spectral correlation for the different PCFS children.

static_doublet_t3_PCFS.one.get_blinking_corrected_PCFS_interferogram()
static_doublet_t3_PCFS.one.get_spectral_correlation_simul('correlations',[0,0,0.01,0])

static_doublet_t3_PCFS.two.get_blinking_corrected_PCFS_interferogram()
static_doublet_t3_PCFS.two.get_spectral_correlation_simul('correlations',[0,0,0.01,0])

static_doublet_t3_PCFS.three.get_blinking_corrected_PCFS_interferogram()
static_doublet_t3_PCFS.three.get_spectral_correlation_simul('correlations',[0,0,0.01,0])

static_doublet_t3_PCFS.four.get_blinking_corrected_PCFS_interferogram()
static_doublet_t3_PCFS.four.get_spectral_correlation_simul('correlations',[0,0,0.01,0])

static_doublet_t3_PCFS.five.get_blinking_corrected_PCFS_interferogram()
static_doublet_t3_PCFS.five.get_spectral_correlation_simul('correlations',[0,0,0.01,0])

static_doublet_t3_PCFS.six.get_blinking_corrected_PCFS_interferogram()
static_doublet_t3_PCFS.six.get_spectral_correlation_simul('correlations',[0,0,0.01,0])

%% get the FWHM for all of these and write them in a matrix.
% This analysis is akin to what is shown in Figure 3 of the manuscript.

two_D_FWHM=zeros(7,length(static_doublet_t3_PCFS.five.tau));
two_D_FWHM(1,:)=static_doublet_t3_PCFS.one.correlations.fit_params(:,3)';
two_D_FWHM(2,:)=static_doublet_t3_PCFS.two.correlations.fit_params(:,3)';
two_D_FWHM(3,:)=static_doublet_t3_PCFS.three.correlations.fit_params(:,3)';
two_D_FWHM(4,:)=static_doublet_t3_PCFS.four.correlations.fit_params(:,3)';
two_D_FWHM(5,:)=static_doublet_t3_PCFS.five.correlations.fit_params(:,3)';
two_D_FWHM(6,:)=static_doublet_t3_PCFS.six.correlations.fit_params(:,3)';
two_D_FWHM(7,:)=static_doublet_t3_PCFS.six.correlations.fit_params(:,3)';

tau_all

tau_plot_2d=[1,100,200,300,400,2000,7000]

figure()
contourf(static_doublet_t3_PCFS.one.tau(5:70),tau_plot_2d,real(2*1000*sqrt(two_D_FWHM(:,5:70))),40)
set(gca,'YScale','log')
%set ( gca, 'ydir', 'reverse' )
set(gca,'XScale','log')
colormap hot
%caxis([0,400])
colorbar
set(gca,'fontsize',18)
ylabel('T [ps]')
xlabel('\tau [ps]')
title('2D PCFS FWHM of two states with different lifetime, spectral diffusion and coherence')

figure()
surf(tau_plot_2d,static_doublet_t3_PCFS.one.tau(5:75),real(2*1000*sqrt(two_D_FWHM(:,5:75)))')
caxis([0,200])
set(gca,'YScale','log')
set(gca,'XScale','log')
ylim([1E5,1.5E11])
colormap hot
colorbar
pbaspect([1,2,1])
set(gca,'fontsize',18)
ylabel('\tau [ps]')
xlabel('T [ps]')
title('2D PCFS - \tau_{c}1=40ps, \tau_{c}(2)=10ps, \tau_{pop}(1)=100ps, \tau_{pop}(2)=1000ps')
hold on

%% also now look at slices and two different points in 2D map for the spectral correlation.
% for the static singlet these should be tau-invariant.
start_color=[0,0,0.5];

figure()
plot(static_doublet_t3_PCFS.one.correlations.zeta,static_doublet_t3_PCFS.one.correlations.spectral_corr(5,:),'Color',start_color +1/14,'linewidth',3)
hold on
plot(static_doublet_t3_PCFS.one.correlations.zeta,static_doublet_t3_PCFS.one.correlations.spectral_corr(30,:),'Color',start_color +2/14,'linewidth',3)
plot(static_doublet_t3_PCFS.one.correlations.zeta,static_doublet_t3_PCFS.one.correlations.spectral_corr(38,:),'Color',start_color +3/14,'linewidth',3)
plot(static_doublet_t3_PCFS.one.correlations.zeta,static_doublet_t3_PCFS.one.correlations.spectral_corr(40,:),'Color',start_color +4/14,'linewidth',3)
plot(static_doublet_t3_PCFS.one.correlations.zeta,static_doublet_t3_PCFS.one.correlations.spectral_corr(42,:),'Color',start_color +5/14,'linewidth',3)
plot(static_doublet_t3_PCFS.one.correlations.zeta,static_doublet_t3_PCFS.one.correlations.spectral_corr(44,:),'Color',start_color +6/14,'linewidth',3)
plot(static_doublet_t3_PCFS.one.correlations.zeta,static_doublet_t3_PCFS.one.correlations.spectral_corr(50,:),'Color',start_color +7/14,'linewidth',3)
xlim([-0.6,0.61])
ylim([-0.1,1.1])
xlabel('\zeta [meV]')
ylabel('p(\zeta)')
set(gca,'fontsize',24)

% plot the spectral correlation for different T.
figure()
subplot(1,2,1)
plot(static_doublet_t3_PCFS.one.correlations.zeta,static_doublet_t3_PCFS.one.correlations.spectral_corr(70,:),'red','linewidth',3)
xlim([-0.4,0.4])
ylim([-0.1,1.1])
subplot(1,2,2)
plot(static_doublet_t3_PCFS.six.correlations.zeta,static_doublet_t3_PCFS.six.correlations.spectral_corr(70,:),'red','linewidth',3)
xlim([-0.4,0.4])
ylim([-0.1,1.1])
xlabel('\zeta [meV]')
set(gca,'fontsize',18)
ylabel('p(\zeta)')

%% representation of the different components of the spectral correlation.

figure()
plot(static_doublet_t3_PCFS.one.correlations.zeta,mean(static_doublet_t3_PCFS.one.correlations.spectral_corr(15:30,:),1),'Color',start_color +1/14,'linewidth',3)
xlim([-0.4,0.4])
ylim([-0.1,1.1])
xlabel('\zeta [meV]')
ylabel('p(\zeta)')
set(gca,'fontsize',24)

figure()
plot(static_doublet_t3_PCFS.three.correlations.zeta,mean(static_doublet_t3_PCFS.three.correlations.spectral_corr(15:30,:),1),'Color',start_color +1/14,'linewidth',3)
xlim([-0.4,0.4])
ylim([-0.1,1.1])
xlabel('\zeta [meV]')
ylabel('p(\zeta)')
set(gca,'fontsize',24)


figure()
plot(static_doublet_t3_PCFS.six.correlations.zeta,mean(static_doublet_t3_PCFS.six.correlations.spectral_corr(15:30,:),1),'Color',start_color +1/14,'linewidth',3)
xlim([-0.4,0.4])
ylim([-0.1,1.1])
xlabel('\zeta [meV]')
ylabel('p(\zeta)')
set(gca,'fontsize',24)




%% perform some global fitting of the 2D PCFS matrix.
% fit each T with the auto-correlation of the linear combination of spectra
% with a given lineshape and linewidth as well as energy separation.
% So that is a three parameter fit without spectral diffusion.


spec_corr_slice=zeros(7,length(static_doublet_t3_PCFS.one.correlations.zeta));
spec_corr_slice(1,:)=static_doublet_t3_PCFS.one.correlations.spectral_corr(30,:)';
spec_corr_slice(2,:)=static_doublet_t3_PCFS.two.correlations.spectral_corr(30,:)';;
spec_corr_slice(3,:)=static_doublet_t3_PCFS.three.correlations.spectral_corr(30,:)';;
spec_corr_slice(4,:)=static_doublet_t3_PCFS.four.correlations.spectral_corr(30,:)';;
spec_corr_slice(5,:)=static_doublet_t3_PCFS.five.correlations.spectral_corr(30,:)';;
spec_corr_slice(6,:)=static_doublet_t3_PCFS.six.correlations.spectral_corr(30,:)';;
spec_corr_slice(7,:)=static_doublet_t3_PCFS.six.correlations.spectral_corr(30,:)';;

% Paper-friendly representation
figure()
surf(static_doublet_t3_PCFS.one.correlations.zeta,tau_plot_2d,spec_corr_slice)
set(gca,'XScale','lin')
set(gca,'YScale','log')
xlim([-0.2,0.2])
caxis([0,1])
colormap hot
xlabel('\zeta [meV]')
ylabel('\tau ps[')
zlabel('p(\zeta) [a.u.]')
set(gca,'fontsize',18)
pbaspect([1,2,1])
title('Static Doublet - \tau_{coherence} = 40ps - \Delta = 50\mu eV')
colorbar


% Paper-friendly representation
figure()
contourf(static_doublet_t3_PCFS.one.correlations.zeta,tau_plot_2d,spec_corr_slice)
set(gca,'XScale','lin')
set(gca,'YScale','log')
xlim([-0.2,0.2])
caxis([0,1])
colormap hot
xlabel('\zeta [meV]')
ylabel('T [ps]')
zlabel('p(\zeta) [a.u.]')
set(gca,'fontsize',18)
set ( gca, 'ydir', 'reverse' )
pbaspect([1,2,1])
title('Relaxation-resolved PCFS')
colorbar

%% plot 2 for the paper. 
figure()
contourf(tau_plot_2d,static_doublet_t3_PCFS.one.correlations.zeta,spec_corr_slice')
set(gca,'YScale','lin')
set(gca,'XScale','log')
ylim([-0.2,0.2])
%xlim([1,8000])
caxis([0,1])
colormap hot
ylabel('\zeta [meV]')
xlabel('T [ps]')
zlabel('p(\zeta) [a.u.]')
set(gca,'fontsize',18)
pbaspect([2,1,1])
title('Relaxation-resolved PCFS')
colorbar

figure()
semilogy(2*sqrt(two_D_FWHM(:,50)),tau_plot_2d,'-o','linewidth',3)
xlabel('FWHM [meV]')
ylabel('\tau [ps]')
set(gca,'fontsize',18)

%% look at slices along the x axis.


spec_corr_slice_x=zeros(7,length(static_doublet_t3_PCFS.one.correlations.zeta));
spec_corr_slice_x(1,:)=static_doublet_t3_PCFS.one.correlations.spectral_corr(50,:)';
spec_corr_slice_x(2,:)=static_doublet_t3_PCFS.two.correlations.spectral_corr(50,:)';;
spec_corr_slice_x(3,:)=static_doublet_t3_PCFS.three.correlations.spectral_corr(50,:)';;
spec_corr_slice_x(4,:)=static_doublet_t3_PCFS.four.correlations.spectral_corr(50,:)';;
spec_corr_slice_x(5,:)=static_doublet_t3_PCFS.five.correlations.spectral_corr(50,:)';;
spec_corr_slice_x(6,:)=static_doublet_t3_PCFS.six.correlations.spectral_corr(50,:)';;
spec_corr_slice_x(7,:)=static_doublet_t3_PCFS.six.correlations.spectral_corr(50,:)';;


figure()
contourf(static_doublet_t3_PCFS.one.tau(1:70),static_doublet_t3_PCFS.one.correlations.zeta,static_doublet_t3_PCFS.one.correlations.spectral_corr(1:70,:)')
set(gca,'XScale','log')
set(gca,'YScale','lin')
ylim([-0.5,0.5])
caxis([0,1])
colormap hot
ylabel('\zeta [meV]')
xlabel('\tau [ps]')
zlabel('p(\zeta) [a.u.]')
set(gca,'fontsize',18)
pbaspect([2,1,1])
title('Spectral Correlation p(\zeta,\tau, T=0-100ps)')
colorbar


figure()
contourf(static_doublet_t3_PCFS.three.tau(1:70),static_doublet_t3_PCFS.three.correlations.zeta,static_doublet_t3_PCFS.three.correlations.spectral_corr(1:70,:)')
set(gca,'XScale','log')
set(gca,'YScale','lin')
ylim([-0.5,0.5])
caxis([0,1])
colormap hot
ylabel('\zeta [meV]')
xlabel('\tau [ps]')
zlabel('p(\zeta) [a.u.]')
set(gca,'fontsize',18)
pbaspect([2,1,1])
title('Spectral Correlation p(\zeta,\tau, T=300-400ps)')
colorbar


figure()
contourf(static_doublet_t3_PCFS.six.tau(1:70),static_doublet_t3_PCFS.six.correlations.zeta,static_doublet_t3_PCFS.six.correlations.spectral_corr(1:70,:)')
set(gca,'XScale','log')
set(gca,'YScale','lin')
ylim([-0.5,0.5])
caxis([0,1])
colormap hot
ylabel('\zeta [meV]')
xlabel('\tau [ps]')
zlabel('p(\zeta) [a.u.]')
set(gca,'fontsize',18)
pbaspect([2,1,1])
title('Spectral Correlation p(\zeta,\tau, T=2000-7000ps)')
colorbar


figure()
subplot(1,2,1)
plot(static_doublet_t3_PCFS.three.correlations.zeta,static_doublet_t3_PCFS.three.correlations.spectral_corr(27,:),'red','linewidth',3)
xlim([-0.4,0.4])
ylim([-0.1,1.1])
xlabel('\zeta [meV]')
ylabel('p(\zeta)')
set(gca,'fontsize',18)
subplot(1,2,2)
plot(static_doublet_t3_PCFS.three.correlations.zeta,static_doublet_t3_PCFS.three.correlations.spectral_corr(70,:),'red','linewidth',3)
xlim([-1,1])
ylim([-0.1,1.5])
xlabel('\zeta [meV]')
set(gca,'fontsize',18)
ylabel('p(\zeta)')

%% Plotting the interferogram

start_color=[0,0,0.5];
positions=[20,45,48,49,50,51];

%100ps
figure()
for ii=1:6
semilogx(static_doublet_t3_PCFS.one.tau,static_doublet_t3_PCFS.one.cross_corr_interferogram(2:end,positions(ii)),'-o','Color', start_color +ii/12,'linewidth',3,'MarkerSize',7)
hold on
end
xlim([2E5,1E11])
ylabel('g_{x}^{(2)}')
xlabel('\tau [ps]')
ylim([0.3,1.5])
set(gca,'fontsize',18)

% 1000ps
figure()
for ii=1:6
semilogx(static_doublet_t3_PCFS.six.tau,static_doublet_t3_PCFS.six.cross_corr_interferogram(2:end,positions(ii)),'-o','Color', start_color +ii/12,'linewidth',3,'MarkerSize',7)
hold on
end
xlim([2E5,1E11])
ylim([0.3,1.5])
ylabel('g_{x}^{(2)}')
xlabel('\tau [ps]')
set(gca,'fontsize',18)

% PCFS interferogram at 100ps and 400ps and 1000ps
positions=[]
c=physconst('LightSpeed')

figure()
plot(static_doublet_t3_PCFS.one.stage_positions/c*1E12,mean(static_doublet_t3_PCFS.one.PCFS_interferogram(5:35,:),1),'o-','linewidth',3,'Color',start_color+1/8)
hold on
%plot(static_doublet_t3_PCFS.two.stage_positions/c*1E12,mean(static_doublet_t3_PCFS.two.PCFS_interferogram(5:35,:),1),'o-','linewidth',3,'Color',start_color+2/8)
plot(static_doublet_t3_PCFS.three.stage_positions/c*1E12,mean(static_doublet_t3_PCFS.three.PCFS_interferogram(5:35,:),1),'o-','linewidth',3,'Color',start_color+3/8)
plot(static_doublet_t3_PCFS.six.stage_positions/c*1E12,mean(static_doublet_t3_PCFS.six.PCFS_interferogram(5:35,:),1),'o-','linewidth',3,'Color',start_color+4/8)
xlim([0,100])
%ylim([0.3,1.5])
ylabel('G^{(2)}(\delta)')
xlabel('time [ps]')
set(gca,'fontsize',18)
legend('T=0-100ps','T=200-300ps','T=2000ps-7000ps')

