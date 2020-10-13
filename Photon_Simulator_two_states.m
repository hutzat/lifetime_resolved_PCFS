%% Simulator for photon-files
% 
% % POISS{Pulse_Number}
% EXP{Lifetime}
% CENTER_COLOR=FUNCTION{SPECTRAL_DIFFUSION, } -> or whatever for spectral diffusion parameters. Could also be == zero.
% COHERENCE_TIME = FUNCTION{LINESHAPE} -> Coherence Time. Here, fluctuations in the coherence time can be introduced. If the coherence time for each photon is constant it means that we only have ONE pure dephasing process.
% 
% write this to a .photon_sim file

% this class assumes Poissonian Photon-Statistics and generates a
% stream_sim file that contains photon records that can then be read by the
% PhotonAnalyzer.m class to convert to .photons files.

function Photon_Simulator_two_states(avg_intensity,simulation_time, bin_width, Sync_rate, tau_1,tau_2, color1, color2, coherence1,coherence2, file_key_out)
%avg_intensity: photon count rate per second.
%bin_widt: serves as a sync divider. If >1, photons are simulated for
%multiple laser pulses.
%Sync_rate: rate of laser excitation in the experiment.
%tau_1: lifetime of state 1
%tau_2: lifetime of state 2
%color1: color of state 1 in nm.
%color2: color of state 2 in nm.
%coherence1: coherence state 1 in ps.
%coherence2: coherence state 2 in ps.
%file_key_out: file_ID for output file.

%convert the average intensity into an average bin_count.
avg_bin_count=avg_intensity/Sync_rate;

%bin_width is in pulses of Sync_rate frequency.
number_of_bins_to_simulate=simulation_time*Sync_rate/bin_width;

%instantiate a vector containing the numeration of laser pulses.
pulse_number_vector=[1:1:number_of_bins_to_simulate]';

% also open a file to write the photons to.
fid_out=fopen(strcat(file_key_out,'.stream_sim'),'w');

% instantiate an array to write the photon records to.
photon_records=zeros(number_of_bins_to_simulate,7);

% populate a collumn of the matrix with a possoinian distribution of
% photons.
photon_records(:,4)=pick_Poisson(avg_bin_count,number_of_bins_to_simulate);

% get the total pulse number for the photon records.
photon_records(:,2)=photon_records(:,4).*pulse_number_vector;

% this decides for each photon wehther it was emitted from state 1 or
% state2
foo=rand(number_of_bins_to_simulate,1);
photon_records(foo>=0.5,7)=1;

%look at state 1
indices = find(photon_records(:,7)==1);
photon_records(indices,3)=pick_Exp(tau_1,length(indices));
photon_records(indices,5)=ones(length(indices),1)*color1;
photon_records(indices,6)=ones(length(indices),1)*coherence1;

%look at state 2
indices = find(photon_records(:,7)==0);
photon_records(indices,3)=pick_Exp(tau_2,length(indices));
photon_records(indices,5)=ones(length(indices),1)*color2;
photon_records(indices,6)=ones(length(indices),1)*coherence2;

% reomve all photon records with zero and two photons per pulse.
indices = find(photon_records(:,4)==1); % find all the pulses that do not contain exactly one photon

% exclude irrelevant information from the records to be written to file.
photon_records_to_write=[photon_records(indices,1:3),photon_records(indices,5:6)];% only keep photon records that actualy contain one photon.

fwrite(fid_out,photon_records_to_write,'float64');

fclose(fid_out)

end

