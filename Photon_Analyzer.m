function Photon_Analyzer(file_key_in,file_key_out,stage_position,dither_amplitude,dither_time_period,Sync_Rate)
%file_key_in: The filename to be read.
%file_key_out: The file_ID for the .photons output file.
%stage_position: the central stage position of the photon-stream.
%dither_amplitude: amplitude of the dither movement in nm.
%dither_time_period: time-period of dither in seconds.
%Sync_Rate: laser repetition rate in Hz.

c = physconst('LightSpeed');


read_file=strcat(file_key_in,'.stream_sim')

fid=fopen(read_file,'r');
fid_out=fopen(strcat(file_key_out,'.photons'),'w');

data=fread(fid,'float64');

photons=reshape(data,[],5);

time_vector=photons(:,2)/Sync_Rate+photons(:,3)/1E12;

% convert the center stage position and dither parameter into a time-delay.
stage_position_time= stage_position/c + dither_amplitude*1E-9/(2*c)*sawtooth(2*pi/dither_time_period*time_vector);

%% calculate the probability of the photon being detected with Detector 0.
% this is done taking the stage_position_time, color and coherence time of
% the photon into account.
p_0 = 0.5 + 0.5*sin(2*pi*stage_position_time./(photons(:,4)*1E-9/c)).*exp(-abs(stage_position_time)./(photons(:,5)/1E12));

% generate a random number to draw the outcomes of photo-detection with
% detector 0 with p_0 probability.
randomnumber=rand(length(time_vector),1);
detectors=photons(:,1);
detectors(randomnumber>=p_0)=1;

% write these photons to file.
photons_to_write=photons;
photons_to_write(:,1)=detectors;
photons_to_write=photons_to_write(:,1:3);

% write .photons file.
 fwrite(fid_out,photons_to_write','ubit64');

fclose(fid)
fclose(fid_out)

end