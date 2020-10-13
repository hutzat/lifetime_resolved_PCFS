function photon_number_bin=pick_Poisson(avg_intensity,vector_length)
% function returns a number of photons for a given bin and for a given
% call. 
% intesity is in counts per second.
% bin_width is in pulses.


photon_number_bin=poissrnd(avg_intensity,vector_length,1);



end
