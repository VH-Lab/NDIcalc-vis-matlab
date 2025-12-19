function sigma = get_sigma(spike_num)
%GET_SIGMA - calculate the standard deviation for the cumulative distribution function
%
% SIGMA = vis.revcorr.get_sigma(SPIKE_NUM)
%
% Inputs:
%  SPIKE_NUM - the number of spikes in the original json file
%
% OUTPUTS:
%  SIGMA - the standard deviation for the cumulative distribution function

sigma = sqrt(0.5)*sqrt(spike_num);
end

