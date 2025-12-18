function p = get_cdf_val(x, spike_num)
%GET_CDF_VAL - generate the cdf value at position x
%
% P = GET_CDF_VAL(X, SPIKE_NUM)
%
% Inputs:
%  X - where the cdf value will be evaluated at
%  SPIKE_NUM - the number of spikes in the original json file
%
% OUTPUTS:
%  P - cumulative distribution function that is evaluated at position x
p = normcdf(x,0,sqrt(0.5)/sqrt(spike_num));

end

