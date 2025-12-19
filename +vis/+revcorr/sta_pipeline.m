function [sta,p_val, rescale, cmap] = sta_pipeline(s,kx_v, ky_v, frameTimes, spiketimes, T_coords, X_coords, Y_coords)
%STA_PIPELINE - generate the sta, p-value array, rescaled matrix
%for plotting and corresponding color map. 
%
% [STA,P_VAL, RESCALE, CMAP] = vis.revcorr.sta_pipeline(S,KX_V, KY_V, FRAMETIMES, SPIKETIMES, T_COORDS, X_COORDS, Y_COORDS)
%
% Inputs:
% s - an array of the values of S for each stimulus
% kx_v - an array of the values of kx for each stimulus
% ky_v - an array of the values of ky for each stimulus
% frameTimes - an array of the time point for each frame
% spiketimes - an array of the time point for each spike
% T_coords - an array of reconstruction time
% X_coords - the size of the Hartley stimulus 
% Y_coords - the size of the Hartley stimulus 

% OUTPUTS:
%  STA - an X_COORDS x X_COORDS x TMAX array, where TMAX = (T_COORDS(end) - T_COORDS(1))/mean(diff(T_COORDS) + 1, the number of timesteps
%    representing the spike-triggered average
%  P_VAL - an array with the same size as STA array, and record the
%   significance of each pixel
%  RESCALE - an array with the same size as STA array. It takes the value
%   of P_VAL array and rescale it for plotting purposes. 
%  CMAP - a 256x3 array representing the color map

deltaT = mean(diff(T_coords));
rf_range = T_coords(end) - T_coords(1);
tmax = round(rf_range/deltaT) + 1;
M = size(X_coords(:), 1);
spike_num = size(spiketimes, 1);
sta = vis.revcorr.generate_STA(s,kx_v, ky_v, frameTimes, spiketimes, rf_range, T_coords, tmax, M);
[p_val_adjusted, p_val, rescale] = vis.revcorr.calc_significance(sta,spike_num);
cmap = vis.revcorr.get_cmap();
end

