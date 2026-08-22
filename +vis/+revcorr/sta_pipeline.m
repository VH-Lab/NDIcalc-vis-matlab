function [sta,p_val, rescale, cmap] = sta_pipeline(s,kx_v, ky_v, frameTimes, spiketimes, T_coords, X_coords, Y_coords, progressFcn)
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
%
%  PROGRESSFCN is an optional handle called as PROGRESSFCN(FRACTION), with
%  FRACTION from 0 to 1, at the boundaries of the two long stages: building
%  the spike-triggered average, and assessing its significance. Both are slow
%  and silent, so a caller with somewhere to display progress can pass one.
%  See ndi.fun.progressReporter. Defaults to doing nothing.

if nargin<9 || isempty(progressFcn),
	progressFcn = @(fraction) [];
end;

deltaT = mean(diff(T_coords));
rf_range = T_coords(end) - T_coords(1);
tmax = round(rf_range/deltaT) + 1;
M = size(X_coords(:), 1);
spike_num = size(spiketimes, 1);
progressFcn(0);
 % generate_STA is the long stage; let it drive the first 70% rather than
 % leaving the bar at 0 for the whole of it.
sta = vis.revcorr.generate_STA(s,kx_v, ky_v, frameTimes, spiketimes, rf_range, T_coords, tmax, M,...
	@(fraction) progressFcn(0.7*fraction));
progressFcn(0.7);
[p_val_adjusted, p_val, rescale] = vis.revcorr.calc_significance(sta,spike_num);
cmap = vis.revcorr.get_cmap();
progressFcn(1);
end

