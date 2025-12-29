function [response, spikeTimes] = calculateHartleyResponse(s, kx_v, ky_v, frameTimes, rf, varargin)
% CALCULATEHARTLEYRESPONSE - compute the response of a model cell to a Hartley stimulus
%
% [RESPONSE, SPIKETIMES] = vis.revcorr.calculateHartleyResponse(S, KX_V, KY_V, FRAMETIMES, RF, ...)
%
% Calculates the response of a model cell with receptive field RF to a Hartley
% stimulus sequence defined by S, KX_V, KY_V, and FRAMETIMES.
%
% Inputs:
%  S - Vector of S parameters for the Hartley stimuli
%  KX_V - Vector of kx parameters
%  KY_V - Vector of ky parameters
%  FRAMETIMES - Vector of frame onset times
%  RF - Receptive field (MxMxTime)
%
% Name-Value Pairs:
%  'Verbose' - (default 0)
%
% Outputs:
%  RESPONSE - The continuous response trace
%  SPIKETIMES - Times of generated spikes (threshold crossing)

Verbose = 0;
vlt.data.assign(varargin{:});

M = size(rf, 1);
num_timesteps_rf = size(rf, 3);

% Parameters from test.m defaults (or inferred)
rf_range = 0.2; % From test.m
deltaT = 0.005; % From test.m
t_start = frameTimes(1) - rf_range;
t_end = frameTimes(end) + rf_range;

% Let's define the simulation time vector
t_sim = t_start : deltaT : t_end;
response = zeros(size(t_sim));

hartley_stimulus_parameters.s = s;
hartley_stimulus_parameters.kx_v = kx_v;
hartley_stimulus_parameters.ky_v = ky_v;

for i = 1:length(t_sim)
    curr_t = t_sim(i);
    % We need the stimulus history relevant to the RF at this time point.
    % RF has duration related to num_timesteps_rf and deltaT?
    % In test.m: rf is 200x200x41. And deltaT=0.005. So 41*0.005 = 0.205s.
    % t_start and t_end passed to hartley_stimulus_resampled_time determine the window.
    % The window should match RF duration.

    t0_window = curr_t;
    t1_window = curr_t + rf_range; % This seems to be "looking forward" or the test.m logic is specific?
    % In test.m: t_start = frameTimes(1) - rf_range : deltaT: ...
    % t_end = t_start + rf_range.
    % So it takes a window of size rf_range starting at t_start.
    % And computes dot product with RF.

    num_steps = size(rf, 3);

    % Calculate response loop
    % Get frames in window [curr_t, curr_t + rf_range]
    [frames_p, frames_t] = vis.revcorr.get_frames(s, kx_v, ky_v, frameTimes, curr_t, curr_t + rf_range);

    % Resample. We need exactly `num_steps` frames to match RF.
    [b, ~] = vis.revcorr.hartley_stimulus_resampled_time(M, frames_p, frames_t, curr_t, curr_t + rf_range, num_steps);

    % Dot product
    product = b .* rf;
    response(i) = sum(product, 'all');
end

% Threshold for spikes
% In test.m: `peak_idx = find(response > 1e6);`
% I'll use 1e6 as threshold.
threshold = 1e6;
spikeTimes = t_sim(response > threshold);

end
