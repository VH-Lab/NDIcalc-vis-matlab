function sta = generate_STA(s,kx_v, ky_v, frameTimes, spiketimes, rf_range, T_coords, tmax, M, progressFcn)
%GENERATE_STA - Generate the Spike-Triggered Average (STA) from Hartley stimulus
%
% STA = vis.revcorr.generate_STA(S, KX_V, KY_V, FRAMETIMES, SPIKETIMES, RF_RANGE, T_COORDS, TMAX, M, PROGRESSFCN)
%
% Inputs:
%  S - the s parameters of the Hartley stimulus
%  KX_V - the kx parameters of the Hartley stimulus
%  KY_V - the ky parameters of the Hartley stimulus
%  FRAMETIMES - the time points of the Hartley stimulus
%  SPIKETIMES - the times of the spikes to trigger on
%  RF_RANGE - (Unused in current implementation logic but kept for signature)
%  T_COORDS - vector of time coordinates relative to spike time for the STA (e.g. [-0.5 ... 0.5])
%  TMAX - number of time bins for the reconstruction
%  M - spatial dimension of the stimulus (MxM)
%  PROGRESSFCN - an optional handle called as PROGRESSFCN(FRACTION), with
%   FRACTION the proportion of spikes accumulated so far. This loop runs for
%   minutes on a real recording, so a caller with somewhere to show progress
%   can pass one. See ndi.fun.progressReporter. Defaults to doing nothing.
%
% Outputs:
%  STA - the calculated Spike-Triggered Average (MxMxTMAX)

%% read data
% Calculate start and end times for each spike based on the provided time coordinates.
% Note: The logic assumes T_coords defines the window around the spike.
% If T_coords = [-0.1, 0, 0.1], then:
% t_start = spike - 0.1
% t_end = spike - (-0.1) = spike + 0.1
% This assumes symmetric or specific ordering.
% Let's stick to the user provided logic:
t_start = spiketimes + T_coords(1);
t_end = spiketimes + T_coords(end);

% The previous code in file was:
% t_start = spiketimes - T_coords(end);
% t_end = spiketimes - T_coords(1);
% This implies T_coords might be "lags" (positive means past?).
% If T_coords is [-0.1 ... 0.1]
% t_start = spike - 0.1
% t_end = spike - (-0.1) = spike + 0.1
% This matches "time around spike".

% However, looking at test.m:
% t_s = t_values(i) + reconstructionT0;
% t_e = t_values(i) + reconstructionT1;
% reconstructionTimeBase = linspace(reconstructionT0, reconstructionT1, ...)
% So passing reconstructionTimeBase as T_coords means T_coords(1) is T0, T_coords(end) is T1.
% So t_s = spike + T0 = spike + T_coords(1).
% t_e = spike + T1 = spike + T_coords(end).

% The original code in generate_STA was:
% t_start = spiketimes - T_coords(end);
% t_end = spiketimes - T_coords(1);
% This is INCONSISTENT with test.m logic if T_coords is the time base.
% I will update generate_STA to match the logic in test.m which is straightforward addition of the window.

if nargin<10,
    progressFcn = [];
end;

%% reconstruction
reconstruction_block = zeros(M, M, tmax);

 % The worker limit comes from the user's NDI parallel preferences; 0 runs
 % this loop serially in the client. See ndi.fun.parallelWorkers.
numberOfWorkers = ndi.fun.parallelWorkers();

numberOfSpikes = size(t_start, 1);

 % Count completed iterations rather than reading the loop index: parfor takes
 % the range in an arbitrary order. See ndi.fun.parforProgress.
[tick, tickEvery] = ndi.fun.parforProgress(numberOfSpikes, progressFcn, numberOfWorkers);

parfor (i = 1:numberOfSpikes, numberOfWorkers)
    t_s = t_start(i);
    t_e = t_end(i);
    [hartley_stimulus_parameters, hartley_stimulus_times] = vis.revcorr.get_frames(s,kx_v, ky_v, frameTimes, t_s, t_e);
    [b,~] = vis.revcorr.hartley_stimulus_resampled_time(M, hartley_stimulus_parameters, hartley_stimulus_times, t_s, t_e, tmax);
    reconstruction_block = reconstruction_block + b;
     % Gated here, not inside tick: parfor rejects tick(i) as a sliced function
     % handle. See ndi.fun.parforProgress.
    if mod(i,tickEvery)==0 || i==numberOfSpikes
        tick();
    end
end
sta = reconstruction_block / size(spiketimes, 1);
end
