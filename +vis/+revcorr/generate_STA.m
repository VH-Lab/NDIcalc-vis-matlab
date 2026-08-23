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
%  T_COORDS - vector of lags, measured backwards from each event: a lag of
%   +0.100 means the stimulus 0.100 s before the event, -0.100 the stimulus
%   0.100 s after it. The reconstruction is returned in reverse-lag order,
%   so plane s holds the stimulus at event-T_COORDS(TMAX+1-s); callers such
%   as ndi.calc.vis.hartley reverse the third dimension so that plane k
%   corresponds to T_COORDS(k).
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
 % T_coords holds lags measured backwards from each event: a lag of +0.100 s
 % means the stimulus that was on the screen 0.100 s BEFORE the event, and a
 % lag of -0.100 s the stimulus 0.100 s after it. The window is therefore
 % taken backwards, from event-T_coords(end) to event-T_coords(1), which
 % fills the reconstruction in reverse-lag order: plane s holds the stimulus
 % at event-T_coords(TMAX+1-s). ndi.calc.vis.hartley reverses the third
 % dimension of the result, after which plane k matches T_coords(k) element
 % for element.
 %
 % Do not "simplify" this to event+T_coords(1) .. event+T_coords(end). That
 % form agrees with this one only when T_coords is symmetric about zero, and
 % neither of the calculator's defaults is; it was in place between c47d6ca
 % and this commit, and it offset the stored time axis by T_coords(1)+
 % T_coords(end). See tests/+vis/+revcorr/test_generate_STA.m, which pins
 % the convention.
t_start = spiketimes - T_coords(end);
t_end = spiketimes - T_coords(1);

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
