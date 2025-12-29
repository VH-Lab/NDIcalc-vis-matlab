function [response, spikeTimes] = calculateHartleyResponse(s, kx_v, ky_v, frameTimes, rf, kwargs)
%CALCULATEHARTLEYRESPONSE - Compute the response of a model cell to a Hartley stimulus
%
% [RESPONSE, SPIKETIMES] = vis.revcorr.calculateHartleyResponse(S, KX_V, KY_V, FRAMETIMES, RF, 'NAME', VALUE, ...)
%
% Calculates the linear response of a model neuron with a spatiotemporal receptive field (RF)
% to a Hartley stimulus sequence. The response is the mean dot product of the stimulus and the
% receptive field at each time point.
%
% Inputs:
%  S - the s parameters of the Hartley stimulus
%  KX_V - the kx parameters of the Hartley stimulus
%  KY_V - the ky parameters of the Hartley stimulus
%  FRAMETIMES - the time points of the Hartley stimulus frames
%  RF - the spatiotemporal receptive field (MxMxNumTimeSteps)
%
% Optional Name-Value pairs:
%  rfDeltaT - the time step of the receptive field (default: 0.005)
%  rfNumTimeSteps - the number of time steps in the receptive field (default: size(RF,3))
%  responseDeltaT - the time step for reconstructing the response (default: 0.01)
%  max_TimeBlock_StartTime - maximum start time to consider for response calculation (default: 500)
%  threshold - threshold for detecting spikes in the response (default: 1)
%  rfTimeRange - the total duration of the receptive field (default: rfDeltaT * rfNumTimeSteps)
%  Verbose - if true, print progress (default: true)
%
% Outputs:
%  RESPONSE - the continuous response trace of the model cell
%  SPIKETIMES - the times where the response exceeds the threshold

arguments
    s
    kx_v
    ky_v
    frameTimes
    rf
    kwargs.rfDeltaT (1,1) double = 0.005
    kwargs.rfNumTimeSteps (1,1) double = size(rf, 3)
    kwargs.responseDeltaT (1,1) double = 0.01
    kwargs.max_TimeBlock_StartTime (1,1) double = 500
    kwargs.threshold (1,1) double = 1
    kwargs.rfTimeRange double = []
    kwargs.Verbose (1,1) logical = true
end

rfDeltaT = kwargs.rfDeltaT;
rfNumTimeSteps = kwargs.rfNumTimeSteps;
responseDeltaT = kwargs.responseDeltaT;
max_TimeBlock_StartTime = kwargs.max_TimeBlock_StartTime;
threshold = kwargs.threshold;
Verbose = kwargs.Verbose;

if isempty(kwargs.rfTimeRange)
    rfTimeRange = rfDeltaT * size(rf,3);
else
    rfTimeRange = kwargs.rfTimeRange;
end

M = size(rf, 1);

% Define response times
responseTimes = frameTimes(1):responseDeltaT:frameTimes(end) + rfTimeRange;

% Filter response times based on max_TimeBlock_StartTime
I = responseTimes < max_TimeBlock_StartTime;
responseTimes = responseTimes(I);

numTimeSteps = size(responseTimes, 2);
response = zeros(numTimeSteps, 1);

% RF needs to be applied backwards in time for convolution/dot product
rf_backwards = rf(:, :, end:-1:1);

% Calculate response
parfor i = 1:numTimeSteps
    % Determine the time interval for the current response point
    % We need stimulus history of length rfTimeRange ending at responseTimes(i)
    t_end = responseTimes(i);
    t_start = t_end - rfTimeRange;

    [hartley_stimulus_parameters, hartley_stimulus_times] = vis.revcorr.get_frames(s, kx_v, ky_v, frameTimes, t_start, t_end);

    % Resample the stimulus to match the RF time steps
    [b, ~] = vis.revcorr.hartley_stimulus_resampled_time(M, hartley_stimulus_parameters, hartley_stimulus_times, t_start, t_end, rfNumTimeSteps);

    product = b .* rf_backwards;
    response(i) = mean(product, 'all');

    if Verbose && mod(i, 1000) == 0
        disp([num2str(100*i/numTimeSteps) '%']);
    end

    if isnan(response(i))
        error('Error. \n NaN at idx %d.', i)
    end
end

% Find spikes
peak_idx = find(response > threshold);
spikeTimes = responseTimes(peak_idx);

end
