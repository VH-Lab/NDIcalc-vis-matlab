function sta = test(kwargs)
% TEST - run the reverse correlation test
%
% STA = vis.revcorr.test('NAME', VALUE, ...)
%
% Optional Name-Value pairs:
%  filename - path to the json file (default: test file in repository)
%  M - spatial dimension (default: 200)
%  reconstructionT0 - start time for reconstruction (default: -0.5)
%  reconstructionT1 - end time for reconstruction (default: 0.5)
%  reconstruction_deltat - time step for reconstruction (default: 0.01)
%  rf_range - time range for RF simulation (default: 0.2)
%  rfDeltaT - time step for RF (default: 0.005)
%  rfNumTimeSteps - number of time steps for RF (default: 42)
%  stimPlotT0 - start time for stimulus plotting (default: -0.500)
%  stimPlotT1 - end time for stimulus plotting (default: 0.5)
%  stimPlotDeltaT - time step for stimulus plotting (default: 0.005)
%  responseDeltaT - the time step for response reconstruction
%  threshold - threshold for detecting spikes (default: 1)
%  Verbose - verbose output (default: true)

arguments
    kwargs.filename (1,:) char = fullfile(ndi.fun.ndiCalcVisPath(),'tests','+ndi','+unittest','+calc','+vis','1_hartley.json');
    kwargs.M (1,1) double = 200;
    kwargs.reconstructionT0 (1,1) double = -0.5;
    kwargs.reconstructionT1 (1,1) double = 0.5;
    kwargs.reconstruction_deltat (1,1) double = 0.01;
    kwargs.rfTimeRange (1,1) double = 0.2;
    kwargs.max_TimeBlock_StartTime (1,1) double = 500;
    kwargs.rfDeltaT (1,1) double = 0.005;
    kwargs.rfNumTimeSteps (1,1) double = 40;
    kwargs.stimPlotT0 (1,1) double = -0.100;
    kwargs.stimPlotT1 (1,1) double = 0.100;
    kwargs.stimPlotDeltaT (1,1) double = 0.01;
    kwargs.responseDeltaT (1,1) double = 0.01;
    kwargs.threshold (1,1) double = 1;
    kwargs.Verbose (1,1) logical = true;
end

% Extract parameters
filename = kwargs.filename;
M = kwargs.M;
reconstructionT0 = kwargs.reconstructionT0;
reconstructionT1 = kwargs.reconstructionT1;
reconstruction_deltat = kwargs.reconstruction_deltat;
rfTimeRange = kwargs.rfTimeRange;
max_TimeBlock_StartTime = kwargs.max_TimeBlock_StartTime;
rfDeltaT = kwargs.rfDeltaT;
rfNumTimeSteps = kwargs.rfNumTimeSteps;
stimPlotT0 = kwargs.stimPlotT0;
stimPlotT1 = kwargs.stimPlotT1;
stimPlotDeltaT = kwargs.stimPlotDeltaT;
responseDeltaT = kwargs.responseDeltaT;
threshold = kwargs.threshold;
Verbose = kwargs.Verbose;

stimReconstructionSteps = ceil((stimPlotT1 - stimPlotT0) / stimPlotDeltaT);

%% main
[rf,rfTimeLags] = vis.revcorr.setRF(M, rfNumTimeSteps, rfDeltaT);
[M,~,rfTimeSteps] = size(rf);
vis.revcorr.stim_plot(rf,[],rfTimeLags);
[s,kx_v, ky_v, frameTimes, ~] = vis.revcorr.json_file_processor(filename);

[response, spiketimes] = vis.revcorr.calculateHartleyResponse(s, kx_v, ky_v, frameTimes, rf, ...
    'rfDeltaT', rfDeltaT, 'rfNumTimeSteps', rfNumTimeSteps, 'responseDeltaT', responseDeltaT, ...
    'max_TimeBlock_StartTime', max_TimeBlock_StartTime, 'threshold', threshold, 'rfTimeRange', rfTimeRange, 'Verbose', Verbose);

% Plot response
responseTimes = frameTimes(1):responseDeltaT:frameTimes(end) + rfTimeRange;
I = responseTimes < max_TimeBlock_StartTime;
responseTimes = responseTimes(I);

figure;
plot(responseTimes(1:length(response)), response, 'b-')
hold on 
% Find indices of spike times in the response vector for plotting 'ro' is harder without exact indices.
% However, spiketimes are the times. We can plot them directly against the threshold/response value?
% Since spiketimes are the times, we need the corresponding response values.
% Let's find indices in our local `responseTimes` that match `spiketimes` (approximately).
peak_vals = interp1(responseTimes(1:length(response)), response, spiketimes, 'nearest');
plot(spiketimes, peak_vals, 'ro')
hold off

for i = 1:5 
    [hartley_stimulus_parameters, hartley_stimulus_times] = vis.revcorr.get_frames(s,kx_v, ky_v, frameTimes, spiketimes(i) + stimPlotT0, spiketimes(i) + stimPlotT1);
    [b,t] = vis.revcorr.hartley_stimulus_resampled_time(M, hartley_stimulus_parameters, hartley_stimulus_times, spiketimes(i) + stimPlotT0, spiketimes(i) + stimPlotT1, stimReconstructionSteps);
    vis.revcorr.stim_plot(b,[],t-spiketimes(i));
end

%% reconstruction
reconstructionTimeBins = (reconstructionT1 - reconstructionT0) / reconstruction_deltat;
reconstructionTimeBase = linspace(reconstructionT0, reconstructionT1, reconstructionTimeBins);

sta = vis.revcorr.generate_STA(s, kx_v, ky_v, frameTimes, spiketimes, rfTimeRange, reconstructionTimeBase, reconstructionTimeBins, M);
vis.revcorr.stim_plot(sta, [], reconstructionTimeBase);

end
