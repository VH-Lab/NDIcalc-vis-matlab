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

stimReconstructionSteps = ceil((stimPlotT1 - stimPlotT0) / stimPlotDeltaT);

%% main
[rf,rfTimeLags] = vis.revcorr.setRF(M, rfNumTimeSteps, rfDeltaT);
[M,~,rfTimeSteps] = size(rf);
vis.revcorr.stim_plot(rf,[],rfTimeLags);
[s,kx_v, ky_v, frameTimes, ~] = vis.revcorr.json_file_processor(filename);
responseTimes = frameTimes(1): responseDeltaT: frameTimes(size(frameTimes,1)) + rfTimeRange;

I = responseTimes < max_TimeBlock_StartTime;
responseTimes = responseTimes(I);

numTimeSteps = size(responseTimes,2);

response = zeros(numTimeSteps,1);
rf_backwards = rf(:,:,end:-1:1);
for i = 1:numTimeSteps
    if mod(i,100) == 0
        disp([num2str(100*i/numTimeSteps) '%']);
    end
    [hartley_stimulus_parameters, hartley_stimulus_times] = vis.revcorr.get_frames(s,kx_v, ky_v, frameTimes, responseTimes(i), responseTimes(i)-rfTimeRange);
    [b,t] = vis.revcorr.hartley_stimulus_resampled_time(M, hartley_stimulus_parameters, hartley_stimulus_times, responseTimes(i), responseTimes(i)-rfTimeRange, rfTimeSteps);

    product = b .* rf_backwards; % rf goes backward in time
    response(i) = mean(product,'all');
    if isnan(response(i))
        error('Error. \n NaN at idx %d.', i)
    end
end

threshold = 1;
peak_idx = find(response > threshold);
t_values = responseTimes(peak_idx);
figure;
plot(responseTimes, response, 'b-')
hold on 
plot(t_values, response(peak_idx), 'ro')
hold off

for i = 1:5 
    [hartley_stimulus_parameters, hartley_stimulus_times] = vis.revcorr.get_frames(s,kx_v, ky_v, frameTimes, t_values(i) + stimPlotT0, t_values(i) + stimPlotT1);
    [b,t] = vis.revcorr.hartley_stimulus_resampled_time(M, hartley_stimulus_parameters, hartley_stimulus_times, t_values(i) + stimPlotT0, t_values(i) + stimPlotT1, stimReconstructionSteps);
    vis.revcorr.stim_plot(b,[],t-t_values(i));
end

%% reconstruction
reconstruction_block = [];

reconstructionTimeBins = (reconstructionT1 - reconstructionT0) / reconstruction_deltat;
reconstructionTimeBase = linspace(reconstructionT0, reconstructionT1, reconstructionTimeBins);

for i = 1:size(t_values, 2)
    t_s = t_values(i) + reconstructionT0;
    t_e = t_values(i) + reconstructionT1;
    cur_tp = [i, t_s, t_e];
    disp(cur_tp);
    [hartley_stimulus_parameters, hartley_stimulus_times] = vis.revcorr.get_frames(s,kx_v, ky_v, frameTimes, t_s, t_e);
    [b,t] = vis.revcorr.hartley_stimulus_resampled_time(M, hartley_stimulus_parameters, hartley_stimulus_times, t_s, t_e, reconstructionTimeBins);
    reconstruction_block = cat(4, reconstruction_block, b);  
end
sta = vis.revcorr.recover_sta(reconstruction_block);
vis.revcorr.stim_plot(sta, [], reconstructionTimeBase);

end
