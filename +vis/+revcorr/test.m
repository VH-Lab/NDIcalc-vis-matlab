function sta = test(kwargs)
% TEST - run the reverse correlation test
%
% STA = vis.revcorr.test('NAME', VALUE, ...)
%
% Optional Name-Value pairs:
%  filename - path to the json file (default: test file in repository)
%  M - spatial dimension (default: 200)
%  reconstruction_range - time range for reconstruction (default: 0.5)
%  reconstruction_t - time step for reconstruction (default: 0.01)
%  rf_range - time range for RF simulation (default: 0.2)

arguments
    kwargs.filename (1,:) char = fullfile(ndi.fun.ndiCalcVisPath(),'tests','+ndi','+unittest','+calc','+vis','1_hartley.json');
    kwargs.M (1,1) double = 200;
    kwargs.reconstruction_range (1,1) double = 0.5;
    kwargs.reconstruction_t (1,1) double = 0.01;
    kwargs.rf_range (1,1) double = 0.2;
    kwargs.max_TimeBlock_StartTime (1,1) double = 500;
end

% Extract parameters
filename = kwargs.filename;
M = kwargs.M;
reconstruction_range = kwargs.reconstruction_range;
reconstruction_t = kwargs.reconstruction_t;
rf_range = kwargs.rf_range;
max_TimeBlock_StartTime = kwargs.max_TimeBlock_StartTime;

deltaT = 0.005;

%% main
rf = vis.revcorr.setRF(M, 42, deltaT);
[M,~,rfTimeSteps] = size(rf);
vis.revcorr.stim_plot(rf);
[s,kx_v, ky_v, frameTimes, spiketimes] = vis.revcorr.json_file_processor(filename);
reconstruction_TimeBlock_StartTimes = frameTimes(1) - rf_range : deltaT: frameTimes(size(frameTimes,1)) + rf_range;
reconstructionTimeBlock_EndTimes = reconstruction_TimeBlock_StartTimes + rf_range;

I = reconstruction_TimeBlock_StartTimes < max_TimeBlock_StartTime;
reconstruction_TimeBlock_StartTimes = reconstruction_TimeBlock_StartTimes(I);
reconstructionTimeBlock_EndTimes = reconstructionTimeBlock_EndTimes(I);

[t0,t1] = vis.revcorr.get_t0(spiketimes,rf_range);

numTimeSteps = size(reconstruction_TimeBlock_StartTimes,2);

response = zeros(numTimeSteps,1);
parfor i = 1:numTimeSteps
    if mod(i,100) == 0
        disp([num2str(100*i/numTimeSteps) '%']);
    end
    cur_tp = [i, reconstruction_TimeBlock_StartTimes(i), reconstructionTimeBlock_EndTimes(i)];
    [hartley_stimulus_parameters, hartley_stimulus_times] = vis.revcorr.get_frames(s,kx_v, ky_v, frameTimes, reconstruction_TimeBlock_StartTimes(i), reconstructionTimeBlock_EndTimes(i));
    [b,t] = vis.revcorr.hartley_stimulus_resampled_time(M, hartley_stimulus_parameters, hartley_stimulus_times, reconstruction_TimeBlock_StartTimes(i), reconstructionTimeBlock_EndTimes(i), rfTimeSteps);

    product = b .* rf;
    response(i) = sum(product,'all');
    if isnan(response(i))
        error('Error. \n NaN at idx %d.', i)
    end
end

peak_idx = find(response > 1e6);
t_values = reconstruction_TimeBlock_StartTimes(peak_idx);
figure;
plot(reconstruction_TimeBlock_StartTimes, response, 'b-')
hold on 
plot(t_values, response(peak_idx), 'ro')
hold off

for i = 1:5 
    t_s = t_values(i);
    cur_tp = [i, t_s, t_s + rf_range];
    [hartley_stimulus_parameters, hartley_stimulus_times] = vis.revcorr.get_frames(s,kx_v, ky_v, frameTimes, t_s, t_s + rf_range);
    [b,t] = vis.revcorr.hartley_stimulus_resampled_time(M, hartley_stimulus_parameters, hartley_stimulus_times, t_s, t_s + rf_range, deltaT);
    vis.revcorr.stim_plot(b)
end

%% reconstruction
reconstruction_block = [];

for i = 1:size(t_values, 2)
% for i = 202:203
    t_s = t_values(i);
    t_e = t_s + reconstruction_range;
    cur_tp = [i, t_s, t_e];
    disp(cur_tp);
    [hartley_stimulus_parameters, hartley_stimulus_times] = vis.revcorr.get_frames(s,kx_v, ky_v, frameTimes, t_s, t_e);
    [b,t] = vis.revcorr.hartley_stimulus_resampled_time(M, hartley_stimulus_parameters, hartley_stimulus_times, t_s, t_e, reconstruction_t);
    reconstruction_block = cat(4, reconstruction_block, b);  
end
sta = vis.revcorr.recover_sta(reconstruction_block);
vis.revcorr.stim_plot(sta);

end
