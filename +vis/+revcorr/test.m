%% set the parameters
filename = "/Users/cxy/Reverse Correlation/data/hartley_data/t00005_leftcortex_1 | 1_hartley.json";
rf_range = 0.2;
deltaT = 0.005;
M = 200;
reconstruction_range = 0.5;
reconstruction_t = 0.01;
%% main
rf = vis.revcorr.setRF();
vis.revcorr.stim_plot(rf);
[s,kx_v, ky_v, frameTimes, spiketimes] = vis.revcorr.json_file_processor(filename);
t_start = frameTimes(1) - rf_range : deltaT: frameTimes(size(frameTimes)) + rf_range;
t_end = t_start + rf_range;
[t0,t1] = vis.revcorr.get_t0(spiketimes,rf_range);

response = [];
for i = 1:size(t_start, 2)
    cur_tp = [i, t_start(i), t_end(i)];
    disp(cur_tp);
    [hartley_stimulus_parameters, hartley_stimulus_times] = vis.revcorr.get_frames(s,kx_v, ky_v, frameTimes, t_start(i), t_end(i));
    [b,t] = vis.revcorr.hartley_stimulus_resampled_time(M, hartley_stimulus_parameters, hartley_stimulus_times, t_start(i), t_end(i), deltaT);

    product = b .* rf;
    response(end+1) = sum(product,'all');
    if isnan(response(end))
        error('Error. \n NaN at idx %d.', i)
        
    end
end

peak_idx = find(response > 1e6);
t_values = t_start(peak_idx);
figure;
plot(t_start, response, 'b-')
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

