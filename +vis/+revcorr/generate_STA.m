function sta = generate_STA(s,kx_v, ky_v, frameTimes, spiketimes, rf_range, T_coords, tmax, M)
%GENERATE_STA Summary of this function goes here
%   Detailed explanation goes here

%% read data
t_start = spiketimes - T_coords(end); 
t_end = spiketimes - T_coords(1);
%% reconstruction
reconstruction_block = zeros(M, M, tmax);

parfor i = 1:size(t_start, 1)
    t_s = t_start(i);
    t_e = t_end(i);
    cur_tp = [i, t_s, t_e];
    disp(cur_tp);
    [hartley_stimulus_parameters, hartley_stimulus_times] = revcorr.get_frames(s,kx_v, ky_v, frameTimes, t_s, t_e);
    [b,t] = revcorr.hartley_stimulus_resampled_time(M, hartley_stimulus_parameters, hartley_stimulus_times, t_s, t_e, tmax);
%     if (i == 1)
%         reconstruction_block = zeros(size(b));
%     end
    reconstruction_block = reconstruction_block + b;  
end
sta = reconstruction_block / size(spiketimes, 1);
% revcorr.stim_plot(sta);
% file_name = strsplit(input_filename, "/");
% file_name = strsplit(file_name(end), ".");
% file_name = file_name(1);
% output_file_name = output_dir + "/" + file_name(end) + ".mat";
% save(output_file_name, 'sta');
end

