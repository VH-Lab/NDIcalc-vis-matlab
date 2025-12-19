% input_filename = "/Users/cxy/Reverse Correlation/data/hartley_data/t00005_leftcortex_5 | 1_hartley.json";
output_dir = "/Users/cxy/Reverse Correlation/matlab data/STA";
input_filename = "/Users/cxy/Reverse Correlation/data/hartley_data/t00005_leftcortex_11 | 1_hartley.json";

% rf_range = 0.5;
% deltaT = 0.01;
% M = 200;

%%
[s,kx_v, ky_v, frameTimes, spiketimes, T_coords, Y_coords] = vis.revcorr.json_file_processor(input_filename);
%%
[sta,p_val, rescale, cmap] = vis.revcorr.sta_pipeline(s,kx_v, ky_v, frameTimes, spiketimes, T_coords, Y_coords, Y_coords);
%%
% vis.revcorr.stim_plot(rescale, cmap);
%%
[X,Y,T] = meshgrid(Y_coords,Y_coords, T_coords);
[a, mu, C] = vis.revcorr.gauss3fit(X, Y, T, sta, 1);
r = a * mvnpdf([x_col y_col t_col], mu, C);
r_on = reshape(r, 200, 200, size(T_coords, 1));
[a, mu, C] = vis.revcorr.gauss3fit(X, Y, T, sta, 0);
r = a * mvnpdf([x_col y_col t_col], mu, C);
r_off = reshape(r, 200, 200, size(T_coords, 1));
sta_fit = cat(3, sta, r_on, r_off);
vis.revcorr.stim_plot(sta_fit);
% imagesc(X(1, :, 1), Y(:, 1, 1)', r_(:, :, 2));