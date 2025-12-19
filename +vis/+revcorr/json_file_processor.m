function [s,kx_v, ky_v, frameTimes, spiketimes, T_coords, Y_coords] = json_file_processor(filename)
% JSON_FILE_PROCESSOR - generate a list of useful data for reverse
% correlation from the json file.
%
% [S,KX_V, KY_V, FRAMETIMES, SPIKETIMES] = vis.revcorr.json_file_processor(FILENAME)
%
% Inputs:
%  FILENAME - the directory of JSON file to be processed
% 
% OUTPUTS:
%  S - the sign for the reconstruction
%  KX_V - the kx vector for the reconstruction
%  KY_V - the ky vector for the reconstruction
%  FRAMETIMES - the time point at which each frame appeared
%  SPIKETIMES - the time that this cell spikes during the experiment

g = jsondecode(vlt.file.textfile2char(filename));
spiketimes = g.spiketimes;
frameTimes = g.frameTimes;
s = g.hartley_numbers.S;
kx_v = g.hartley_numbers.KXV;
ky_v = g.hartley_numbers.KYV;
T_coords = g.reconstruction_properties.T_coords;
Y_coords = g.reconstruction_properties.Y_coords;

end