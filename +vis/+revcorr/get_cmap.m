function cmap = get_cmap()
%GET_CMAP Summary - generate the color map for plotting the p_value
%
% CMAP = GET_CMAP();
%
% Returns a 256-entry color map where the first 64 entries are blue
% going from strong (values near 0) to weak (entry 64), the middle values
% are 1, and the high values from 193 to 256 are increasing intensity of
% red, such that 256 has the highest red intensity.
%

N = 256;
cmap = zeros(N, 3);
cmap(1:64, 3) = linspace(1, 0, 64);
cmap(1:64, 2) = 0; 
cmap(1:64, 1) = 0; 
cmap(193:256, 2) = 0; 
cmap(193:256, 1) = linspace(0, 1, 64);
cmap(193:256, 3) = 0;

end

