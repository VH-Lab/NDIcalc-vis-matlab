function [rescale] = rescale_p_image(cdf_image, varargin)
%RESCALE_P_IMAGE - rescale the P-value image so that only very significant points are plotted
%
% RESCALE = vis.revcorr.rescale_p_image(CDF_IMAGE, ...)
%
% Rescales an image of CDF values so that highly significant pixels are colorized.
% 
% This command takes names/value pairs that modify its behavior:
% -----------------------------------------------------------------------|
% | Parameter (default)         | Description                            |
% -----------------------------------------------------------------------|
% | low_cut_off (log10(2*1e-3)  | Cut off below which P values are shown |
% |                             |  or 1-P values are shown               |
% | high_cut_off (log10(2*1e-10)| P-value where maximum brightness is    |
% |                             |  used.                                 |
% | colortableentries (256)     | Number of color table entries          |
% | low_entries ([1 64])        | The low entries for low p_values       |
% | high_entries ([193 256])    | The high entries for low 1-p_values    |
% -----------------------------------------------------------------------|
%

low_cut_off = log10(2*1e-5);
high_cut_off = log10(2*1e-10);
colortableentries = 256;
low_entries = [0 64];
high_entries = [193 256];

vlt.data.assign(varargin{:});

rescale = 128+zeros(size(cdf_image));

indexes_high = find(cdf_image >= 0.5);
indexes_low = find(cdf_image < 0.5);

cdf_adjust = cdf_image;

 % log transform
cdf_adjust(indexes_high) = -1 * log10(2 * (1-cdf_image(indexes_high)));
cdf_adjust(indexes_low)  = log10(2*(cdf_image(indexes_low)));

indexes_high = find(cdf_adjust >= -low_cut_off);
rescale(indexes_high) = vlt.math.rescale(cdf_adjust(indexes_high), [-low_cut_off -high_cut_off], high_entries);

indexes_low = find(cdf_adjust <= low_cut_off);
rescale(indexes_low) = vlt.math.rescale(cdf_adjust(indexes_low), [high_cut_off low_cut_off], low_entries);


