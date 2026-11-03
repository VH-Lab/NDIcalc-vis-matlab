function [bandwidth, low, high] = compute_bandwidth(x, y)
% COMPUTE_BANDWIDTH - Compute tuning bandwidth in octaves
%
%   [BANDWIDTH, LOW, HIGH] = vis.oridir.index.compute_bandwidth(X, Y)
%
%   Calculates the bandwidth (in octaves) of a response curve defined by
%   stimulus values X and mean response values Y.
%
%   Inputs:
%     X - Array of independent variable values (e.g. stimulus values,
%         spatial/temporal frequencies, or angles).
%     Y - Array of response mean values corresponding to X.
%
%   Outputs:
%     BANDWIDTH - Bandwidth computed in octaves as log2(HIGH / LOW), where
%                 LOW and HIGH are the values of X on either side of the peak
%                 where the response first falls to half of its maximum
%                 height. If either LOW or HIGH cannot be determined (e.g.,
%                 the response does not fall below half maximum), BANDWIDTH
%                 is returned as Inf.
%     LOW, HIGH - The half-height points (NaN when not found). See
%                 vis.oridir.index.compute_fwhh.
%
%   See also: vis.oridir.index.compute_fwhh

[~, low, high] = vis.oridir.index.compute_fwhh(x, y);

if isnan(low) || isnan(high)
    bandwidth = Inf;
else
    bandwidth = log2(high/low);
end
