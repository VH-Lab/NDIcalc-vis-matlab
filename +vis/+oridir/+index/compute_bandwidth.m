function [bandwidth, low, high] = compute_bandwidth(x, y, options)
% COMPUTE_BANDWIDTH - Compute tuning bandwidth in octaves
%
%   [BANDWIDTH, LOW, HIGH] = vis.oridir.index.compute_bandwidth(X, Y)
%   [BANDWIDTH, LOW, HIGH] = vis.oridir.index.compute_bandwidth(X, Y, 'Period', P)
%
%   Calculates the bandwidth (in octaves) of a response curve defined by
%   stimulus values X and mean response values Y.
%
%   Inputs:
%     X - Array of independent variable values (e.g. stimulus values,
%         spatial/temporal frequencies, or angles).
%     Y - Array of response mean values corresponding to X.
%
%   Optional name/value inputs:
%     Period - The period of X, for circular X such as directions
%              (default Inf, not circular). See vis.oridir.index.compute_fwhh.
%
%   Outputs:
%     BANDWIDTH - Bandwidth computed in octaves as log2(HIGH / LOW), where
%                 LOW and HIGH are the values of X on either side of the peak
%                 where the response first falls to half of its maximum
%                 height, walking outward from the peak. If either LOW or
%                 HIGH cannot be determined (e.g., the response does not
%                 fall below half maximum), BANDWIDTH is returned as Inf.
%                 With a finite Period, LOW and HIGH are the values reported
%                 modulo the Period, so when the half-height points lie on
%                 either side of the wrap point LOW is greater than HIGH and
%                 BANDWIDTH is negative.
%     LOW, HIGH - The half-height points (NaN when not found). See
%                 vis.oridir.index.compute_fwhh.
%
%   See also: vis.oridir.index.compute_fwhh

arguments
    x {mustBeNumeric}
    y {mustBeNumeric}
    options.Period (1,1) double {mustBePositive} = Inf
end

[~, low, high] = vis.oridir.index.compute_fwhh(x, y, 'Period', options.Period);

if isnan(low) || isnan(high)
    bandwidth = Inf;
else
    bandwidth = log2(high/low);
end
