function [low,maxv,high]= compute_halfwidth_interp( x, y )

% COMPUTE_HALFWIDTH_INTERP - find the half-maximum points either side of the peak
%     [LOW,MAX,HIGH] = vis.frequency.compute_halfwidth_interp( X, Y )
%
%     Returns MAX, the X position where Y attains its maximum value, and
%     LOW < MAX and HIGH > MAX, where the curve falls to half of its maximum.
%
%     Starting at the peak, the curve is walked outward toward lower X and
%     toward higher X, and on each side the first point where it falls to
%     half of the maximum is taken; the crossing is located by linear
%     interpolation between the two samples that bracket it. Because each
%     walk stops at its first crossing, a response that rises and falls again
%     farther from the peak does not affect the result.
%
%     Returns NaN for LOW and/or HIGH when the curve does not fall to half of
%     its maximum on that side (or when the maximum is not positive).
%     X is not treated as circular.
%
%     Earlier versions interpolated onto a 500-point grid and took the grid
%     point nearest half height anywhere on each side of the peak, which on a
%     curve that goes up and down could pick a point far from the peak.
%
%     See also: vis.oridir.index.compute_fwhh

x = x(:);
y = y(:);

[~, imax] = max(y);
maxv = x(imax);

[~, low, high] = vis.oridir.index.compute_fwhh(x, y);
