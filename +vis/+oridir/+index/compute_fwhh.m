function [fwhh, low, high] = compute_fwhh(x, y, options)
% COMPUTE_FWHH - Compute full width at half height of a response curve
%
%   [FWHH, LOW, HIGH] = vis.oridir.index.compute_fwhh(X, Y)
%   [FWHH, LOW, HIGH] = vis.oridir.index.compute_fwhh(X, Y, 'Period', P)
%
%   Calculates the full width at half height of a response curve defined by
%   stimulus values X and mean response values Y.
%
%   Starting at the maximum response, the curve is walked outward toward
%   lower X and toward higher X, and on each side the first point where it
%   falls to half of the maximum is taken. The crossing is located by linear
%   interpolation between the two samples that bracket it. Because each walk
%   stops at its first crossing, the result is well defined even when the
%   response rises and falls again farther from the peak (for example, the
%   response to the null direction).
%
%   If 'Period' is given (for example 360 for directions in degrees), X is
%   treated as circular: the data are laid out three times, at X-P, X and
%   X+P, and the walks start from the peak in the middle copy, so that each
%   walk can continue across the wrap point for up to a full period. LOW and
%   HIGH are then reported modulo P, while FWHH is measured before that
%   wrapping, so a curve peaked at 0 with half-height points at 340 and 20
%   has FWHH = 40, LOW = 340 and HIGH = 20. By default the Period is Inf and
%   X is not treated as circular.
%
%   Inputs:
%     X - Array of independent variable values (e.g. directions, in degrees).
%     Y - Array of response mean values corresponding to X.
%
%   Optional name/value inputs:
%     Period - The period of X (default Inf, not circular).
%
%   Outputs:
%     FWHH - The full width at half height, in the units of X. If either
%            LOW or HIGH cannot be determined, FWHH is Inf.
%     LOW  - The value of X below the peak where the response first falls
%            to half of its maximum, or NaN if it does not fall that far.
%     HIGH - The value of X above the peak where the response first falls
%            to half of its maximum, or NaN if it does not fall that far.
%
%   If the maximum response is not positive, LOW and HIGH are NaN.
%
%   See also: vis.oridir.index.compute_bandwidth

arguments
    x {mustBeNumeric}
    y {mustBeNumeric}
    options.Period (1,1) double {mustBePositive} = Inf
end

P = options.Period;

x = x(:);
y = y(:);

good = ~isnan(x) & ~isnan(y);
x = x(good);
y = y(good);

if isfinite(P)
    % put every sample in [0, P) and drop repeats such as 0 and P
    x = mod(x, P);
    [x, ia] = unique(x, 'stable');
    y = y(ia);
end

[x, order] = sort(x);
y = y(order);

low = NaN;
high = NaN;

if ~isempty(y)
    [ymax, imax] = max(y);
    if ymax > 0
        half = ymax / 2;

        if isfinite(P)
            n = numel(x);
            x = [x - P; x; x + P];
            y = [y; y; y];
            imax = imax + n; % start from the peak in the middle copy
        end

        % walk from the peak toward lower x
        k = find(y(1:imax) <= half, 1, 'last');
        if ~isempty(k)
            low = x(k) + (half - y(k)) * (x(k+1) - x(k)) / (y(k+1) - y(k));
        end

        % walk from the peak toward higher x
        k = find(y(imax:end) <= half, 1, 'first');
        if ~isempty(k)
            k = k + imax - 1;
            high = x(k-1) + (half - y(k-1)) * (x(k) - x(k-1)) / (y(k) - y(k-1));
        end
    end
end

if isnan(low) || isnan(high)
    fwhh = Inf;
else
    fwhh = high - low;
end

if isfinite(P)
    low = mod(low, P);
    high = mod(high, P);
end
