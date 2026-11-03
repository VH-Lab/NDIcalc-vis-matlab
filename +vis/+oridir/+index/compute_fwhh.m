function [fwhh, low, high] = compute_fwhh(x, y)
% COMPUTE_FWHH - Compute full width at half height of a response curve
%
%   [FWHH, LOW, HIGH] = vis.oridir.index.compute_fwhh(X, Y)
%
%   Calculates the full width at half height of a response curve defined by
%   stimulus values X and mean response values Y.
%
%   Starting from the maximum response, the curve is followed outward in each
%   direction to the first point where it falls to half of the maximum. The
%   crossing is located by linear interpolation between the two samples that
%   bracket it. Because the search stops at the first crossing, a second peak
%   elsewhere in the curve (such as the response to the null direction) does
%   not affect the result.
%
%   Inputs:
%     X - Array of independent variable values (e.g. directions, in degrees).
%     Y - Array of response mean values corresponding to X.
%
%   Outputs:
%     FWHH - The full width at half height, HIGH - LOW, in the units of X.
%            If either LOW or HIGH cannot be determined, FWHH is Inf.
%     LOW  - The value of X below the peak where the response first falls to
%            half of its maximum, or NaN if it does not fall that far within
%            the range of X.
%     HIGH - The value of X above the peak where the response first falls to
%            half of its maximum, or NaN if it does not fall that far within
%            the range of X.
%
%   The search covers the range of X only; it does not wrap around the circle.
%   If the maximum response is not positive, LOW and HIGH are NaN.
%
%   See also: vis.oridir.index.compute_bandwidth

x = x(:);
y = y(:);

good = ~isnan(x) & ~isnan(y);
x = x(good);
y = y(good);

[x, order] = sort(x);
y = y(order);

low = NaN;
high = NaN;

if ~isempty(y)
    [ymax, imax] = max(y);
    if ymax > 0
        half = ymax / 2;

        k = find(y(1:imax) <= half, 1, 'last');
        if ~isempty(k)
            low = x(k) + (half - y(k)) * (x(k+1) - x(k)) / (y(k+1) - y(k));
        end

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
