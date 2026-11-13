function out = xi_test(xi_boot, xi_hat, xi0, options)
% XI_TEST - bootstrap test of whether the coupling index xi differs from a value
%
%   OUT = vis.speed.bootstrap.xi_test(XI_BOOT, XI_HAT, XI0, ...)
%
%   Given the bootstrap estimates XI_BOOT of the spatiotemporal coupling index
%   xi (column 3 of vis.speed.bootstrap.fit), the point estimate XI_HAT, and a
%   hypothesized value XI0 (e.g. 0 for separable tuning, 1 for speed tuning),
%   returns a two-sided p-value for H0: xi == XI0 and a confidence interval.
%
%   No bootstrap fit will equal XI0 exactly; the tests ask whether XI0 is a
%   plausible value of xi given the spread of the bootstrap estimates.
%
%   Methods ('method'):
%     'percentile' - p = min(1, 2*min(1+#{xi* <= XI0}, 1+#{xi* >= XI0}) / (S+1)),
%                    i.e. XI0 is rejected at level alpha when it lies outside
%                    the (1-alpha) percentile interval of XI_BOOT.
%     'basic'      - treats xi* - XI_HAT as the sampling error of XI_HAT,
%                    giving the interval [2*XI_HAT - q_hi, 2*XI_HAT - q_lo]
%                    (q = percentiles of XI_BOOT). The p-value is the
%                    percentile p-value evaluated at the reflected value
%                    2*XI_HAT - XI0, so XI0 is rejected at level alpha exactly
%                    when it lies outside the (1-alpha) basic interval.
%     'bc'         - bias-corrected percentile: z0 = PHI^-1(F(XI_HAT)), where
%                    F is the empirical CDF of XI_BOOT; the interval uses the
%                    percentiles PHI(2*z0 +/- z_{alpha/2}), and
%                    p = 2*(1 - PHI(|PHI^-1(F(XI0)) - 2*z0|)). Corrects for a
%                    bootstrap distribution whose median is shifted from XI_HAT.
%   (S+1 rather than S in the denominators keeps p > 0 for finite S.)
%
%   Inputs:
%     XI_BOOT - vector of bootstrap estimates of xi (NaNs, from failed fits,
%               are ignored).
%     XI_HAT  - point estimate of xi from the observed data.
%     XI0     - hypothesized value (scalar), or a vector of values (the
%               p-values and output fields are then vectors).
%
%   Name-Value Pair Arguments:
%     method          - 'percentile' (default), 'basic', or 'bc'
%     confidenceLevel - percent coverage of the interval (default 95)
%
%   Output OUT, a struct:
%     p        - two-sided p-value(s) for H0: xi == XI0
%     ci       - [low high] confidence interval for xi
%     method, confidenceLevel, numBootstrap (non-NaN estimates used), xi0
%
%   Example:
%     out0 = vis.speed.bootstrap.xi_test(B(:,3), info.pointEstimate(3), 0);
%     isInseparable = out0.p < 0.05;
%
%   See also: vis.speed.bootstrap.fit, vis.speed.bootstrap.classify

    arguments
        xi_boot double
        xi_hat (1,1) double
        xi0 double
        options.method (1,1) string {mustBeMember(options.method, ["percentile", "basic", "bc"])} = "percentile"
        options.confidenceLevel (1,1) double {mustBeGreaterThan(options.confidenceLevel, 0), mustBeLessThan(options.confidenceLevel, 100)} = 95
    end

    x = xi_boot(:);
    x = x(~isnan(x));
    S = numel(x);
    alpha = 1 - options.confidenceLevel / 100;

    out.method = char(options.method);
    out.confidenceLevel = options.confidenceLevel;
    out.numBootstrap = S;
    out.xi0 = xi0;

    if S == 0
        out.p = nan(size(xi0));
        out.ci = [NaN NaN];
        return;
    end

    p = nan(size(xi0));
    switch options.method
        case "percentile"
            for i = 1:numel(xi0)
                nLow = sum(x <= xi0(i));
                nHigh = sum(x >= xi0(i));
                p(i) = min(1, 2 * min(1 + nLow, 1 + nHigh) / (S + 1));
            end
            out.ci = localQuantile(x, [alpha/2, 1 - alpha/2]);
        case "basic"
            for i = 1:numel(xi0)
                v = 2*xi_hat - xi0(i);    % reflected value
                p(i) = min(1, 2 * min(1 + sum(x <= v), 1 + sum(x >= v)) / (S + 1));
            end
            q = localQuantile(x, [alpha/2, 1 - alpha/2]);
            out.ci = [2*xi_hat - q(2), 2*xi_hat - q(1)];
        case "bc"
            z0 = localNorminv(localEcdf(x, xi_hat, S));
            zA = localNorminv(1 - alpha/2);
            out.ci = localQuantile(x, [localNormcdf(2*z0 - zA), localNormcdf(2*z0 + zA)]);
            for i = 1:numel(xi0)
                zq = localNorminv(localEcdf(x, xi0(i), S)) - 2*z0;
                p(i) = min(1, 2 * (1 - localNormcdf(abs(zq))));
            end
    end
    out.p = p;

end

function F = localEcdf(x, v, S)
% Mid-rank empirical CDF, kept inside [1/(S+1), S/(S+1)] so its inverse
% normal is finite.
    F = (sum(x < v) + 0.5 * sum(x == v)) / S;
    F = min(max(F, 1/(S+1)), S/(S+1));
end

function q = localQuantile(x, p)
% Percentiles of x at probabilities p (linear interpolation, like prctile).
    x = sort(x(:));
    n = numel(x);
    if n == 1
        q = repmat(x, size(p));
        return;
    end
    pos = 1 + (n - 1) * min(max(p, 0), 1);
    lo = floor(pos);
    hi = min(lo + 1, n);
    w = pos - lo;
    q = x(lo).' .* (1 - w) + x(hi).' .* w;
    q = reshape(q, size(p));
end

function z = localNorminv(p)
    z = -sqrt(2) * erfcinv(2 * p);
end

function c = localNormcdf(z)
    c = 0.5 * erfc(-z / sqrt(2));
end
