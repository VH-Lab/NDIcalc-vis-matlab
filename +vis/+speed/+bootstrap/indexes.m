function out = indexes(B, pointEstimate, options)
% INDEXES - spatiotemporal coupling indexes and classifications from bootstrap fits
%
%   OUT = vis.speed.bootstrap.indexes(B, POINTESTIMATE, ...)
%
%   Computes every coupling index used to describe a site from the bootstrap
%   fits of vis.speed.bootstrap.fit, without refitting, so the indexes (and
%   their thresholds) can be recomputed later from stored fits.
%
%   Inputs:
%     B             - S x 9 bootstrap matrix from vis.speed.bootstrap.fit
%                     (columns: A zeta xi sigma_sf sigma_tf sf0 tf0 sse r_squared);
%                     NaN rows (failed fits) are ignored.
%     POINTESTIMATE - 7-element parameters of the fit to the observed data
%                     [A zeta xi sigma_sf sigma_tf sf0 tf0].
%
%   Name-Value Pair Arguments (thresholds; all recorded in OUT.parameters).
%   The four threshold parameters (deltaInseparable, deltaSpeedTuned,
%   tfShiftOctaves, speedConstancyOctaves) may be VECTORS, so several values
%   can be evaluated at once; the corresponding outputs are then logical
%   (or numeric) row vectors with one element per threshold value.
%     alpha                 - significance level. (0.05) Two-sided tests use
%                             the (1-alpha) interval; equivalence and
%                             minimum-effect tests use the (1-2*alpha)
%                             interval (two one-sided tests at level alpha).
%     method                - interval/p-value method for xi:
%                             'percentile' (default), 'basic' or 'bc'
%                             (vis.speed.bootstrap.xi_test). Derived
%                             quantities always use percentile intervals.
%     deltaInseparable      - minimum coupling magnitude for the
%                             minimum-effect inseparability rule. (0.2)
%     deltaSpeedTuned       - equivalence margin around xi = 1 for the
%                             equivalence speed-tuning rule. (0.5)
%     tfShiftOctaves        - passband inseparability: the preferred temporal
%                             frequency must shift by MORE than this many
%                             octaves across the cell's SF passband. (1)
%     speedConstancyOctaves - passband speed tuning: the preferred speed must
%                             change by LESS than this many octaves across the
%                             cell's SF passband. (1)
%     min_xi, max_xi        - fit range used for xi, to flag samples at a
%                             bound. (-4, 4)
%     sfRange               - [lowest highest] stimulus spatial frequency, to
%                             flag fits whose sf0 sits at an edge of the
%                             sampled range. ([] = not computed)
%
%   Background for the passband rules. In the Priebe et al. (2006) model the
%   preferred temporal frequency varies as SF^xi, so across a span of B
%   octaves of spatial frequency
%       preferred TF shifts by       xi * B      octaves, and
%       preferred speed changes by  (xi - 1) * B octaves.
%   B is taken as the cell's SF passband, the full width at half maximum of its
%   SF tuning: B = 2*sqrt(2*ln 2) * sigma_sf / log10(2) octaves (sigma_sf is in
%   log10 units). Each bootstrap sample uses its own xi and sigma_sf, so the
%   uncertainty in B is included.
%
%   Output OUT, a struct:
%     .xi            estimate, bootstrapMedian, se, ci, ci90, ciWidth,
%                    p0, p1 (p-values for xi == 0 and xi == 1),
%                    probPositive, probNegative (fractions of bootstrap xi
%                    > 0 and < 0), fracAtBound
%     .inseparable   differentFrom0  - p0 < alpha
%                    minimumEffect   - (1-2a) interval of xi entirely outside
%                                      [-deltaInseparable, deltaInseparable]
%                    passband        - (1-2a) interval of xi*B entirely outside
%                                      [-tfShiftOctaves, tfShiftOctaves]
%     .speedTuned    priebe          - p0 < alpha and p1 >= alpha (95% interval
%                                      excludes 0, includes 1; Priebe et al. 2006)
%                    equivalence     - differentFrom0 and (1-2a) interval of xi
%                                      inside [1-deltaSpeedTuned, 1+deltaSpeedTuned]
%                    passband        - differentFrom0 and (1-2a) interval of
%                                      (xi-1)*B inside
%                                      [-speedConstancyOctaves, speedConstancyOctaves]
%                    probNearOne     - fraction of bootstrap xi with
%                                      |xi - 1| <= deltaSpeedTuned
%     .tfShiftOctaves      estimate (xi*B), ci, ci90
%     .speedChangeOctaves  estimate ((xi-1)*B), ci, ci90
%     .passbandOctaves     estimate (B), ci
%     .preferredSpeed      estimate (tf0/sf0, deg/s), ci
%     .quality       numBootstrap (valid samples), medianR2,
%                    fracSF0AtLowEdge, fracSF0AtHighEdge (NaN without sfRange)
%     .parameters    the thresholds used
%
%   See also: vis.speed.bootstrap.fit, vis.speed.bootstrap.xi_test

    arguments
        B (:,9) double
        pointEstimate double
        options.alpha (1,1) double {mustBeGreaterThan(options.alpha, 0), mustBeLessThan(options.alpha, 0.5)} = 0.05
        options.method (1,1) string {mustBeMember(options.method, ["percentile", "basic", "bc"])} = "percentile"
        options.deltaInseparable (1,:) double {mustBeNonnegative, mustBeNonempty} = 0.2
        options.deltaSpeedTuned (1,:) double {mustBeNonnegative, mustBeNonempty} = 0.5
        options.tfShiftOctaves (1,:) double {mustBeNonnegative, mustBeNonempty} = 1
        options.speedConstancyOctaves (1,:) double {mustBeNonnegative, mustBeNonempty} = 1
        options.min_xi (1,1) double = -4
        options.max_xi (1,1) double = 4
        options.sfRange double = []
    end

    pe = pointEstimate(:);
    if numel(pe) ~= 7
        error('vis:speed:bootstrap:indexes:badPointEstimate', 'POINTESTIMATE must have 7 elements.');
    end
    ok = all(~isnan(B(:, 1:7)), 2);
    B = B(ok, :);
    S = size(B, 1);
    a = options.alpha;
    bwFactor = 2 * sqrt(2 * log(2)) / log10(2);   % FWHM in octaves per unit sigma_sf

    xiStar = B(:, 3);
    bwStar = bwFactor * B(:, 4);
    xiHat = pe(3);
    bwHat = bwFactor * pe(4);

    % --- xi -----------------------------------------------------------------
    t2 = vis.speed.bootstrap.xi_test(xiStar, xiHat, [0 1], 'method', options.method, ...
        'confidenceLevel', 100 * (1 - a));
    t1 = vis.speed.bootstrap.xi_test(xiStar, xiHat, 0, 'method', options.method, ...
        'confidenceLevel', 100 * (1 - 2 * a));
    out.xi.estimate = xiHat;
    out.xi.bootstrapMedian = localMedian(xiStar);
    out.xi.se = localStd(xiStar);
    out.xi.ci = t2.ci;
    out.xi.ci90 = t1.ci;
    out.xi.ciWidth = t2.ci(2) - t2.ci(1);
    out.xi.p0 = t2.p(1);
    out.xi.p1 = t2.p(2);
    out.xi.probPositive = localFrac(xiStar > 0, S);
    out.xi.probNegative = localFrac(xiStar < 0, S);
    out.xi.fracAtBound = localFrac(abs(xiStar - options.min_xi) < 1e-3 | ...
        abs(xiStar - options.max_xi) < 1e-3, S);

    % --- derived quantities (percentile intervals over the bootstrap) -------
    shiftStar = xiStar .* bwStar;
    speedStar = (xiStar - 1) .* bwStar;
    out.tfShiftOctaves = localSummary(xiHat * bwHat, shiftStar, a);
    out.speedChangeOctaves = localSummary((xiHat - 1) * bwHat, speedStar, a);
    out.passbandOctaves = localSummary(bwHat, bwStar, a);
    out.preferredSpeed = localSummary(pe(7) / pe(6), B(:, 7) ./ B(:, 6), a);

    % --- inseparability ------------------------------------------------------
    dI = options.deltaInseparable;
    TI = options.tfShiftOctaves;
    out.inseparable.differentFrom0 = out.xi.p0 < a;
    out.inseparable.minimumEffect = localOutside(out.xi.ci90, dI);
    out.inseparable.passband = localOutside(out.tfShiftOctaves.ci90, TI);

    % --- speed tuning ---------------------------------------------------------
    dS = options.deltaSpeedTuned;
    TS = options.speedConstancyOctaves;
    ins = out.inseparable.differentFrom0;
    out.speedTuned.priebe = ins && out.xi.p1 >= a;
    out.speedTuned.equivalence = ins & localInside(out.xi.ci90, 1 - dS, 1 + dS);
    out.speedTuned.passband = ins & localInside(out.speedChangeOctaves.ci90, -TS, TS);
    out.speedTuned.probNearOne = arrayfun(@(d) localFrac(abs(xiStar - 1) <= d, S), dS);

    % --- quality --------------------------------------------------------------
    out.quality.numBootstrap = S;
    out.quality.medianR2 = localMedian(B(:, 9));
    if isempty(options.sfRange)
        out.quality.fracSF0AtLowEdge = NaN;
        out.quality.fracSF0AtHighEdge = NaN;
    else
        lo = min(options.sfRange); hi = max(options.sfRange);
        out.quality.fracSF0AtLowEdge = localFrac(abs(log10(B(:, 6)) - log10(lo)) < 1e-3, S);
        out.quality.fracSF0AtHighEdge = localFrac(abs(log10(B(:, 6)) - log10(hi)) < 1e-3, S);
    end

    out.parameters = struct('alpha', a, 'method', char(options.method), ...
        'deltaInseparable', dI, 'deltaSpeedTuned', dS, 'tfShiftOctaves', TI, ...
        'speedConstancyOctaves', TS, 'min_xi', options.min_xi, 'max_xi', options.max_xi, ...
        'sfRange', options.sfRange);

end

% -------------------------------------------------------------------------

function s = localSummary(estimate, samples, a)
    s.estimate = estimate;
    s.ci = localQuantile(samples, [a/2, 1 - a/2]);
    s.ci90 = localQuantile(samples, [a, 1 - a]);
end

function b = localOutside(ci, d)
% interval entirely above d or entirely below -d (d may be a vector)
    b = all(isfinite(ci)) & (ci(1) > d | ci(2) < -d);
end

function b = localInside(ci, lo, hi)
% interval inside [lo, hi] (lo, hi may be vectors)
    b = all(isfinite(ci)) & ci(1) >= lo & ci(2) <= hi;
end

function f = localFrac(tf, S)
    if S == 0
        f = NaN;
    else
        f = sum(tf) / S;
    end
end

function m = localMedian(x)
    x = x(~isnan(x));
    if isempty(x), m = NaN; else, m = median(x); end
end

function s = localStd(x)
    x = x(~isnan(x));
    if numel(x) < 2, s = NaN; else, s = std(x); end
end

function q = localQuantile(x, p)
    x = sort(x(~isnan(x)));
    n = numel(x);
    if n == 0
        q = nan(size(p));
        return;
    elseif n == 1
        q = repmat(x, size(p));
        return;
    end
    pos = 1 + (n - 1) * min(max(p, 0), 1);
    lo = floor(pos);
    hi = min(lo + 1, n);
    w = pos - lo;
    q = reshape(x(lo).' .* (1 - w) + x(hi).' .* w, size(p));
end
