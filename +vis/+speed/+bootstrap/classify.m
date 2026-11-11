function out = classify(xi_boot, xi_hat, options)
% CLASSIFY - classify spatiotemporal coupling from bootstrap estimates of xi
%
%   OUT = vis.speed.bootstrap.classify(XI_BOOT, XI_HAT, ...)
%
%   Classifies a site from the bootstrap distribution of its coupling index xi
%   (column 3 of vis.speed.bootstrap.fit) using vis.speed.bootstrap.xi_test:
%
%     inseparable - xi differs from 0: p(xi == 0) < alpha.
%     speedTuned  - EQUIVALENCE criterion: inseparable AND the
%                   (1 - 2*alpha) confidence interval lies entirely inside
%                   [1 - delta, 1 + delta]. A poorly constrained site has a
%                   wide interval and cannot pass, unlike the "fails to reject
%                   xi == 1" rule. (Using the 1 - 2*alpha interval is the
%                   standard two one-sided tests (TOST) construction at level
%                   alpha.)
%     speedTunedFailToReject - the older rule, for comparison: inseparable
%                   AND p(xi == 1) >= alpha.
%
%   Name-Value Pair Arguments:
%     method - 'percentile' (default), 'basic', or 'bc' (see xi_test)
%     alpha  - significance level (default 0.05)
%     delta  - equivalence margin around xi = 1 (default 0.5)
%
%   Output OUT, a struct:
%     inseparable, speedTuned, speedTunedFailToReject (logical)
%     p0, p1          - p-values for xi == 0 and xi == 1
%     ci              - (1 - alpha) confidence interval for xi
%     ciEquivalence   - (1 - 2*alpha) interval used for the equivalence test
%     method, alpha, delta
%
%   See also: vis.speed.bootstrap.xi_test, vis.speed.bootstrap.fit

    arguments
        xi_boot double
        xi_hat (1,1) double
        options.method (1,1) string {mustBeMember(options.method, ["percentile", "basic", "bc"])} = "percentile"
        options.alpha (1,1) double {mustBeGreaterThan(options.alpha, 0), mustBeLessThan(options.alpha, 0.5)} = 0.05
        options.delta (1,1) double {mustBeNonnegative} = 0.5
    end

    t = vis.speed.bootstrap.xi_test(xi_boot, xi_hat, [0 1], ...
        'method', options.method, 'confidenceLevel', 100 * (1 - options.alpha));
    tEq = vis.speed.bootstrap.xi_test(xi_boot, xi_hat, 1, ...
        'method', options.method, 'confidenceLevel', 100 * (1 - 2 * options.alpha));

    out.p0 = t.p(1);
    out.p1 = t.p(2);
    out.ci = t.ci;
    out.ciEquivalence = tEq.ci;
    out.inseparable = out.p0 < options.alpha;
    out.speedTuned = out.inseparable && ...
        tEq.ci(1) >= 1 - options.delta && tEq.ci(2) <= 1 + options.delta;
    out.speedTunedFailToReject = out.inseparable && out.p1 >= options.alpha;
    out.method = char(options.method);
    out.alpha = options.alpha;
    out.delta = options.delta;

end
