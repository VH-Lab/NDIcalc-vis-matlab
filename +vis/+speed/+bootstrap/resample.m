function meanResp = resample(trials, stream, options)
% RESAMPLE - one bootstrap resample of the per-condition mean responses
%
%   MEANRESP = vis.speed.bootstrap.resample(TRIALS, STREAM, ...)
%
%   Draws one nonparametric bootstrap resample of a tuning curve's condition
%   means from its single-trial responses.
%
%   For each condition k with n_k >= 2 trials, the residuals about the
%   condition mean are scaled by sqrt(n_k/(n_k-1)) and n_k of them are drawn
%   with replacement:
%       MEANRESP(k) = mean_k + mean(drawn scaled residuals).
%   Resampling n trials with replacement underestimates the variance of a mean
%   by the factor (n-1)/n (0.8 for n = 5); the scaling removes that bias, so
%   the resampled means have, on average, the unbiased sample variance s^2/n.
%   Set 'scaleResiduals' to false for the classical (unscaled) bootstrap,
%   which is identical to redrawing the raw trials with replacement.
%
%   A condition with a single trial carries no information about its noise and
%   returns that trial unchanged; a condition with no trials returns NaN.
%
%   Inputs:
%     TRIALS - nCond x nTrials matrix of single-trial responses (one row per
%              stimulus condition; NaN where a condition has fewer trials).
%     STREAM - a RandStream (e.g. from vis.randomstream) used for every draw.
%
%   Name-Value Pair Arguments:
%     scaleResiduals - true (default) to scale residuals by sqrt(n/(n-1))
%
%   Output:
%     MEANRESP - nCond x 1 resampled condition means.
%
%   Example:
%     trials = [1 2 3 4 5; 10 10 10 10 10];
%     s = vis.randomstream(1);
%     m = vis.speed.bootstrap.resample(trials, s);   % m(2) is always 10
%
%   See also: vis.speed.bootstrap.fit, vis.randomstream

    arguments
        trials double
        stream (1,1) RandStream
        options.scaleResiduals (1,1) logical = true
    end

    nCond = size(trials, 1);
    meanResp = nan(nCond, 1);
    for k = 1:nCond
        x = trials(k, :);
        x = x(~isnan(x));
        n = numel(x);
        if n == 0
            continue;
        end
        m = mean(x);
        if n == 1
            meanResp(k) = m;
            continue;
        end
        r = x - m;
        if options.scaleResiduals
            r = r * sqrt(n / (n - 1));
        end
        draws = r(randi(stream, n, n, 1));
        meanResp(k) = m + mean(draws);
    end

end
