function [B, info] = fit(sf, tf, trials, numBootstrap, options)
% FIT - bootstrap the free speed-tuning fit over single trials
%
%   [B, INFO] = vis.speed.bootstrap.fit(SF, TF, TRIALS, NUMBOOTSTRAP, ...)
%
%   Repeats the free Priebe et al. (2006) fit (vis.speed.fit, xi free in
%   [min_xi, max_xi]) on NUMBOOTSTRAP bootstrap resamples of the condition
%   means (vis.speed.bootstrap.resample), and returns every fit.
%
%   Each resample is exactly one call to vis.speed.fit, so the bootstrap fits
%   are the same fit as the point estimate. Starting points follow
%   vis.speed.fit_priebe_triplet (ndi.calc.vis.speed_tuning): the observed
%   means are first fit with xi fixed at 0 (vis.speed.fit_nospeed) and at 1
%   (vis.speed.fit_fullspeed), and those two fits are fixed starts of the
%   point estimate, in addition to its 'pointEstimateStartPoints' random
%   starts. Each bootstrap fit starts from the point estimate and the two
%   nested fits, in addition to its 'numberStartPoints' random starts, so
%   every resample tries the xi = 0 and xi = 1 solutions as well as the
%   observed data's own.
%
%   All resamples are drawn first, serially, from one private random stream,
%   so the result does not depend on whether the fits then run serially or in
%   parallel ('useParallel').
%
%   Inputs:
%     SF, TF        - spatial (cyc/deg) and temporal (Hz) frequency of each
%                     condition (nCond values each).
%     TRIALS        - nCond x nTrials single-trial responses (NaN padded),
%                     rows matched to SF/TF.
%     NUMBOOTSTRAP  - number of bootstrap resamples S.
%
%   Name-Value Pair Arguments:
%     min_xi, max_xi    - range of xi for the free fit. (Default: -4, 4; use
%                         a range that contains the values to be tested, e.g.
%                         0 and 1, away from its edges)
%     PointEstimate     - 7x1 parameters of the fit to the observed means. If
%                         empty (default), it is computed here with
%                         vis.speed.fit, 'pointEstimateStartPoints' random
%                         starts and the two nested fits as fixed starts.
%     pointEstimateStartPoints - random starts for the point estimate. (40)
%     numberStartPoints - random starts per bootstrap fit, in addition to the
%                         point estimate and the two nested fits. (Default: 5;
%                         use 40 to match the point estimate)
%     scaleResiduals    - passed to vis.speed.bootstrap.resample. (true)
%     RandomSeed        - seed for the resampling stream (see
%                         vis.randomstream). (Default: 'shuffle')
%     useParallel       - use parfor (with whatever pool is open). (false)
%
%   Outputs:
%     B    - NUMBOOTSTRAP x 9 matrix, one row per resample:
%              columns 1-7: [A zeta xi sigma_sf sigma_tf sf0 tf0]
%              column 8:    sum of squared errors of that fit
%              column 9:    r-squared of that fit
%            A row is NaN if that fit failed.
%     INFO - struct with fields:
%              columns        - names of the 9 columns of B
%              pointEstimate  - 7x1 parameters of the observed-data fit
%              pointEstimateSSE, pointEstimateR2
%              nestedStarts   - 7x2 fits of the observed means with xi = 0
%                               and xi = 1 (fixed starts of every fit; a
%                               column is dropped if that fit failed)
%              resampledMeans - nCond x NUMBOOTSTRAP resampled condition means
%              min_xi, max_xi, numberStartPoints, scaleResiduals
%
%   Example:
%     [SF, TF] = meshgrid([0.05 0.1 0.2 0.4 0.8], [1 2 4 8 16]);
%     mu = vis.speed.tuningfunc(SF(:), TF(:), [10 0 0.5 0.25 0.4 0.2 4]);
%     trials = mu + 2*randn(numel(mu), 5);
%     [B, info] = vis.speed.bootstrap.fit(SF(:), TF(:), trials, 200);
%     histogram(B(:,3));    % bootstrap distribution of xi
%
%   See also: vis.speed.fit, vis.speed.bootstrap.resample,
%             vis.speed.bootstrap.xi_test, vis.speed.bootstrap.classify

    arguments
        sf double
        tf double
        trials double
        numBootstrap (1,1) double {mustBeInteger, mustBePositive}
        options.min_xi (1,1) double = -4
        options.max_xi (1,1) double = 4
        options.PointEstimate double = []
        options.pointEstimateStartPoints (1,1) double {mustBeInteger, mustBePositive} = 40
        options.numberStartPoints (1,1) double {mustBeInteger, mustBePositive} = 5
        options.scaleResiduals (1,1) logical = true
        options.RandomSeed = 'shuffle'
        options.useParallel (1,1) logical = false
    end

    sf = sf(:); tf = tf(:);
    nCond = numel(sf);
    if numel(tf) ~= nCond || size(trials, 1) ~= nCond
        error('vis:speed:bootstrap:fit:sizeMismatch', ...
            'SF, TF and the rows of TRIALS must all have one entry per condition.');
    end

    observedMean = mean(trials, 2, 'omitnan');

    % --- Nested fits (xi = 0 and xi = 1), as in vis.speed.fit_priebe_triplet;
    % they are fixed starting points of the point estimate and of every
    % bootstrap fit.
    pNo = vis.speed.fit_nospeed(sf, tf, observedMean, ...
        'numberStartPoints', options.pointEstimateStartPoints);
    pFull = vis.speed.fit_fullspeed(sf, tf, observedMean, ...
        'numberStartPoints', options.pointEstimateStartPoints);
    nestedStarts = zeros(7, 0);
    for p = {pNo, pFull}
        if numel(p{1}) == 7 && all(isfinite(p{1}(:)))
            nestedStarts(:, end+1) = p{1}(:); %#ok<AGROW>
        end
    end
    % keep xi inside the bounds of the free fit
    nestedStarts(3, :) = min(max(nestedStarts(3, :), options.min_xi), options.max_xi);

    % --- Point estimate (the regular fit to the observed means) -----------
    if isempty(options.PointEstimate)
        [pe, peSSE, peR2] = vis.speed.fit(sf, tf, observedMean, options.min_xi, options.max_xi, ...
            'numberStartPoints', options.pointEstimateStartPoints, ...
            'SpecificStartPoint', nestedStarts);
        if numel(pe) ~= 7 || any(~isfinite(pe(:)))
            error('vis:speed:bootstrap:fit:pointEstimateFailed', ...
                ['The fit to the observed means did not converge from any of its %d ' ...
                'starting points, so there is no point estimate to bootstrap.'], ...
                options.pointEstimateStartPoints);
        end
    else
        pe = options.PointEstimate(:);
        if numel(pe) ~= 7 || any(~isfinite(pe))
            error('vis:speed:bootstrap:fit:badPointEstimate', 'PointEstimate must have 7 finite elements.');
        end
        R = vis.speed.tuningfunc(sf, tf, pe);
        peSSE = sum((observedMean - R).^2);
        sst = sum((observedMean - mean(observedMean)).^2);
        peR2 = 1 - peSSE / max(sst, eps);
    end
    % Keep the starting point inside the bounds of the bootstrap fits.
    peStart = pe;
    peStart(3) = min(max(peStart(3), options.min_xi), options.max_xi);
    fixedStarts = [peStart nestedStarts];

    % --- Draw every resample first, serially, from one private stream ----
    stream = vis.randomstream(options.RandomSeed);
    resampledMeans = nan(nCond, numBootstrap);
    for s = 1:numBootstrap
        resampledMeans(:, s) = vis.speed.bootstrap.resample(trials, stream, ...
            'scaleResiduals', options.scaleResiduals);
    end

    % --- One free fit per resample ------------------------------------------
    B = nan(numBootstrap, 9);
    min_xi = options.min_xi;
    max_xi = options.max_xi;
    nStarts = options.numberStartPoints;
    if options.useParallel
        nWorkers = Inf;
    else
        nWorkers = 0;
    end
    parfor (s = 1:numBootstrap, nWorkers)
        B(s, :) = localFitOne(sf, tf, resampledMeans(:, s), min_xi, max_xi, nStarts, fixedStarts);
    end

    info.columns = {'A', 'zeta', 'xi', 'sigma_sf', 'sigma_tf', 'sf0', 'tf0', 'sse', 'r_squared'};
    info.pointEstimate = pe;
    info.pointEstimateSSE = peSSE;
    info.pointEstimateR2 = peR2;
    info.nestedStarts = nestedStarts;
    info.resampledMeans = resampledMeans;
    info.min_xi = options.min_xi;
    info.max_xi = options.max_xi;
    info.numberStartPoints = options.numberStartPoints;
    info.scaleResiduals = options.scaleResiduals;

end

function row = localFitOne(sf, tf, r, min_xi, max_xi, nStarts, fixedStarts)
% One free fit of one resample; NaN row on failure.
    row = nan(1, 9);
    ok = ~isnan(r);
    if sum(ok) < 8
        return;
    end
    try
        [p, sse, r2] = vis.speed.fit(sf(ok), tf(ok), r(ok), min_xi, max_xi, ...
            'numberStartPoints', nStarts, 'SpecificStartPoint', fixedStarts);
        % A fit with no usable start returns no parameters (or NaN); keep
        % the NaN row so one bad resample cannot stop the whole run.
        if numel(p) == 7 && all(isfinite(p(:))) && isscalar(sse) && isscalar(r2)
            row = [p(:).' sse r2];
        end
    catch
        row = nan(1, 9);
    end
end
