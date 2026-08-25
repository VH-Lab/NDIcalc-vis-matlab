function [parameters, sse, r_squared] = fit_fullspeed(sf, tf, r, options)
% FIT_FULLSPEED - Fit a 2D Gaussian with speed parameter (xi=1).
%
%  [PARAMETERS, SSE, R_SQUARED] = vis.speed.fit_fullspeed(SF, TF, R, ...)
%
%  Fits a Gaussian to a set of responses where the speed parameter 'xi'
%  is constrained to be 1. This function uses a harmonized multi-start
%  approach that is identical to `vis.speed.fit` to ensure a fair comparison
%  for nested F-tests.
%
%  Inputs:
%    sf - An array of spatial frequency values.
%    tf - An array of temporal frequency values.
%    r  - An array of measured responses.
%
%  Name-Value Pair Arguments:
%     numberStartPoints  - The number of initial starting points to try. (Default: 40)
%     SpecificStartPoint - A 7xN matrix specifying specific starting
%                          points to include in the search.
%     RandomSeed         - Seed for the random components of the starting
%                          points. (Default: 'shuffle', a fresh draw on every
%                          call.) The draws use a private stream, so they do
%                          not disturb the caller's global random stream.
%                          Pass a non-negative integer to pin the search to
%                          one reproducible set of starts, for instance when
%                          reproducing a particular run while debugging.
%                          See vis.randomstream.
%
%  Outputs:
%    PARAMETERS - A 7x1 vector with the best-fit parameters.
%    SSE        - The total sum of squared errors for the best fit.
%    R_SQUARED  - The R-squared value (coefficient of determination).
%
%  Parameters:
%    [A, zeta, xi, sigma_sf, sigma_tf, sf0, tf0], where xi is always 1.
%
%  See also: vis.speed.fit, vis.speed.fit_nospeed, lsqcurvefit
%

    arguments
        sf
        tf
        r
        options.numberStartPoints (1,1) double {mustBeInteger, mustBePositive} = 40
        options.SpecificStartPoint (7,:) double = []
        options.RandomSeed = 'shuffle'
    end

    sf = sf(:); tf = tf(:); r = r(:); % Ensure column vectors

    % --- Multi-start Optimization Setup ---
    % Define hard constraints, identical to vis.speed.fit but with xi fixed at 1.
    % sigma_sf, sigma_tf, and zeta are scaled by log(10) so that the search
    % space matches the prior natural-log formulation now that tuningfunc
    % uses log10.
    Upper = [2 * max(r); 3 * log(10); 1; 4 / log(10); 4 / log(10); max(sf); max(tf)];
    Lower = [0; -3 * log(10); 1; 0.1 / log(10); 0.1 / log(10); min(sf); min(tf)];

    StartPoint = zeros(7, options.numberStartPoints);

    % The random components of the start points are drawn fresh on every
    % call, from a private stream that leaves the caller's global stream
    % undisturbed. They are deliberately not seeded: the search must reach
    % the same optimum whatever starts it is given, and a run that only
    % agrees with a stored answer for one pinned set of starts would be
    % evidence of nothing. Pass RandomSeed to pin them while debugging.
    startStream = vis.randomstream(options.RandomSeed);

    for i = 1:options.numberStartPoints
        rule_idx = mod(i-1, 10) + 1;

        % Generate random components
        zeta_val = (Upper(2) - Lower(2)) * rand(startStream) + Lower(2);
        % xi_val is fixed to 1 (Upper(3)=1, Lower(3)=1)
        xi_val = (Upper(3) - Lower(3)) * rand(startStream) + Lower(3);
        sigma_sf_val = (Upper(4) - Lower(4)) * rand(startStream) + Lower(4);
        sigma_tf_val = (Upper(5) - Lower(5)) * rand(startStream) + Lower(5);

        A_val = (0.2 * rule_idx) * max(r);
        sf0_val = (0.1 * rule_idx) * max(sf);
        tf0_val = (0.1 * rule_idx) * max(tf);

        StartPoint(:, i) = [A_val; zeta_val; xi_val; sigma_sf_val; sigma_tf_val; sf0_val; tf0_val];
    end

    if ~isempty(options.SpecificStartPoint)
        StartPoint = [StartPoint, options.SpecificStartPoint];
    end

    % --- Fitting Loop ---
    optimal_error = Inf;
    P_optimal = [];
    opts = optimset('Display', 'off');

    for i = 1:size(StartPoint, 2)
        [P_current, error_current] = lsqcurvefit(@(P, SF) vis.speed.tuningfunc(SF, tf, P), StartPoint(:, i), sf, r, Lower, Upper, opts);
        if (error_current < optimal_error)
            optimal_error = error_current;
            P_optimal = P_current;
        end
    end

    parameters = P_optimal;
    sse = optimal_error;

    % --- Calculate Goodness of Fit ---
    sst = sum((r - mean(r)).^2);
    if sst == 0
        r_squared = 1;
    else
        r_squared = 1 - (sse / sst);
    end
end
