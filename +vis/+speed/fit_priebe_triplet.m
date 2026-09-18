function [f, f_no_speed, f_fullspeed, stats] = fit_priebe_triplet(sf, tf, r, min_xi, max_xi, options)
% FIT_PRIEBE_TRIPLET - fit the three nested Priebe speed-tuning models
%
%  [F, F_NO_SPEED, F_FULLSPEED, STATS] = ...
%       vis.speed.fit_priebe_triplet(SF, TF, R, MIN_XI, MAX_XI, ...)
%
%  Runs the three fits the nested F test for speed tuning requires, exactly as
%  ndi.calc.vis.speed_tuning does, so the two speed-tuning calculators cannot
%  drift in how they fit:
%    * F           - the free fit, xi free in [MIN_XI, MAX_XI]  (vis.speed.fit)
%    * F_NO_SPEED  - xi constrained to 0               (vis.speed.fit_nospeed)
%    * F_FULLSPEED - xi constrained to 1              (vis.speed.fit_fullspeed)
%  Each returned parameter vector is the 7x1 Priebe vector
%  [A, zeta, xi, sigma_sf, sigma_tf, sf0, tf0]. The constrained arms carry xi
%  pinned to 0 or 1 in element 3.
%
%  The constrained fits are computed first and handed to the free fit as its
%  SpecificStartPoint, so the free fit never does worse than either constrained
%  arm -- the property the nested F test relies on.
%
%  STATS is a struct:
%    .sse, .sse_no_speed, .sse_fullspeed
%    .r_squared, .r_squared_no_speed, .r_squared_fullspeed
%    .partial_r2_no_speed, .partial_r2_fullspeed
%    .nested_F_no_speed_p_value, .nested_F_fullspeed_p_value
%
%  Name/value options:
%    FitRandomSeed - forwarded to vis.speed.fit / fit_nospeed / fit_fullspeed
%                    as their RandomSeed. Default [] (do not forward, so each
%                    fit uses its own private 'shuffle' stream). Per AGENTS.md,
%                    self-tests must not pin this; it exists so a real analysis
%                    can be reproduced when a caller chooses to.
%
%  See also: vis.speed.fit, vis.speed.fit_nospeed, vis.speed.fit_fullspeed,
%            vis.speed.speed_nested_f, ndi.calc.vis.speed_tuning,
%            ndi.calc.vis.speed_tuning_bootstrap

    arguments
        sf (:,1) double
        tf (:,1) double
        r  (:,1) double
        min_xi (1,1) double = 0
        max_xi (1,1) double = 1
        options.FitRandomSeed = []
    end

    fitArgs = {};
    if ~isempty(options.FitRandomSeed)
        fitArgs = {'RandomSeed', options.FitRandomSeed};
    end

    % fit with the speed parameter set to 0
    [f_no_speed,  sse_no_speed,  r2_no_speed]  = vis.speed.fit_nospeed(sf, tf, r, fitArgs{:});
    % fit with the speed parameter set to 1 (full speed)
    [f_fullspeed, sse_fullspeed, r2_fullspeed] = vis.speed.fit_fullspeed(sf, tf, r, fitArgs{:});
    % free fit, seeded from the two constrained fits
    [f, sse, r_squared] = vis.speed.fit(sf, tf, r, min_xi, max_xi, ...
        'SpecificStartPoint', [f_no_speed f_fullspeed], fitArgs{:});

    if sse_no_speed ~= 0
        partial_r2_no_speed = (sse_no_speed - sse) / sse_no_speed;
    else
        partial_r2_no_speed = 0;
    end
    if sse_fullspeed ~= 0
        partial_r2_fullspeed = (sse_fullspeed - sse) / sse_fullspeed;
    else
        partial_r2_fullspeed = 0;
    end

    num_responses = numel(r);

    stats = struct( ...
        'sse',                        sse, ...
        'sse_no_speed',               sse_no_speed, ...
        'sse_fullspeed',              sse_fullspeed, ...
        'r_squared',                  r_squared, ...
        'r_squared_no_speed',         r2_no_speed, ...
        'r_squared_fullspeed',        r2_fullspeed, ...
        'partial_r2_no_speed',        partial_r2_no_speed, ...
        'partial_r2_fullspeed',       partial_r2_fullspeed, ...
        'nested_F_no_speed_p_value',  vis.speed.speed_nested_f(num_responses, sse, sse_no_speed), ...
        'nested_F_fullspeed_p_value', vis.speed.speed_nested_f(num_responses, sse, sse_fullspeed));

end % fit_priebe_triplet()
