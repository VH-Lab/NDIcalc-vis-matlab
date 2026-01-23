function p_value = speed_nested_f(n_obs, rss_full, rss_reduced)
% SPEED.FIT.SPEED_NESTED_F  Perform F-test on Speed vs. No-Speed models.
%
%   P = SPEED.FIT.SPEED_NESTED_F(N_OBS, RSS_FULL, RSS_REDUCED)
%
%   Compares a full model (Speed, 7 params) against a reduced model
%   (No Speed, 6 params) using a nested F-test.
%
%   Inputs:
%     n_obs       : Number of data points (response vector length).
%     rss_full    : Residual Sum of Squares (SSE) of the full model (with speed).
%     rss_reduced : Residual Sum of Squares (SSE) of the reduced model (no speed).
%
%   Outputs:
%     p_value     : The probability of observing this F-statistic if the
%                   null hypothesis (reduced model) were true.

    % --- 1. Define Model Complexity ---
    p_full = 7;      % Parameters in the Speed model
    p_reduced = 6;   % Parameters in the No-Speed model

    % --- 2. Calculate Degrees of Freedom ---
    % df1: Difference in parameter counts (Numerator DF)
    df1 = p_full - p_reduced;

    % df2: Degrees of freedom of the residuals (Denominator DF)
    df2 = n_obs - p_full;

    % --- 3. Compute F-Statistic ---
    % F = (Improvement / df1) / (ErrorVariance / df2)

    delta_rss = rss_reduced - rss_full;

    % Sanity check: Full model should generally have lower error than reduced.
    if delta_rss < 0
        warning('Full model has higher error than reduced model. Check fits.');
    end

    numerator = delta_rss / df1;
    denominator = rss_full / df2;

    F_stat = numerator / denominator;

    % --- 4. Compute P-Value ---
    % Probability of observing a value > F_stat in the F-distribution
    p_value = 1 - fcdf(F_stat, df1, df2);

end
