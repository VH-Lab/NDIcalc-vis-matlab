function [parameters, sse, r_squared] = fit_fullspeed(sf, tf, r)
% FIT_FULLSPEED - Fit a 2D Gaussian with speed parameter (xi=1).
%
%  [PARAMETERS, SSE, R_SQUARED] = vis.speed.fit_fullspeed(SF, TF, R)
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
sf = sf(:); tf = tf(:); r = r(:); % Ensure column vectors
% --- Multi-start Optimization Setup ---
% Define hard constraints, identical to vis.speed.fit but with xi fixed at 1
Upper = [2 * max(r); 3; 1; 4; 4; max(sf); max(tf)];
Lower = [0; -3; 1; 0.1; 0.1; min(sf); min(tf)];
% Use the same systematic, semi-random multi-start approach as vis.speed.fit
% to ensure a fair comparison for the nested F-test.
zeta_rand = (Upper(2) - Lower(2)).*rand(10, 1) + Lower(2);
xi_rand = (Upper(3) - Lower(3)).*rand(10, 1) + Lower(3); % This will be all ones
sigma_sf_rand = (Upper(4) - Lower(4)).*rand(10, 1) + Lower(4);
sigma_tf_rand = (Upper(5) - Lower(5)).*rand(10, 1) + Lower(5);
StartPoint = [[0.2 * max(r); zeta_rand(1); xi_rand(1); sigma_sf_rand(1); sigma_tf_rand(1); 0.1 * max(sf); 0.1 * max(tf)], ...
    [0.4 * max(r); zeta_rand(2); xi_rand(2); sigma_sf_rand(2); sigma_tf_rand(2); 0.2 * max(sf); 0.2 * max(tf)], ...
    [0.6 * max(r); zeta_rand(3); xi_rand(3); sigma_sf_rand(3); sigma_tf_rand(3); 0.3 * max(sf); 0.3 * max(tf)], ...
    [0.8 * max(r); zeta_rand(4); xi_rand(4); sigma_sf_rand(4); sigma_tf_rand(4); 0.4 * max(sf); 0.4 * max(tf)], ...
    [max(r); zeta_rand(5); xi_rand(5); sigma_sf_rand(5); sigma_tf_rand(5); 0.5 * max(sf); 0.5 * max(tf)], ...
    [1.2 * max(r); zeta_rand(6); xi_rand(6); sigma_sf_rand(6); sigma_tf_rand(6); 0.6 * max(sf); 0.6 * max(tf)], ...
    [1.4 * max(r); zeta_rand(7); xi_rand(7); sigma_sf_rand(7); sigma_tf_rand(7); 0.7 * max(sf); 0.7 * max(tf)], ...
    [1.6 * max(r); zeta_rand(8); xi_rand(8); sigma_sf_rand(8); sigma_tf_rand(8); 0.8 * max(sf); 0.8 * max(tf)], ...
    [1.8 * max(r); zeta_rand(9); xi_rand(9); sigma_sf_rand(9); sigma_tf_rand(9); 0.9 * max(sf); 0.9 * max(tf)], ...
    [2 * max(r); zeta_rand(10); xi_rand(10); sigma_sf_rand(10); sigma_tf_rand(10); max(sf); max(tf)]];
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
