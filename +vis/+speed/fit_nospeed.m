function [parameters, sse, r_squared] = fit_nospeed(sf, tf, r)
% FIT_NOSPEED Fit a two-dimensional Gaussian function with no speed parameter.
%
%  [PARAMETERS, SSE, R_SQUARED] = FIT_NOSPEED(SF, TF, R)
%
%  Fits a Gaussian to a set of responses where the speed parameter 'xi'
%  is constrained to be 0.
%
%  Inputs:
%    SF - An array of spatial frequency values.
%    TF - An array of temporal frequency values.
%    R  - An array of measured responses.
%
%  Outputs:
%    PARAMETERS - An array with the best-fit parameters.
%    SSE        - The total sum of squared errors for the fit.
%    R_SQUARED  - The R-squared value (coefficient of determination),
%                 indicating the goodness of fit.
%
%  Parameters:
%    |-----------|-------------------------------|
%    | A         | Peak response of the neuron   |
%    | zeta      | Skew of temporal freq. curve  |
%    | xi        | Speed parameter (fixed at 0)  |
%    | sigma_sf  | Tuning width (spatial freq.)  |
%    | sigma_tf  | Tuning width (temporal freq.) |
%    | sf0       | Preferred spatial frequency   |
%    | tf0       | Preferred temporal frequency  |
%    |-----------|-------------------------------|
%
%  See also: lsqcurvefit
%
    random_starts = 10;

    sf = sf(:); tf = tf(:); r = r(:); % Make them column vectors

    % Define constraints on the fit
    Upper = [2 * max(r); 3; 0; 4; 4; max(sf); max(tf)]; % xi is fixed at 0
    Lower = [0; -3; 0; 0.1; 0.1; min(sf); min(tf)];     % xi is fixed at 0

    % Use lsqcurvefit with a multi-start optimization
    optimal_error = Inf;
    P_optimal = NaN(7, 1);
    opts = optimset('Display','off'); % Suppress solver output

    for i = 1:random_starts
        % Generate random starting points within the bounds for each parameter
        StartPoint = Lower + rand(size(Lower)).*(Upper - Lower);
        StartPoint(3) = 0; % Ensure xi starts at 0

        % Perform the curve fit for this starting point
        [current_params, current_error] = lsqcurvefit(@(P, SF) vis.speed.tuningfunc(SF, tf, P), StartPoint, sf, r, Lower, Upper, opts);

        % Keep the best result found so far
        if (current_error < optimal_error)
            optimal_error = current_error;
            P_optimal = current_params;
        end
    end
    
    parameters = P_optimal;
    sse = optimal_error;

    % Calculate R-squared
    sst = sum((r - mean(r)).^2);
    if sst == 0  % Handle the case of zero variance in the data
        r_squared = 1; % If sse is also 0, this is a perfect fit
    else
        r_squared = 1 - (sse / sst);
    end

end
