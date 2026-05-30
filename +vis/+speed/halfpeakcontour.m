function [X, Y, peak_sf, peak_tf] = halfpeakcontour(P, options)
% HALFPEAKCONTOUR - Compute the half-peak contour of a speed tuning function
%
%   [X, Y] = vis.speed.halfpeakcontour(P)
%   [X, Y, PEAK_SF, PEAK_TF] = vis.speed.halfpeakcontour(P, ...)
%
%   Given a set of parameters P for a speed tuning function (see
%   vis.speed.tuningfunc), compute the coordinates of the contour at which
%   the response falls to half of its peak value (the half-peak, or
%   half-maximum, contour) surrounding the preferred spatial and temporal
%   frequency.
%
%   The contour is computed by sweeping rays outward from the peak in
%   log10(SF)/log10(TF) space, which is the space in which the tuning
%   function is naturally symmetric. One point is returned for each whole
%   degree of angle, so X and Y are 360x1 column vectors. X contains the
%   spatial frequency (cycles/degree) coordinates and Y contains the
%   temporal frequency (Hz) coordinates of the contour.
%
%   Inputs:
%     P  - A 7-element parameter vector [A, zeta, xi, sigma_sf, sigma_tf,
%          sf0, tf0] as used by vis.speed.tuningfunc.
%
%   Name-Value Pair Arguments:
%     MaxRadius  - The maximum distance (in log10 units) to search outward
%                  from the peak along each ray. (Default: 10)
%     RadiusStep - The step size (in log10 units) used when searching for
%                  the half-peak crossing along each ray. (Default: 1e-3)
%
%   Outputs:
%     X        - 360x1 vector of spatial frequency coordinates (c/deg).
%     Y        - 360x1 vector of temporal frequency coordinates (Hz).
%     PEAK_SF  - The preferred spatial frequency (cycles/degree), equal to
%                the sf0 parameter, P(6).
%     PEAK_TF  - The preferred temporal frequency (Hz), equal to the tf0
%                parameter, P(7).
%
%   The i-th point (X(i),Y(i)) lies along the ray at angle (i-1) degrees,
%   measured counter-clockwise in log10(SF)/log10(TF) space, where 0 degrees
%   points toward increasing spatial frequency. If the response never falls
%   to half its peak along a ray within MaxRadius, the corresponding
%   coordinates are returned as NaN.
%
%   Example:
%     P = [10, 0, 0, 0.5, 0.5, 0.1, 4]; % A simple symmetric tuning curve
%     [X, Y] = vis.speed.halfpeakcontour(P);
%     figure;
%     loglog(X, Y, 'k-'); hold on;
%     loglog(P(6), P(7), 'r+'); % mark the peak
%     xlabel('Spatial frequency (c/deg)');
%     ylabel('Temporal frequency (Hz)');
%
%   See also: vis.speed.tuningfunc, vis.speed.plottuning
%

    arguments
        P (1,:) double {mustBeReal}
        options.MaxRadius (1,1) double {mustBePositive} = 10
        options.RadiusStep (1,1) double {mustBePositive} = 1e-3
    end

    if numel(P) ~= 7
        error('vis:speed:halfpeakcontour:badParameters', ...
            'P must be a 7-element parameter vector [A zeta xi sigma_sf sigma_tf sf0 tf0].');
    end

    peak_sf = P(6); % sf0, the preferred spatial frequency
    peak_tf = P(7); % tf0, the preferred temporal frequency

    if peak_sf <= 0 || peak_tf <= 0
        error('vis:speed:halfpeakcontour:badPeak', ...
            'The preferred frequencies sf0 (P(6)) and tf0 (P(7)) must be positive.');
    end

    % Work in log10 space, where the tuning function is naturally symmetric.
    center = [log10(peak_sf), log10(peak_tf)];

    % The peak response is attained at (sf0, tf0).
    peak_response = vis.speed.tuningfunc(peak_sf, peak_tf, P);
    half_response = peak_response / 2;

    % Set up the radial search grid (shared across all directions).
    r = (0:options.RadiusStep:options.MaxRadius).';
    nAngles = 360;

    X = nan(nAngles, 1);
    Y = nan(nAngles, 1);

    if ~(peak_response > 0)
        % A non-positive peak has no meaningful half-peak contour.
        return;
    end

    for i = 1:nAngles
        theta = (i - 1) * pi / 180; % angle in radians, one point per degree
        dir = [cos(theta), sin(theta)];

        % Coordinates along this ray in log10 space.
        logSF = center(1) + r * dir(1);
        logTF = center(2) + r * dir(2);

        Rvals = vis.speed.tuningfunc(10.^logSF, 10.^logTF, P);

        % Find the first point along the ray at or below the half-peak level.
        below = find(Rvals <= half_response, 1, 'first');

        if isempty(below) || below == 1
            % No crossing found (or the peak itself is already <= half): leave NaN.
            continue;
        end

        % Linearly interpolate (in radius) between the bracketing samples to
        % refine the location of the half-peak crossing.
        rLow = r(below - 1);
        rHigh = r(below);
        Rlow = Rvals(below - 1);
        Rhigh = Rvals(below);

        if Rhigh == Rlow
            rCross = rHigh;
        else
            rCross = rLow + (half_response - Rlow) * (rHigh - rLow) / (Rhigh - Rlow);
        end

        X(i) = 10.^(center(1) + rCross * dir(1));
        Y(i) = 10.^(center(2) + rCross * dir(2));
    end

end
