classdef test_halfpeakcontour < matlab.unittest.TestCase
% TEST_HALFPEAKCONTOUR - Unit tests for vis.speed.halfpeakcontour

    methods (Test)

        function testOutputSizeAndPeak(testCase)
            % The contour must contain exactly 360 points (one per degree)
            % and report the peak at sf0/tf0.

            % P = [A, zeta, xi, sigma_sf, sigma_tf, sf0, tf0]
            P = [10, 0, 0, 0.5, 0.5, 0.1, 4];

            [X, Y, peak_sf, peak_tf] = vis.speed.halfpeakcontour(P);

            testCase.verifySize(X, [360 1], 'X must be a 360x1 column vector.');
            testCase.verifySize(Y, [360 1], 'Y must be a 360x1 column vector.');

            testCase.verifyEqual(peak_sf, P(6), 'peak_sf must equal sf0 (P(6)).');
            testCase.verifyEqual(peak_tf, P(7), 'peak_tf must equal tf0 (P(7)).');
        end

        function testResponseIsHalfPeak(testCase)
            % Evaluating the tuning function at each contour point should
            % return exactly half of the peak response.

            P = [10, 0, 0, 0.5, 0.5, 0.1, 4];

            [X, Y] = vis.speed.halfpeakcontour(P);

            % No NaNs expected for a well-behaved symmetric tuning curve.
            testCase.verifyFalse(any(isnan(X)), 'Unexpected NaN in X for a symmetric curve.');
            testCase.verifyFalse(any(isnan(Y)), 'Unexpected NaN in Y for a symmetric curve.');

            peak_response = vis.speed.tuningfunc(P(6), P(7), P);
            half_response = peak_response / 2;

            R_contour = vis.speed.tuningfunc(X(:), Y(:), P);

            testCase.verifyEqual(R_contour, repmat(half_response, 360, 1), ...
                'AbsTol', 1e-2, 'Response along the contour must be half the peak.');
        end

        function testSymmetricCaseIsCircleInLogSpace(testCase)
            % For an isotropic Gaussian (zeta=0, xi=0, sigma_sf=sigma_tf),
            % the half-peak contour is a circle in log10 space of radius
            % sigma*sqrt(2*ln2).

            sigma = 0.5;
            P = [10, 0, 0, sigma, sigma, 0.1, 4];

            [X, Y, peak_sf, peak_tf] = vis.speed.halfpeakcontour(P);

            center = [log10(peak_sf), log10(peak_tf)];
            log_radius = sqrt((log10(X) - center(1)).^2 + (log10(Y) - center(2)).^2);

            expected_radius = sigma * sqrt(2 * log(2));

            testCase.verifyEqual(log_radius, repmat(expected_radius, 360, 1), ...
                'AbsTol', 1e-3, 'Symmetric contour must be a circle in log10 space.');
        end

        function testContourEnclosesPeak(testCase)
            % The peak response must exceed the response everywhere on the
            % contour (the contour surrounds the peak).

            P = [8, 0.2, 0.7, 0.6, 0.4, 0.2, 2];

            [X, Y] = vis.speed.halfpeakcontour(P);

            peak_response = vis.speed.tuningfunc(P(6), P(7), P);

            valid = ~isnan(X);
            R_contour = vis.speed.tuningfunc(X(valid), Y(valid), P);

            testCase.verifyTrue(all(R_contour < peak_response), ...
                'All contour responses must be below the peak response.');
        end

        function testSkewedAndSpeedTunedCurve(testCase)
            % A curve with non-zero skew (zeta) and speed (xi) should still
            % produce a finite half-peak contour whose responses are half
            % the peak.

            P = [5, 0.15, 1, 0.5, 0.5, 0.1, 4];

            [X, Y] = vis.speed.halfpeakcontour(P);

            valid = ~isnan(X) & ~isnan(Y);
            testCase.verifyTrue(all(valid), 'Expected a complete contour for this curve.');

            peak_response = vis.speed.tuningfunc(P(6), P(7), P);
            half_response = peak_response / 2;

            R_contour = vis.speed.tuningfunc(X, Y, P);
            testCase.verifyEqual(R_contour, repmat(half_response, 360, 1), ...
                'AbsTol', 1e-2, 'Response along the skewed contour must be half the peak.');
        end

        function testCustomRadiusOptions(testCase)
            % The MaxRadius and RadiusStep options should be accepted and
            % produce a valid contour.

            P = [10, 0, 0, 0.5, 0.5, 0.1, 4];

            [X, Y] = vis.speed.halfpeakcontour(P, 'MaxRadius', 6, 'RadiusStep', 5e-3);

            testCase.verifySize(X, [360 1], 'X size must remain 360x1.');
            testCase.verifySize(Y, [360 1], 'Y size must remain 360x1.');

            R_contour = vis.speed.tuningfunc(X, Y, P);
            half_response = vis.speed.tuningfunc(P(6), P(7), P) / 2;
            testCase.verifyEqual(R_contour, repmat(half_response, 360, 1), ...
                'AbsTol', 5e-2, 'Coarser step should still recover the half-peak level.');
        end

        function testInvalidParameterLength(testCase)
            % A parameter vector that is not length 7 must error.

            badP = [10, 0, 0, 0.5, 0.5, 0.1]; % only 6 elements

            testCase.verifyError(@() vis.speed.halfpeakcontour(badP), ...
                'vis:speed:halfpeakcontour:badParameters');
        end

        function testNonPositivePeakFrequencies(testCase)
            % Non-positive preferred frequencies must error.

            P = [10, 0, 0, 0.5, 0.5, 0, 4]; % sf0 = 0

            testCase.verifyError(@() vis.speed.halfpeakcontour(P), ...
                'vis:speed:halfpeakcontour:badPeak');
        end

    end
end
