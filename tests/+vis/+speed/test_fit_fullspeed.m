classdef test_fit_fullspeed < matlab.unittest.TestCase

    methods (Test)

        function testFitFullSpeed(testCase)
            % Test that vis.speed.fit_fullspeed correctly constrains xi to 1
            % and finds reasonable parameters for synthetic data.

            % 1. Generate synthetic data with xi = 1
            sf_vals = logspace(log10(0.01), log10(1), 5);
            tf_vals = logspace(log10(0.5), log10(16), 5);
            [SF, TF] = meshgrid(sf_vals, tf_vals);
            sf = SF(:);
            tf = TF(:);

            % True parameters: [A, zeta, xi, sigma_sf, sigma_tf, sf0, tf0]
            true_xi = 1;
            true_params = [10; 0; true_xi; 1; 1; 0.1; 4];

            % Generate response
            r = vis.speed.tuningfunc(sf, tf, true_params);

            % Add a tiny bit of noise to avoid perfect fit issues/warnings,
            % but keep it small to verify accuracy.
            rng(1);
            r_noisy = r + 0.01 * randn(size(r));
            r_noisy(r_noisy<0) = 0;

            % 2. Call fit_fullspeed
            [fitted_params, sse, r_squared] = vis.speed.fit_fullspeed(sf, tf, r_noisy);

            % 3. Verifications

            % Check that xi is exactly 1
            xi_idx = 3;
            testCase.verifyEqual(fitted_params(xi_idx), 1, 'AbsTol', 1e-10, ...
                'The fitted speed parameter xi must be exactly 1.');

            % Check that R-squared is high (the fit should be good)
            testCase.verifyGreaterThan(r_squared, 0.95, ...
                'The R-squared value should be high for data generated with xi=1.');

            % Check that other parameters are reasonably close to true parameters
            % (Allowing some tolerance due to noise and optimization landscape)
            % Parameter vector: [A, zeta, xi, sigma_sf, sigma_tf, sf0, tf0]

            % Verify A (Amplitude)
            testCase.verifyEqual(fitted_params(1), true_params(1), 'RelTol', 0.2, 'Amplitude A estimation failed.');

            % Verify sf0 (Spatial Frequency Center)
            testCase.verifyEqual(fitted_params(6), true_params(6), 'RelTol', 0.2, 'sf0 estimation failed.');

            % Verify tf0 (Temporal Frequency Center)
            testCase.verifyEqual(fitted_params(7), true_params(7), 'RelTol', 0.2, 'tf0 estimation failed.');

        end

        function testFitFullSpeedStructure(testCase)
             % Test that the output structure and types are correct
             % We need at least 7 data points to avoid the warning about trust-region-reflective
             % requiring at least as many equations as variables.

             sf = [0.1; 0.2; 0.4; 0.8; 0.1; 0.2; 0.4; 0.8];
             tf = [1; 1; 1; 1; 2; 2; 2; 2];
             r = [5; 10; 5; 2; 6; 11; 6; 3]; % Dummy data, 8 points

             [params, sse, rsq] = vis.speed.fit_fullspeed(sf, tf, r);

             testCase.verifySize(params, [7 1], 'Parameters must be a 7x1 vector.');
             testCase.verifyClass(sse, 'double', 'SSE must be a double.');
             testCase.verifyClass(rsq, 'double', 'R-squared must be a double.');
             testCase.verifyEqual(params(3), 1, 'xi must be 1 even for dummy data.');
        end

    end
end
