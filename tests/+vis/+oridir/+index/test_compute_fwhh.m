classdef test_compute_fwhh < matlab.unittest.TestCase

    methods (Test)

        function testGaussianWidth(testCase)
            % A Gaussian peaked at 4 with exp(-(x-4)^2/4) reaches half height
            % at 4 +/- sqrt(4*ln(2)), so its full width is 2*sqrt(4*ln(2)).
            x = linspace(1, 10, 100);
            y = exp(-((x - 4).^2) / 4);

            [fwhh, low, high] = vis.oridir.index.compute_fwhh(x, y);

            expected_low = 4 - sqrt(4 * log(2));
            expected_high = 4 + sqrt(4 * log(2));
            testCase.verifyEqual(low, expected_low, 'AbsTol', 0.05);
            testCase.verifyEqual(high, expected_high, 'AbsTol', 0.05);
            testCase.verifyEqual(fwhh, expected_high - expected_low, 'AbsTol', 0.1);
            testCase.verifyEqual(fwhh, high - low, 'AbsTol', 1e-12);
        end

        function testWidthDoesNotDependOnPeakLocation(testCase)
            % The same tuning shape at two preferred directions has the same width.
            angles = 0:30:330;
            shape = @(pref) 10 * exp(-(angles - pref).^2 / (2 * 25^2));

            fwhh_a = vis.oridir.index.compute_fwhh(angles, shape(90));
            fwhh_b = vis.oridir.index.compute_fwhh(angles, shape(240));

            testCase.verifyTrue(isfinite(fwhh_a));
            testCase.verifyEqual(fwhh_a, fwhh_b, 'AbsTol', 1);
        end

        function testIgnoresNullDirectionPeak(testCase)
            % A direction-tuned curve whose null response is half the peak
            % must still be measured at the crossings next to the peak.
            angles = 0:10:350;
            sig = 20;
            y = 20 * exp(-(angles - 90).^2 / (2 * sig^2)) + ...
                10 * exp(-(angles - 270).^2 / (2 * sig^2));

            [fwhh, low, high] = vis.oridir.index.compute_fwhh(angles, y);

            expected_half_width = sig * sqrt(2 * log(2));
            testCase.verifyEqual(low, 90 - expected_half_width, 'AbsTol', 1);
            testCase.verifyEqual(high, 90 + expected_half_width, 'AbsTol', 1);
            testCase.verifyEqual(fwhh, 2 * expected_half_width, 'AbsTol', 2);
        end

        function testInfiniteWhenNeverDropsBelowHalfMax(testCase)
            x = 1:10;
            y = [5 6 7 8 9 10 9 8 7 6];

            [fwhh, low, high] = vis.oridir.index.compute_fwhh(x, y);

            testCase.verifyEqual(fwhh, Inf);
            testCase.verifyTrue(isnan(low) || isnan(high));
        end

    end
end
