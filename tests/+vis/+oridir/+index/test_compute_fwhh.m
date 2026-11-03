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

        function testStopsAtFirstCrossingFromPeak(testCase)
            % After falling below half height the response rises above it
            % again; the walk from the peak must stop at the first crossing.
            x = 0:10;
            y = [0 2 9 1 3 10 4 1 8 2 0];

            [fwhh, low, high] = vis.oridir.index.compute_fwhh(x, y);

            % half height is 5: crossings between x=4 (3) and x=5 (10),
            % and between x=5 (10) and x=6 (4)
            testCase.verifyEqual(low, 4 + (5 - 3) / (10 - 3), 'AbsTol', 1e-12);
            testCase.verifyEqual(high, 5 + (5 - 10) / (4 - 10), 'AbsTol', 1e-12);
            testCase.verifyEqual(fwhh, high - low, 'AbsTol', 1e-12);
        end

        function testWrapsAroundWithPeriod(testCase)
            % A peak at 0 degrees has its lower half-height point near 360.
            angles = 0:30:330;
            sig = 25;
            shape = @(pref) 10 * exp(-vlt.math.angdiff(angles - pref).^2 / (2 * sig^2));

            [fwhh0, low0, high0] = vis.oridir.index.compute_fwhh(angles, shape(0), 'Period', 360);
            fwhh180 = vis.oridir.index.compute_fwhh(angles, shape(180), 'Period', 360);

            testCase.verifyTrue(isfinite(fwhh0));
            testCase.verifyEqual(fwhh0, fwhh180, 'AbsTol', 1e-9);
            testCase.verifyGreaterThan(low0, 300);
            testCase.verifyLessThan(high0, 60);
            testCase.verifyEqual(mod(high0 - low0, 360), fwhh0, 'AbsTol', 1e-9);
        end

        function testPeakAtLastSampleWithPeriod(testCase)
            % A peak at the last sample, 330 degrees, must be able to walk
            % past it on the high side, between 330 and 360.
            angles = 0:30:330;
            sig = 25;
            shape = @(pref) 10 * exp(-vlt.math.angdiff(angles - pref).^2 / (2 * sig^2));

            [fwhh, low, high] = vis.oridir.index.compute_fwhh(angles, shape(330), 'Period', 360);
            fwhh180 = vis.oridir.index.compute_fwhh(angles, shape(180), 'Period', 360);

            testCase.verifyEqual(fwhh, fwhh180, 'AbsTol', 1e-9);
            testCase.verifyGreaterThan(low, 270);
            testCase.verifyLessThan(low, 330);
            testCase.verifyGreaterThan(high, 330);
        end

        function testNoWrapByDefault(testCase)
            % Without a Period, a peak at the first sample has no low point.
            angles = 0:30:330;
            y = 10 * exp(-vlt.math.angdiff(angles).^2 / (2 * 25^2));

            [fwhh, low] = vis.oridir.index.compute_fwhh(angles, y);

            testCase.verifyEqual(fwhh, Inf);
            testCase.verifyTrue(isnan(low));
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
