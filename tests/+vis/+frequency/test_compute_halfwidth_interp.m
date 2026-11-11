classdef test_compute_halfwidth_interp < matlab.unittest.TestCase

    methods (Test)

        function testSmoothPeak(testCase)
            % A Gaussian in log frequency peaked at 1 cycle/degree.
            sf = logspace(log10(0.05), log10(1.2), 9) * 2;
            y = 10 * exp(-(log2(sf) - 0).^2 / 2);

            [low, maxv, high] = vis.frequency.compute_halfwidth_interp(sf, y);

            [~, imax] = max(y);
            testCase.verifyEqual(maxv, sf(imax));
            testCase.verifyLessThan(low, maxv);
            testCase.verifyGreaterThan(high, maxv);
            % linear interpolation puts each point exactly at half height
            testCase.verifyEqual(interp1(sf, y, low), max(y)/2, 'AbsTol', 1e-9);
            testCase.verifyEqual(interp1(sf, y, high), max(y)/2, 'AbsTol', 1e-9);
        end

        function testStopsAtFirstCrossingFromPeak(testCase)
            % After falling below half height the response rises above it
            % again; the walk from the peak must stop at the first crossing.
            x = [0.05 0.1 0.2 0.3 0.4 0.6 0.8 1.0];
            y = [  4   1   3  10   2   4.9 0.5 4.8];

            [low, maxv, high] = vis.frequency.compute_halfwidth_interp(x, y);

            testCase.verifyEqual(maxv, 0.3);
            testCase.verifyEqual(low, 0.2 + (5 - 3) / (10 - 3) * 0.1, 'AbsTol', 1e-12);
            testCase.verifyEqual(high, 0.3 + (5 - 10) / (2 - 10) * 0.1, 'AbsTol', 1e-12);
        end

        function testNaNWhenNotReached(testCase)
            % Low-pass: never falls to half on the low side.
            x = [0.05 0.1 0.2 0.4 0.8];
            y = [10 9 8 3 1];

            [low, maxv, high] = vis.frequency.compute_halfwidth_interp(x, y);

            testCase.verifyEqual(maxv, 0.05);
            testCase.verifyTrue(isnan(low));
            testCase.verifyTrue(isfinite(high));
        end

    end
end
