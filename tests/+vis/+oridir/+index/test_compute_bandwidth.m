classdef test_compute_bandwidth < matlab.unittest.TestCase

    methods (Test)

        function testFiniteBandwidth(testCase)
            % Test computing bandwidth for a Gaussian-like response curve
            % where low and high half-max points are well defined and > 0.
            x = linspace(1, 10, 100);
            % Peak at x = 4, dropping below half-max on both sides (e.g. at x=2 and x=8)
            y = exp(-((x - 4).^2) / 4);

            bw = vis.oridir.index.compute_bandwidth(x, y);

            testCase.verifyTrue(isfinite(bw), 'Bandwidth should be finite for a well-behaved peak.');
            testCase.verifyGreaterThan(bw, 0, 'Bandwidth in octaves should be positive.');

            % Expected low and high can be estimated from y = exp(-(x-4)^2 / 4) = 0.5
            % (x-4)^2 / 4 = ln(2) => x-4 = +/- sqrt(4 * ln(2)) ~= +/- 1.6651
            % low ~= 2.3349, high ~= 5.6651
            % expected bw = log2(5.6651 / 2.3349) ~= 1.2788
            expected_low = 4 - sqrt(4 * log(2));
            expected_high = 4 + sqrt(4 * log(2));
            expected_bw = log2(expected_high / expected_low);

            testCase.verifyEqual(bw, expected_bw, 'AbsTol', 0.05, ...
                'Bandwidth should match log2(high/low) calculated from half-max cutoff points.');
        end

        function testInfiniteBandwidthWhenNeverDropsBelowHalfMax(testCase)
            % Test when response curve does not fall below half maximum on one or both sides.
            x = 1:10;
            y = [5 6 7 8 9 10 9 8 7 6]; % min is 5, max is 10, half-max is 5, so min never drops below 5

            bw = vis.oridir.index.compute_bandwidth(x, y);

            testCase.verifyEqual(bw, Inf, 'Bandwidth should be Inf when response does not drop below half max.');
        end

        function testInfiniteBandwidthWhenFlatResponse(testCase)
            % Test flat response curve
            x = 1:10;
            y = ones(1, 10);

            bw = vis.oridir.index.compute_bandwidth(x, y);

            testCase.verifyEqual(bw, Inf, 'Bandwidth should be Inf for flat response curve.');
        end

    end
end
