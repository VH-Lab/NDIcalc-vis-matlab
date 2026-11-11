classdef test_xi_test < matlab.unittest.TestCase
    % Tests for vis.speed.bootstrap.xi_test, using constructed bootstrap
    % distributions (no fitting).

    properties
        X      % deterministic "bootstrap" sample: normal quantiles, mean 0.5, sd 0.1
    end

    methods (TestMethodSetup)
        function makeSample(testCase)
            S = 1999;
            q = ((1:S).' - 0.5) / S;
            testCase.X = 0.5 + 0.1 * (-sqrt(2) * erfcinv(2 * q));
        end
    end

    methods (Test)

        function testPercentile(testCase)
            x = testCase.X;
            out = vis.speed.bootstrap.xi_test(x, 0.5, [0.5 0 0.5 + 0.1*1.959964]);
            testCase.verifyEqual(out.p(1), 1, 'AbsTol', 1e-3);        % at the center
            testCase.verifyEqual(out.p(2), 2 / (numel(x) + 1), 'AbsTol', 1e-12); % far outside
            testCase.verifyEqual(out.p(3), 0.05, 'AbsTol', 0.005);    % at the 97.5th percentile
            testCase.verifyEqual(out.ci, [0.5 - 0.196 0.5 + 0.196], 'AbsTol', 0.003);
            testCase.verifyEqual(out.numBootstrap, numel(x));
        end

        function testBasicSymmetric(testCase)
            % For a symmetric distribution centered on xi_hat, basic = percentile
            x = testCase.X;
            outB = vis.speed.bootstrap.xi_test(x, 0.5, 0.5 + 0.1*1.959964, 'method', 'basic');
            outP = vis.speed.bootstrap.xi_test(x, 0.5, 0.5 + 0.1*1.959964, 'method', 'percentile');
            testCase.verifyEqual(outB.p, outP.p, 'AbsTol', 0.005);
            testCase.verifyEqual(outB.ci, outP.ci, 'AbsTol', 0.003);
        end

        function testBasicShifted(testCase)
            % Basic bootstrap reflects the spread about xi_hat: a sample
            % centered at 0.5 with xi_hat = 0.7 puts the interval around 0.9,
            % and 0.9 (the center of that interval) is not rejected.
            out = vis.speed.bootstrap.xi_test(testCase.X, 0.7, 0.9, 'method', 'basic');
            testCase.verifyEqual(mean(out.ci), 0.9, 'AbsTol', 0.003);
            testCase.verifyEqual(out.p, 1, 'AbsTol', 1e-3);
        end

        function testBiasCorrected(testCase)
            % With xi_hat at the median, z0 = 0 and BC equals percentile.
            x = testCase.X;
            outBC = vis.speed.bootstrap.xi_test(x, 0.5, [0.3 0.5 + 0.1*1.959964], 'method', 'bc');
            outP = vis.speed.bootstrap.xi_test(x, 0.5, [0.3 0.5 + 0.1*1.959964]);
            testCase.verifyEqual(outBC.ci, outP.ci, 'AbsTol', 0.003);
            testCase.verifyEqual(outBC.p(2), 0.05, 'AbsTol', 0.005);
            testCase.verifyEqual(outBC.p(1), 0.0455, 'AbsTol', 0.005);  % 2 sd below
            outFar = vis.speed.bootstrap.xi_test(x, 0.5, 0, 'method', 'bc');
            testCase.verifyLessThan(outFar.p, 0.01);
            % A point estimate above the median shifts the BC interval upward
            outShift = vis.speed.bootstrap.xi_test(x, 0.55, 0.5, 'method', 'bc');
            testCase.verifyGreaterThan(outShift.ci(1), outP.ci(1));
        end

        function testNaNsIgnoredAndEmpty(testCase)
            x = [testCase.X; NaN(10, 1)];
            out = vis.speed.bootstrap.xi_test(x, 0.5, 0);
            testCase.verifyEqual(out.numBootstrap, numel(testCase.X));
            outEmpty = vis.speed.bootstrap.xi_test(NaN(5, 1), 0.5, 0);
            testCase.verifyTrue(isnan(outEmpty.p));
            testCase.verifyTrue(all(isnan(outEmpty.ci)));
        end

        function testConfidenceLevel(testCase)
            out90 = vis.speed.bootstrap.xi_test(testCase.X, 0.5, 0, 'confidenceLevel', 90);
            testCase.verifyEqual(out90.ci, [0.5 - 0.1645 0.5 + 0.1645], 'AbsTol', 0.003);
        end

    end
end
