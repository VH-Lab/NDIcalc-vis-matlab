classdef test_resample < matlab.unittest.TestCase
    % Tests for vis.speed.bootstrap.resample

    methods (Test)

        function testConstantTrialsReturnMean(testCase)
            % A condition whose trials are all equal has zero residuals, so
            % every resample returns its mean exactly.
            trials = [3 3 3 3 3; 7 7 7 7 7; -2 -2 -2 -2 -2];
            s = vis.randomstream('shuffle');
            for k = 1:20
                m = vis.speed.bootstrap.resample(trials, s);
                testCase.verifyEqual(m, [3; 7; -2], 'AbsTol', 1e-12);
            end
        end

        function testOutputSizeAndMissingTrials(testCase)
            % NaN-padded rows: a single trial is returned unchanged, no trials
            % gives NaN, and the output is nCond x 1.
            trials = [1 2 3 NaN NaN; 5 NaN NaN NaN NaN; NaN NaN NaN NaN NaN];
            m = vis.speed.bootstrap.resample(trials, vis.randomstream('shuffle'));
            testCase.verifySize(m, [3 1]);
            testCase.verifyEqual(m(2), 5);
            testCase.verifyTrue(isnan(m(3)));
            % Row 1 can only take values reachable from 3 scaled draws
            testCase.verifyGreaterThanOrEqual(m(1), 2 - sqrt(3/2) - 1e-12);
            testCase.verifyLessThanOrEqual(m(1), 2 + sqrt(3/2) + 1e-12);
        end

        function testMeanAndVarianceOfResamples(testCase)
            % Over many resamples, the resampled means are centered on the
            % condition mean, and (with scaled residuals) their variance is
            % the unbiased s^2/n; unscaled, it is (n-1)/n of that.
            x = [0 1 4 9 16];             % mean 6, unbiased var 42.5
            n = numel(x);
            N = 20000;
            s = vis.randomstream('shuffle');
            mScaled = zeros(N, 1);
            mRaw = zeros(N, 1);
            for k = 1:N
                mScaled(k) = vis.speed.bootstrap.resample(x, s);
                mRaw(k) = vis.speed.bootstrap.resample(x, s, 'scaleResiduals', false);
            end
            target = var(x) / n;          % 8.5
            testCase.verifyEqual(mean(mScaled), mean(x), 'AbsTol', 0.15);
            testCase.verifyEqual(var(mScaled), target, 'RelTol', 0.06);
            testCase.verifyEqual(var(mRaw), target * (n - 1) / n, 'RelTol', 0.06);
        end

        function testReproducibleWithSeed(testCase)
            % The same seed gives the same resample; the draws come only from
            % the stream passed in.
            trials = randn(10, 5);
            m1 = vis.speed.bootstrap.resample(trials, vis.randomstream(7));
            m2 = vis.speed.bootstrap.resample(trials, vis.randomstream(7));
            testCase.verifyEqual(m1, m2);
        end

    end
end
