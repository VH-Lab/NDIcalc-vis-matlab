classdef test_fit < matlab.unittest.TestCase
    % Tests for vis.speed.bootstrap.fit

    properties
        SF
        TF
    end

    methods (TestMethodSetup)
        function makeGrid(testCase)
            [S, T] = meshgrid([0.05 0.1 0.2 0.4 0.8], [1 2 4 8 16]);
            testCase.SF = S(:);
            testCase.TF = T(:);
        end
    end

    methods (Test)

        function testNoiselessRecoversParameters(testCase)
            % With identical trials every resample equals the true surface,
            % so every bootstrap fit recovers the true parameters.
            P = [10 0 0.5 0.25 0.4 0.2 4];
            mu = vis.speed.tuningfunc(testCase.SF, testCase.TF, P);
            trials = repmat(mu, 1, 5);
            S = 4;
            [B, info] = vis.speed.bootstrap.fit(testCase.SF, testCase.TF, trials, S, ...
                'PointEstimate', P(:), 'numberStartPoints', 2);
            testCase.verifySize(B, [S 9]);
            testCase.verifyEqual(numel(info.columns), 9);
            testCase.verifyEqual(info.columns{3}, 'xi');
            testCase.verifySize(info.resampledMeans, [numel(mu) S]);
            testCase.verifyEqual(B(:, 3), repmat(P(3), S, 1), 'AbsTol', 1e-2);
            testCase.verifyEqual(B(:, 6), repmat(P(6), S, 1), 'RelTol', 1e-2);
            testCase.verifyLessThan(B(:, 8), 1e-6);   % sse
        end

        function testNoisyCellCentersOnTruth(testCase)
            % Noisy cell with known xi: the bootstrap distribution of xi is
            % centered near the truth and has non-zero spread.
            P = [20 0 0.5 0.25 0.4 0.2 4];
            mu = vis.speed.tuningfunc(testCase.SF, testCase.TF, P);
            trials = mu + 1.0 * randn(numel(mu), 5);
            [B, info] = vis.speed.bootstrap.fit(testCase.SF, testCase.TF, trials, 30, ...
                'numberStartPoints', 3);
            xi = B(:, 3);
            testCase.verifyLessThan(abs(median(xi, 'omitnan') - P(3)), 0.25);
            testCase.verifyGreaterThan(std(xi, 'omitnan'), 0);
            testCase.verifyLessThan(abs(info.pointEstimate(3) - P(3)), 0.25);
        end

        function testResamplingReproducibleWithSeed(testCase)
            % The same RandomSeed draws the same resamples (the fits' own
            % start points stay unseeded, per AGENTS.md).
            P = [10 0 0 0.25 0.4 0.2 4];
            mu = vis.speed.tuningfunc(testCase.SF, testCase.TF, P);
            trials = mu + randn(numel(mu), 5);
            [~, i1] = vis.speed.bootstrap.fit(testCase.SF, testCase.TF, trials, 3, ...
                'PointEstimate', P(:), 'numberStartPoints', 1, 'RandomSeed', 11);
            [~, i2] = vis.speed.bootstrap.fit(testCase.SF, testCase.TF, trials, 3, ...
                'PointEstimate', P(:), 'numberStartPoints', 1, 'RandomSeed', 11);
            testCase.verifyEqual(i1.resampledMeans, i2.resampledMeans);
        end

        function testSizeMismatchErrors(testCase)
            testCase.verifyError(@() vis.speed.bootstrap.fit(testCase.SF, testCase.TF(1:end-1), ...
                ones(numel(testCase.SF), 5), 2), 'vis:speed:bootstrap:fit:sizeMismatch');
        end

    end
end
