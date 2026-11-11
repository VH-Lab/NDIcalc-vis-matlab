classdef test_indexes < matlab.unittest.TestCase
    % Tests for vis.speed.bootstrap.indexes, using constructed bootstrap
    % matrices (no fitting): xi follows exact normal quantiles, the other
    % parameters are fixed.

    methods (Static)
        function B = makeB(xiCenter, xiSD, sigma_sf)
            % S x 9 bootstrap matrix with the given xi distribution and a
            % fixed SF tuning width (sigma_sf, log10 units).
            S = 1999;
            q = ((1:S).' - 0.5) / S;
            xi = xiCenter + xiSD * (-sqrt(2) * erfcinv(2 * q));
            o = ones(S, 1);
            B = [10*o, 0*o, xi, sigma_sf*o, 0.4*o, 0.2*o, 4*o, 1*o, 0.9*o];
        end
        function pe = makePE(xi, sigma_sf)
            pe = [10 0 xi sigma_sf 0.4 0.2 4];
        end
    end

    properties (Constant)
        SIG2 = 2 * log10(2) / (2 * sqrt(2 * log(2)));   % sigma_sf giving a 2-octave passband
    end

    methods (Test)

        function testSeparable(testCase)
            s = testCase.SIG2;
            out = vis.speed.bootstrap.indexes(testCase.makeB(0.02, 0.1, s), testCase.makePE(0.02, s));
            testCase.verifyFalse(out.inseparable.differentFrom0);
            testCase.verifyFalse(out.inseparable.minimumEffect);
            testCase.verifyFalse(out.inseparable.passband);
            testCase.verifyFalse(out.speedTuned.priebe);
            testCase.verifyFalse(out.speedTuned.equivalence);
            testCase.verifyFalse(out.speedTuned.passband);
            testCase.verifyEqual(out.passbandOctaves.estimate, 2, 'AbsTol', 1e-10);
        end

        function testWellConstrainedSpeedTuned(testCase)
            s = testCase.SIG2;
            out = vis.speed.bootstrap.indexes(testCase.makeB(1, 0.1, s), testCase.makePE(1, s));
            testCase.verifyTrue(out.inseparable.differentFrom0);
            testCase.verifyTrue(out.speedTuned.priebe);
            testCase.verifyTrue(out.speedTuned.equivalence);
            testCase.verifyTrue(out.speedTuned.passband);
            testCase.verifyGreaterThan(out.speedTuned.probNearOne, 0.99);
            testCase.verifyEqual(out.preferredSpeed.estimate, 20, 'AbsTol', 1e-10);
        end

        function testPoorlyConstrainedOnlyPriebe(testCase)
            % A wide distribution excluding 0 and including 1 passes the
            % Priebe rule only.
            s = testCase.SIG2;
            out = vis.speed.bootstrap.indexes(testCase.makeB(1.4, 0.6, s), testCase.makePE(1.4, s));
            testCase.verifyTrue(out.speedTuned.priebe);
            testCase.verifyFalse(out.speedTuned.equivalence);
            testCase.verifyFalse(out.speedTuned.passband);
            testCase.verifyGreaterThan(out.xi.ciWidth, 2);
        end

        function testIntermediateCoupling(testCase)
            % xi ~ 0.5 (tight), 2-octave passband: preferred TF shifts ~1 octave.
            s = testCase.SIG2;
            B = testCase.makeB(0.5, 0.08, s);
            pe = testCase.makePE(0.5, s);
            out = vis.speed.bootstrap.indexes(B, pe);
            testCase.verifyTrue(out.inseparable.differentFrom0);
            testCase.verifyTrue(out.inseparable.minimumEffect);            % |xi| > 0.2
            testCase.verifyFalse(out.inseparable.passband);                % shift ~1 octave, not > 1
            testCase.verifyTrue(vis.speed.bootstrap.indexes(B, pe, 'tfShiftOctaves', 0.5).inseparable.passband);
            testCase.verifyFalse(out.speedTuned.priebe);
            testCase.verifyEqual(out.tfShiftOctaves.estimate, 1, 'AbsTol', 1e-10);
            testCase.verifyEqual(out.speedChangeOctaves.estimate, -1, 'AbsTol', 1e-10);
        end

        function testSmallCouplingFailsMinimumEffect(testCase)
            s = testCase.SIG2;
            out = vis.speed.bootstrap.indexes(testCase.makeB(0.15, 0.03, s), testCase.makePE(0.15, s));
            testCase.verifyTrue(out.inseparable.differentFrom0);
            testCase.verifyFalse(out.inseparable.minimumEffect);
            testCase.verifyTrue(vis.speed.bootstrap.indexes(testCase.makeB(0.15, 0.03, s), ...
                testCase.makePE(0.15, s), 'deltaInseparable', 0.05).inseparable.minimumEffect);
        end

        function testNarrowPassbandSpeedConstancy(testCase)
            % With a 1-octave passband, xi = 0.3 keeps preferred speed within
            % 0.7 octaves: passband speed-tuned, though far from xi = 1.
            s = testCase.SIG2 / 2;
            out = vis.speed.bootstrap.indexes(testCase.makeB(0.3, 0.05, s), testCase.makePE(0.3, s));
            testCase.verifyEqual(out.passbandOctaves.estimate, 1, 'AbsTol', 1e-10);
            testCase.verifyTrue(out.speedTuned.passband);
            testCase.verifyFalse(out.speedTuned.equivalence);
            testCase.verifyFalse(out.speedTuned.priebe);
        end

        function testQualityFlagsAndNaNs(testCase)
            s = testCase.SIG2;
            B = testCase.makeB(0.5, 0.1, s);
            B(1:100, 3) = 4;                       % samples at the xi bound
            B(1:200, 6) = 0.04;                    % sf0 at the lowest stimulus SF
            B(end-9:end, :) = NaN;                 % failed fits
            out = vis.speed.bootstrap.indexes(B, testCase.makePE(0.5, s), 'sfRange', [0.04 1.25]);
            nValid = size(B, 1) - 10;
            testCase.verifyEqual(out.quality.numBootstrap, nValid);
            testCase.verifyEqual(out.xi.fracAtBound, 100 / nValid, 'AbsTol', 1e-12);
            testCase.verifyEqual(out.quality.fracSF0AtLowEdge, 200 / nValid, 'AbsTol', 1e-12);
            testCase.verifyEqual(out.quality.fracSF0AtHighEdge, 0);
            outNoRange = vis.speed.bootstrap.indexes(B, testCase.makePE(0.5, s));
            testCase.verifyTrue(isnan(outNoRange.quality.fracSF0AtLowEdge));
        end

        function testThresholdVectors(testCase)
            % Vectors of thresholds give one classification per value.
            s = testCase.SIG2;
            out = vis.speed.bootstrap.indexes(testCase.makeB(0.15, 0.03, s), testCase.makePE(0.15, s), ...
                'deltaInseparable', [0.05 0.2], 'deltaSpeedTuned', [0.25 0.5 1], ...
                'tfShiftOctaves', [0.1 1], 'speedConstancyOctaves', [0.5 2]);
            testCase.verifyEqual(out.inseparable.minimumEffect, [true false]);
            testCase.verifyEqual(out.inseparable.passband, [true false]);   % shift ~0.3 octaves
            testCase.verifySize(out.speedTuned.equivalence, [1 3]);
            testCase.verifyEqual(out.speedTuned.equivalence, [false false true]);   % 90% CI ~[0.10 0.20]
            testCase.verifyEqual(out.speedTuned.passband, [false true]);     % speed change ~-1.7 octaves
            testCase.verifySize(out.speedTuned.probNearOne, [1 3]);
        end

        function testParametersRecorded(testCase)
            s = testCase.SIG2;
            out = vis.speed.bootstrap.indexes(testCase.makeB(0, 0.1, s), testCase.makePE(0, s), ...
                'deltaInseparable', 0.3, 'deltaSpeedTuned', 0.25, 'tfShiftOctaves', 2, ...
                'speedConstancyOctaves', 0.5, 'alpha', 0.01);
            testCase.verifyEqual(out.parameters.deltaInseparable, 0.3);
            testCase.verifyEqual(out.parameters.deltaSpeedTuned, 0.25);
            testCase.verifyEqual(out.parameters.tfShiftOctaves, 2);
            testCase.verifyEqual(out.parameters.speedConstancyOctaves, 0.5);
            testCase.verifyEqual(out.parameters.alpha, 0.01);
        end

    end
end
