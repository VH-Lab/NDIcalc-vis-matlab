classdef test_classify < matlab.unittest.TestCase
    % Tests for vis.speed.bootstrap.classify, using constructed bootstrap
    % distributions of xi (no fitting).

    methods (Static)
        function x = sample(center, sd)
            S = 1999;
            q = ((1:S).' - 0.5) / S;
            x = center + sd * (-sqrt(2) * erfcinv(2 * q));
        end
    end

    methods (Test)

        function testSeparable(testCase)
            x = testCase.sample(0.02, 0.1);
            out = vis.speed.bootstrap.classify(x, 0.02);
            testCase.verifyFalse(out.inseparable);
            testCase.verifyFalse(out.speedTuned);
            testCase.verifyFalse(out.speedTunedFailToReject);
        end

        function testWellConstrainedSpeedTuned(testCase)
            x = testCase.sample(1.0, 0.1);
            out = vis.speed.bootstrap.classify(x, 1.0);
            testCase.verifyTrue(out.inseparable);
            testCase.verifyTrue(out.speedTuned);
            testCase.verifyTrue(out.speedTunedFailToReject);
        end

        function testPoorlyConstrainedIsNotSpeedTuned(testCase)
            % Wide distribution that excludes 0 but is consistent with 1:
            % the old rule calls it speed-tuned, the equivalence rule does not.
            x = testCase.sample(1.4, 0.6);
            out = vis.speed.bootstrap.classify(x, 1.4);
            testCase.verifyTrue(out.inseparable);
            testCase.verifyTrue(out.speedTunedFailToReject);
            testCase.verifyFalse(out.speedTuned);
        end

        function testIntermediateCoupling(testCase)
            % xi ~ 0.5, tight: inseparable, differs from 1, not speed-tuned
            x = testCase.sample(0.5, 0.08);
            out = vis.speed.bootstrap.classify(x, 0.5);
            testCase.verifyTrue(out.inseparable);
            testCase.verifyLessThan(out.p1, 0.05);
            testCase.verifyFalse(out.speedTuned);
        end

        function testDeltaControlsEquivalence(testCase)
            x = testCase.sample(1.0, 0.2);   % 90% CI about [0.67 1.33]
            testCase.verifyTrue(vis.speed.bootstrap.classify(x, 1, 'delta', 0.5).speedTuned);
            testCase.verifyFalse(vis.speed.bootstrap.classify(x, 1, 'delta', 0.2).speedTuned);
        end

    end
end
