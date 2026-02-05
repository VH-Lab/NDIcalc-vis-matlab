classdef test_tuningfunc < matlab.unittest.TestCase

    methods (Test)

        function testTuningFunc(testCase)
            % Test that vis.speed.tuningfunc calculates responses correctly

            % Define inputs
            [SF, TF] = meshgrid([0.05 0.1 0.2], [1 2 4]);

            % Parameters: [A, zeta, xi, sigma_sf, sigma_tf, sf0, tf0]
            % Case 1: Simple Gaussian at sf0, tf0 (xi=0, zeta=0)
            A = 10;
            zeta = 0;
            xi = 0;
            sigma_sf = 1;
            sigma_tf = 1;
            sf0 = 0.1;
            tf0 = 2;

            P = [A, zeta, xi, sigma_sf, sigma_tf, sf0, tf0];

            R = vis.speed.tuningfunc(SF, TF, P);

            % Verification
            testCase.verifySize(R, size(SF), 'Output size must match input meshgrid size.');

            % Check peak response location
            % Ideally at sf=0.1, tf=2, the response should be roughly A * exp(0) * (exp(0) - exp(...))
            % Let's calculate one point manually

            sf_val = 0.1;
            tf_val = 2;
            logtfpsf = xi * (log10(sf_val)-log10(sf0)) + log10(tf0); % should be log10(2)

            term1 = exp( (-(log10(sf_val)-log10(sf0)).^2) / (2*sigma_sf*sigma_sf) ); % exp(0) = 1
            term2 = exp( -(log10(tf_val)-logtfpsf).^2 ./ (2.*(sigma_tf+zeta*(log10(tf_val)-logtfpsf)).^2) ); % exp(0) = 1
            term3 = exp(-1/(zeta.^2)); % For zeta=0, this term is tricky.

            % Wait, looking at the code:
            % exp(-1/(zeta.^2))
            % If zeta is 0, this is exp(-Inf) = 0.

            % Let's verify what the code actually does when zeta is 0.
            % MATLAB handles division by zero -> Inf. exp(-Inf) -> 0.
            % So the offset term is 0.

            expected_peak = A * 1 * (1 - 0); % Should be A = 10.

            % Find index in R where SF=0.1 and TF=2
            idx = find(SF==0.1 & TF==2);
            testCase.verifyEqual(R(idx), expected_peak, 'AbsTol', 1e-6, 'Peak response calculation incorrect.');

            % Case 2: Check shapes with a non-zero xi and zeta
            zeta = 0.1;
            xi = 1;
            P2 = [A, zeta, xi, sigma_sf, sigma_tf, sf0, tf0];
            R2 = vis.speed.tuningfunc(SF, TF, P2);

            testCase.verifySize(R2, size(SF), 'Output size incorrect for P2.');
            testCase.verifyTrue(all(isfinite(R2(:))), 'Output contains non-finite values.');
        end

        function testInputShapes(testCase)
             % Verify it handles column vectors vs matrices
             sf = [0.1; 0.2];
             tf = [1; 2];
             P = [10, 0, 0, 1, 1, 0.1, 2];

             R = vis.speed.tuningfunc(sf, tf, P);
             testCase.verifySize(R, [2 1], 'Column vector input should return column vector output.');
        end
    end
end
