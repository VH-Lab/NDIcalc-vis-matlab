classdef test_fit < matlab.unittest.TestCase

    methods (Test)

        function testFit(testCase)
            % Test vis.speed.fit

            % Generate data
            [SF, TF] = meshgrid([0.05 0.1 0.2], [1 2 4]);
            sf = SF(:); tf = TF(:);

            % True params: xi = 0.5 (between 0 and 1)
            true_params = [10; 0; 0.5; 1; 1; 0.1; 2];
            r = vis.speed.tuningfunc(sf, tf, true_params);

            % Fit
            % usage: fit(SF, TF, R, min_xi, max_xi)
            [params, sse, r2] = vis.speed.fit(sf, tf, r, 0, 1);

            % Verify
            testCase.verifyEqual(params(3), 0.5, 'RelTol', 0.1, 'xi should be close to 0.5');
            testCase.verifyGreaterThan(r2, 0.9, 'R-squared should be high');
        end

        function testFitNospeed(testCase)
            % Test vis.speed.fit_nospeed (forces xi=0)

            % Generate data with xi=0
            [SF, TF] = meshgrid([0.05 0.1 0.2], [1 2 4]);
            sf = SF(:); tf = TF(:);

            true_params = [10; 0; 0; 1; 1; 0.1; 2];
            r = vis.speed.tuningfunc(sf, tf, true_params);

            [params, sse, r2] = vis.speed.fit_nospeed(sf, tf, r);

            % Verify xi is 0
            testCase.verifyEqual(params(3), 0, 'AbsTol', 1e-10, 'xi must be 0 for fit_nospeed');
            testCase.verifyGreaterThan(r2, 0.9, 'R-squared should be high');
        end

    end
end
