classdef test_fit_nospeed < matlab.unittest.TestCase

    methods (Test)

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
