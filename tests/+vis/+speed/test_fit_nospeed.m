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

        function testFitNospeedOptions(testCase)
            % Test options for fit_nospeed

            [SF, TF] = meshgrid([0.05 0.1 0.2], [1 2 4]);
            sf = SF(:); tf = TF(:);
            true_params = [10; 0; 0; 1; 1; 0.1; 2];
            r = vis.speed.tuningfunc(sf, tf, true_params);

            % Test numberStartPoints
            [params, ~, ~] = vis.speed.fit_nospeed(sf, tf, r, 'numberStartPoints', 10);
            testCase.verifyEqual(params(3), 0, 'AbsTol', 1e-10, 'xi must be 0');

             % Test SpecificStartPoint
            startPt = [9; 0.1; 0; 1.1; 1.1; 0.15; 2.2]; % xi must be 0 or it might drift?
            % The function enforces constraints so even if start point has xi!=0 it should project or error?
            % lsqcurvefit respects bounds. Lower(3)=0, Upper(3)=0.
            % But providing a start point that violates bounds might be an issue for lsqcurvefit depending on algorithm?
            % Let's provide a valid start point.

            [params2, ~, ~] = vis.speed.fit_nospeed(sf, tf, r, 'SpecificStartPoint', startPt);
            testCase.verifyEqual(params2(3), 0, 'AbsTol', 1e-10, 'xi must be 0');
        end

        function testFitNospeedWeights(testCase)
            % A zero weight removes a corrupted point from the fit, and the
            % returned SSE is the weighted SSE.

            [SF, TF] = meshgrid([0.05 0.1 0.2 0.4], [1 2 4 8]);
            sf = SF(:); tf = TF(:);
            true_params = [10; 0; 0; 0.3; 0.4; 0.1; 2];
            r = vis.speed.tuningfunc(sf, tf, true_params);
            r_bad = r;
            r_bad(1) = r_bad(1) + 50; % corrupt one point
            w = ones(size(r));
            w(1) = 0;

            [params, sse, r2] = vis.speed.fit_nospeed(sf, tf, r_bad, 'Weights', w);
            testCase.verifyEqual(params(3), 0, 'AbsTol', 1e-10, 'xi must be 0');
            testCase.verifyLessThan(sse, 1e-3, 'Weighted SSE should ignore the zero-weight point');
            testCase.verifyGreaterThan(r2, 0.99, 'Weighted R-squared should be high');
            testCase.verifyEqual(params(6), 0.1, 'RelTol', 0.05, 'sf0 should be recovered');

            testCase.verifyError(@() vis.speed.fit_nospeed(sf, tf, r, 'Weights', w(1:end-1)), ...
                'vis:speed:fit:badWeights');
        end

    end
end
