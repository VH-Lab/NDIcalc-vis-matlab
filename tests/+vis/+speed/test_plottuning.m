classdef test_plottuning < matlab.unittest.TestCase

    methods (Test)

        function testPlotTuningRuns(testCase)
            % Test that vis.speed.plottuning runs without error

            % Mock data
            [SF, TF] = meshgrid([0.1 0.2], [1 2]);
            R = rand(size(SF));

            % Create a figure to plot into (invisible)
            f = figure('Visible', 'off');
            cleanup = onCleanup(@() close(f));

            try
                vis.speed.plottuning(SF, TF, R, 'do_surf', 0);
                % If it runs without error, good.

                vis.speed.plottuning(SF, TF, R, 'do_surf', 1);
                 % If it runs without error, good.
            catch e
                testCase.verifyFail(['plottuning failed with error: ' e.message]);
            end

        end
    end
end
