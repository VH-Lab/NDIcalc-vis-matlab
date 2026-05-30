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

        function testGrayscaleColormapAndContour(testCase)
            % Test that plottuning draws a grayscale heat map and overlays
            % the half-height contour when given model parameters.

            [SF, TF] = meshgrid([0.05 0.1 0.2 0.4], [1 2 4 8]);
            P = [10, 0, 0, 0.5, 0.5, 0.1, 4];
            R = vis.speed.tuningfunc(SF, TF, P);

            f = figure('Visible', 'off');
            cleanup = onCleanup(@() close(f));

            vis.speed.plottuning(SF, TF, R, 'half_contour_params', P);

            % Locate the axes that contains the surface heat map.
            ax = findobj(f, 'Type', 'axes');
            surfax = [];
            for k = 1:numel(ax)
                if ~isempty(findobj(ax(k), 'Type', 'surface'))
                    surfax = ax(k);
                    break;
                end
            end
            testCase.verifyNotEmpty(surfax, 'Expected an axes containing a surface heat map.');

            % The colormap should be grayscale (black to white).
            testCase.verifyEqual(get(surfax, 'Colormap'), gray(256), 'AbsTol', 1e-12, ...
                'Heat map colormap should ramp from black to white.');

            % The half-height contour should be drawn as a line on the heat map.
            contour_lines = findobj(surfax, 'Type', 'line');
            testCase.verifyNotEmpty(contour_lines, ...
                'Expected the half-height contour to be overlaid on the heat map.');
        end

        function testCustomColormap(testCase)
            % Test that a custom colormap is honored.

            [SF, TF] = meshgrid([0.1 0.2 0.4], [1 2 4]);
            R = rand(size(SF));

            f = figure('Visible', 'off');
            cleanup = onCleanup(@() close(f));

            customMap = parula(64);
            vis.speed.plottuning(SF, TF, R, 'surf_colormap', customMap);

            ax = findobj(f, 'Type', 'axes');
            surfax = [];
            for k = 1:numel(ax)
                if ~isempty(findobj(ax(k), 'Type', 'surface'))
                    surfax = ax(k);
                    break;
                end
            end
            testCase.verifyNotEmpty(surfax, 'Expected an axes containing a surface heat map.');
            testCase.verifyEqual(get(surfax, 'Colormap'), customMap, 'AbsTol', 1e-12, ...
                'Custom colormap should be applied to the heat map.');
        end
    end
end
