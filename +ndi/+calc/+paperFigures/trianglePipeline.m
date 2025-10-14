function fh = trianglePipeline(options)
% TRIANGLEPIPELINE - Create a figure demonstrating a processing pipeline
%
% FH = NDI.CALC.PAPERFIGURES.TRIANGLEPIPELINE(...)
%
% Creates a figure that illustrates a data processing pipeline using
% ndi.calc.paperFigures.triangleNode objects.
%
% This function can be called with name-value pairs:
%
% | Name-Value Pair | Description | Default |
% |---|---|---|
% | 'columnSpacing' | Horizontal distance between node centers | 10 |
% | 'rowSpacing' | Vertical distance between node centers | 5 |
% | 'nodeWidth' | The width of the nodes | 3 |
%
    arguments
        options.columnSpacing (1,1) {mustBeNumeric} = 10;
        options.rowSpacing (1,1) {mustBeNumeric} = 5;
        options.nodeWidth (1,1) {mustBeNumeric} = 3;
    end

    fh = figure();
    hold on;

    % Define positions
    c1 = 0;
    c2 = c1 + options.columnSpacing;
    c3 = c2 + options.columnSpacing;

    r1 = 0;
    r2 = r1 - options.rowSpacing;
    r3 = r2 - options.rowSpacing;

    node_w = options.nodeWidth;

    % Row 1
    spikes = ndi.calc.paperFigures.triangleNode(...
        'shape', 'rectangle', ...
        'numberOfInputs', 0, ...
        'name', 'spikes', ...
        'position', [c1 r1], ...
        'width', node_w);

    node2 = ndi.calc.paperFigures.triangleNode(...
        'shape', 'triangle', ...
        'numberOfInputs', 2, ...
        'name', 'Select neurons', ...
        'position', [c2 r1], ...
        'width', node_w);

    node3 = ndi.calc.paperFigures.triangleNode(...
        'shape', 'triangle', ...
        'numberOfInputs', 2, ...
        'name', 'Select epoch', ...
        'position', [c3 r1], ...
        'width', node_w);

    % Row 2
    index_node = ndi.calc.paperFigures.triangleNode(...
        'shape', 'rectangle', ...
        'name', 'index_12345', ...
        'position', [c1 r2], ...
        'width', node_w);

    epoch_node = ndi.calc.paperFigures.triangleNode(...
        'shape', 'rectangle', ...
        'name', 'epoch_xyz', ...
        'position', [c2 r2], ...
        'width', node_w);

    % Row 3
    stim_node = ndi.calc.paperFigures.triangleNode(...
        'shape', 'rectangle', ...
        'name', 'stimulator', ...
        'position', [c1 r3], ...
        'width', node_w);

    node7 = ndi.calc.paperFigures.triangleNode(...
        'shape', 'triangle', ...
        'numberOfInputs', 2, ...
        'name', 'Select epoch', ...
        'position', [c2 r3], ...
        'width', node_w);

    % Draw connections

    % Spikes -> Node 2 (upper)
    ndi.calc.paperFigures.triangleNode.plot_hvh_line(spikes.outputPort, node2.inputPorts(1,:));

    % Node 2 -> Node 3 (upper)
    ndi.calc.paperFigures.triangleNode.plot_hvh_line(node2.outputPort, node3.inputPorts(1,:));

    % index_12345 -> Node 2 (lower)
    ndi.calc.paperFigures.triangleNode.plot_hvh_line(index_node.outputPort, node2.inputPorts(2,:));

    % epoch_xyz -> Node 3 (lower)
    ndi.calc.paperFigures.triangleNode.plot_hvh_line(epoch_node.outputPort, node3.inputPorts(2,:));

    % epoch_xyz -> Node 7 (lower)
    ndi.calc.paperFigures.triangleNode.plot_hvh_line(epoch_node.outputPort, node7.inputPorts(2,:));

    % stimulator -> Node 7 (upper)
    ndi.calc.paperFigures.triangleNode.plot_hvh_line(stim_node.outputPort, node7.inputPorts(1,:));

    hold off;
    axis equal;
    axis off;

end % trianglePipeline()