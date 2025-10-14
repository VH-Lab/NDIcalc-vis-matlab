function [fh, nodes] = trianglePipeline(options)
% TRIANGLEPIPELINE - Create a figure demonstrating a processing pipeline
%
% [FH, NODES] = NDI.CALC.PAPERFIGURES.TRIANGLEPIPELINE(...)
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
% Outputs:
% | Name | Description |
% |---|---|
% | FH | The figure handle. |
% | NODES | A cell array of the ndi.calc.paperFigures.triangleNode objects. |
%
    arguments
        options.columnSpacing (1,1) {mustBeNumeric} = 30;
        options.rowSpacing (1,1) {mustBeNumeric} = 10;
        options.nodeWidth (1,1) {mustBeNumeric} = 3;
    end

    fh = figure();
    hold on;

    nodes = {};

    % Define positions
    c1 = 0;
    c2 = c1 + options.columnSpacing;
    c3 = c2 + options.columnSpacing;
    c4 = c3 + options.columnSpacing;
    c5 = c4 + options.columnSpacing;
    c6 = c5 + options.columnSpacing;

    r1 = 0;
    r2 = r1 - options.rowSpacing;
    r3 = r2 - options.rowSpacing;
    r4 = r3 - options.rowSpacing;
    r5 = r4 - options.rowSpacing;
    r6 = r5 - options.rowSpacing;
    r7 = r6 - options.rowSpacing;

    node_w = options.nodeWidth;

    % Row 1
    spikes = ndi.calc.paperFigures.triangleNode(...
        'shape', 'rectangle', ...
        'numberOfInputs', 0, ...
        'name', 'spikes', ...
        'position', [c1 r1], ...
        'width', node_w*3, ...
        'height', node_w*2);
    nodes{end+1} = spikes;

    node2 = ndi.calc.paperFigures.triangleNode(...
        'shape', 'triangle', ...
        'numberOfInputs', 2, ...
        'name', {'Select','neurons'}, ...
        'position', [c2 r1], ...
        'width', node_w*3, ...
        'height', node_w*3);
    nodes{end+1} = node2;

    node3 = ndi.calc.paperFigures.triangleNode(...
        'shape', 'triangle', ...
        'numberOfInputs', 2, ...
        'name', {'Select','epoch'}, ...
        'position', [c3 r1], ...
        'width', node_w*3, ...
        'height', node_w*3);
    nodes{end+1} = node3;

    % Row 2
    index_node = ndi.calc.paperFigures.triangleNode(...
        'shape', 'rectangle', ...
        'name', 'index_12345', ...
        'position', [c1 r2], ...
        'width', node_w*3, ...
        'height', node_w*2);
    nodes{end+1} = index_node;

    epoch_node = ndi.calc.paperFigures.triangleNode(...
        'shape', 'rectangle', ...
        'name', 'epoch_xyz', ...
        'position', [c2 r2], ...
        'width', node_w*3, ...
        'height', node_w*2);
    nodes{end+1} = epoch_node;

    % Row 3
    stim_node = ndi.calc.paperFigures.triangleNode(...
        'shape', 'rectangle', ...
        'name', 'stimulator', ...
        'position', [c1 r3], ...
        'width', node_w*3, ...
        'height', node_w*2);
    nodes{end+1} = stim_node;

    node7 = ndi.calc.paperFigures.triangleNode(...
        'shape', 'triangle', ...
        'numberOfInputs', 2, ...
        'name', {'Select','epoch'}, ...
        'position', [c2 r3], ...
        'width', node_w*3, ...
        'height', node_w*3);
    nodes{end+1} = node7;

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

    % Add the "Sort stimuli" nodes
    rows = [r4, r5, r6, r7];
    for r = 1:numel(rows)
        sort_node = ndi.calc.paperFigures.triangleNode(...
            'shape', 'triangle', ...
            'numberOfInputs', 1, ...
            'name', {'Sort','stimuli'}, ...
            'position', [c3 rows(r)], ...
            'width', node_w*3, ...
            'height', node_w*3);
        nodes{end+1} = sort_node;
        ndi.calc.paperFigures.triangleNode.plot_hvh_line(node7.outputPort, sort_node.inputPorts(1,:));
    end

    hold off;
    axis equal;
    axis off;

end % trianglePipeline()