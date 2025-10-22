classdef triangleNode < handle
    % TRIANGLENODE - A graphical object for creating triangle/rectangle nodes in a figure.
    %
    % This class creates a graphical node, shaped as a triangle or a rectangle,
    % commonly used in diagrams like schematics. It allows for customization of
    % its appearance and behavior, such as size, number of inputs, and title.
    % The object is a handle, so any changes to its properties will automatically
    % update the drawing in the figure.

    properties
        width = 3;
        height = 3;
        numberOfInputs = 2;
        name = '';
        titleLocation = 'middle';
        titleFontName = 'Helvetica';
        titleFontSize = 0.5;
        show = true;
        shape = 'triangle';
        position = [0 0]; % Center position [x, y]
        inputTerminalColor = [1 0 0]; % red
        outputTerminalColor = [0 1 0]; % green
        fillColor = 'none';

        % Graphics handles
        shapeHandle
        inputLines
        outputLine
        titleHandle

        % Input/Output port locations
        inputPorts
        outputPort
    end

    properties (Access = private)
        isConstructing = false;
    end

    methods
        function obj = triangleNode(options)
            % TRIANGLENODE - Constructor for the triangleNode class.
            %
            %   OBJ = NDI.CALC.PAPERFIGURES.TRIANGLENODE(Name, Value, ...)
            %
            %   Creates a triangleNode object with specified properties.
            %
            %   Optional Name-Value pair arguments:
            %   'width' - The width of the node's base. Default is 3.
            %   'numberOfInputs' - The number of input lines. Default is 2.
            %   'name' - A title for the node. Default is empty.
            %   'titleLocation' - 'middle', 'above', or 'below'. Default is 'middle'.
            %   'titleFontName' - Font name for the title. Default is 'Helvetica'.
            %   'titleFontSize' - Font size for the title. Default is 12.
            %   'show' - Whether to plot the node. Default is true.
            %   'shape' - 'triangle' or 'rectangle'. Default is 'triangle'.
            %   'position' - The [x,y] center position of the node. Default is [0 0].
            %
            arguments
                options.width (1,1) {mustBeNumeric} = 3;
                options.height (1,1) {mustBeNumeric} = 3;
                options.numberOfInputs (1,1) {mustBeNumeric} = 2;
                options.name {mustBeA(options.name,{'char','cell'})} = '';
                options.titleLocation (1,:) char {mustBeMember(options.titleLocation,{'middle','above','below'})} = 'middle';
                options.titleFontName (1,:) char = 'Helvetica';
                options.titleFontSize (1,1) {mustBeNumeric} = 0.5;
                options.show (1,1) {mustBeNumericOrLogical} = true;
                options.shape (1,:) char {mustBeMember(options.shape,{'triangle','rectangle'})} = 'triangle';
                options.position (1,2) {mustBeNumeric} = [0 0];
                options.inputTerminalColor = [1 0 0];
                options.outputTerminalColor = [0 1 0];
                options.fillColor = 'none';
            end

            obj.isConstructing = true;

            % Assign properties from arguments block
            obj.width = options.width;
            obj.height = options.height;
            obj.numberOfInputs = options.numberOfInputs;
            obj.name = options.name;
            obj.titleLocation = options.titleLocation;
            obj.titleFontName = options.titleFontName;
            obj.titleFontSize = options.titleFontSize;
            obj.show = options.show;
            obj.shape = options.shape;
            obj.position = options.position;
            obj.inputTerminalColor = options.inputTerminalColor;
            obj.outputTerminalColor = options.outputTerminalColor;

            if strcmp(options.fillColor,'none')
                if strcmpi(obj.shape,'triangle')
                    obj.fillColor = [1 1 0.9];
                else
                    obj.fillColor = [0.9 0.9 1];
                end
            else
                obj.fillColor = options.fillColor;
            end

            obj.isConstructing = false;

            % Explicitly plot the node
            obj.plotNode();
        end

        function set.width(obj, val)
            obj.width = val;
            if ~obj.isConstructing, obj.plotNode(); end
        end

        function set.height(obj, val)
            obj.height = val;
            if ~obj.isConstructing, obj.plotNode(); end
        end

        function set.numberOfInputs(obj, val)
            obj.numberOfInputs = val;
            if ~obj.isConstructing, obj.plotNode(); end
        end

        function set.name(obj, val)
            obj.name = val;
            if ~obj.isConstructing, obj.plotNode(); end
        end

        function set.titleLocation(obj, val)
            obj.titleLocation = val;
            if ~obj.isConstructing, obj.plotNode(); end
        end

        function set.titleFontName(obj, val)
            obj.titleFontName = val;
            if ~obj.isConstructing, obj.plotNode(); end
        end

        function set.titleFontSize(obj, val)
            obj.titleFontSize = val;
            if ~obj.isConstructing, obj.plotNode(); end
        end

        function set.show(obj, val)
            obj.show = val;
            if ~obj.isConstructing, obj.plotNode(); end
        end

        function set.shape(obj, val)
            obj.shape = val;
            if ~obj.isConstructing, obj.plotNode(); end
        end

        function set.position(obj, val)
            obj.position = val;
            if ~obj.isConstructing, obj.plotNode(); end
        end

        function set.inputTerminalColor(obj, val)
            obj.inputTerminalColor = val;
            if ~obj.isConstructing, obj.plotNode(); end
        end

        function set.outputTerminalColor(obj, val)
            obj.outputTerminalColor = val;
            if ~obj.isConstructing, obj.plotNode(); end
        end

        function set.fillColor(obj, val)
            obj.fillColor = val;
            if ~obj.isConstructing, obj.plotNode(); end
        end

        function delete(obj)
            % DELETE - cleans up graphics objects
            % Per user request, we are no longer deleting the graphics handles
            % when the object is deleted.
            %
            % if isgraphics(obj.shapeHandle), delete(obj.shapeHandle); end
            % if isgraphics(obj.titleHandle), delete(obj.titleHandle); end
            % delete(obj.inputLines(isgraphics(obj.inputLines)));
            % if isgraphics(obj.outputLine), delete(obj.outputLine); end
        end

    end

    methods (Static)
        function plot_hvh_line(start_point, end_point, color)
            % PLOT_HVH_LINE - plots a horizontal-vertical-horizontal line

            x = [start_point(1), ...
                 start_point(1) + (end_point(1)-start_point(1))/2, ...
                 start_point(1) + (end_point(1)-start_point(1))/2, ...
                 end_point(1)];

            y = [start_point(2), ...
                 start_point(2), ...
                 end_point(2), ...
                 end_point(2)];

            plot(x,y,'Color',color);
        end
    end

    methods (Access = private)
        function plotNode(obj)
            % PLOTNODE - Draws or updates the node's graphical representation.
            disp(['Plotting node: ' obj.name]);

            % Clean up old graphics
            if isgraphics(obj.shapeHandle), delete(obj.shapeHandle); end
            if isgraphics(obj.titleHandle), delete(obj.titleHandle); end
            if isgraphics(obj.inputLines), delete(obj.inputLines); end
            if isgraphics(obj.outputLine), delete(obj.outputLine); end

            if ~obj.show, return; end % Do nothing if not shown

            % Node geometry
            w = obj.width;
            h = obj.height;
            pos = obj.position;

            if strcmpi(obj.shape, 'triangle')
                % Triangle vertices centered at obj.position
                vertices = [pos(1)-w/2, pos(2)-h/2; pos(1)-w/2, pos(2)+h/2; pos(1)+w/2, pos(2)];
                obj.shapeHandle = patch('Vertices', vertices, 'Faces', [1 2 3], 'FaceColor', obj.fillColor, 'EdgeColor', 'k', 'LineWidth', 2);
            else % Rectangle
                x = pos(1) - w/2;
                y = pos(2) - h/2;
                vertices = [x, y; x+w, y; x+w, y+h; x, y+h];
                obj.shapeHandle = patch('Vertices', vertices, 'Faces', [1 2 3 4], 'FaceColor', obj.fillColor, 'EdgeColor', 'k', 'LineWidth', 2);
            end

            % Input and Output ports
            % Output port
            obj.outputPort = [pos(1) + w/2, pos(2)];

            % Input ports
            obj.inputPorts = zeros(obj.numberOfInputs, 2);
            if obj.numberOfInputs > 0
                y_inputs = linspace(pos(2)-h/2, pos(2)+h/2, obj.numberOfInputs+2);
                y_inputs = y_inputs(2:end-1); % Remove top and bottom points
                for i = 1:obj.numberOfInputs
                    obj.inputPorts(i,:) = [pos(1)-w/2, y_inputs(i)];
                end
            end

            % Draw input/output lines
            obj.inputLines = gobjects(obj.numberOfInputs, 1);
            for i=1:obj.numberOfInputs
                p_start = obj.inputPorts(i,:);
                p_end = [p_start(1)-w/4, p_start(2)];
                obj.inputLines(i) = plot([p_start(1) p_end(1)], [p_start(2) p_end(2)], 'Color', obj.inputTerminalColor, 'LineWidth', 1);
            end

            p_start = obj.outputPort;
            p_end = [p_start(1)+w/4, p_start(2)];
            obj.outputLine = plot([p_start(1) p_end(1)], [p_start(2) p_end(2)], 'Color', obj.outputTerminalColor, 'LineWidth', 1);

            % Title
            title_pos = [0, 0];
            switch obj.titleLocation
                case 'middle'
                    title_pos = [pos(1), pos(2)];
                    if strcmpi(obj.shape, 'triangle')
                        title_pos(1) = pos(1) - w/6;
                    end
                case 'above'
                    title_pos = [pos(1), pos(2) + h/2 + h*0.2];
                case 'below'
                    title_pos = [pos(1), pos(2) - h/2 - h*0.2];
            end

            obj.titleHandle = text(title_pos(1), title_pos(2), obj.name);
            set(obj.titleHandle, 'FontUnits', 'centimeters');
            set(obj.titleHandle, ...
                'FontName', obj.titleFontName, ...
                'FontSize', obj.titleFontSize, ...
                'HorizontalAlignment', 'center', ...
                'Interpreter', 'none');

        end
    end
end