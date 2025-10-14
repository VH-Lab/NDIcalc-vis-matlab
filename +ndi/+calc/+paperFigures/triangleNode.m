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
        numberOfInputs = 2;
        name = '';
        titleLocation = 'middle';
        titleFontName = 'Helvetica';
        titleFontSize = 12;
        show = true;
        shape = 'triangle';
        position = [0 0]; % Center position [x, y]

        % Graphics handles
        shapeHandle
        inputLines
        outputLine
        titleHandle

        % Input/Output port locations
        inputPorts
        outputPort
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
                options.numberOfInputs (1,1) {mustBeNumeric} = 2;
                options.name (1,:) char = '';
                options.titleLocation (1,:) char {mustBeMember(options.titleLocation,{'middle','above','below'})} = 'middle';
                options.titleFontName (1,:) char = 'Helvetica';
                options.titleFontSize (1,1) {mustBeNumeric} = 12;
                options.show (1,1) {mustBeNumericOrLogical} = true;
                options.shape (1,:) char {mustBeMember(options.shape,{'triangle','rectangle'})} = 'triangle';
                options.position (1,2) {mustBeNumeric} = [0 0];
            end

            % Assign properties from arguments
            obj.width = options.width;
            obj.numberOfInputs = options.numberOfInputs;
            obj.name = options.name;
            obj.titleLocation = options.titleLocation;
            obj.titleFontName = options.titleFontName;
            obj.titleFontSize = options.titleFontSize;
            obj.shape = options.shape;
            obj.position = options.position;

            % Defer plotting until 'show' is set
            obj.show = options.show;

            % Explicitly plot the node
            obj.plotNode();
        end

        function set.width(obj, val)
            obj.width = val;
            obj.plotNode();
        end

        function set.numberOfInputs(obj, val)
            obj.numberOfInputs = val;
            obj.plotNode();
        end

        function set.name(obj, val)
            obj.name = val;
            obj.plotNode();
        end

        function set.titleLocation(obj, val)
            obj.titleLocation = val;
            obj.plotNode();
        end

        function set.titleFontName(obj, val)
            obj.titleFontName = val;
            obj.plotNode();
        end

        function set.titleFontSize(obj, val)
            obj.titleFontSize = val;
            obj.plotNode();
        end

        function set.show(obj, val)
            obj.show = val;
            obj.plotNode();
        end

        function set.shape(obj, val)
            obj.shape = val;
            obj.plotNode();
        end

        function set.position(obj, val)
            obj.position = val;
            obj.plotNode();
        end

        function delete(obj)
            % DELETE - cleans up graphics objects
            if isgraphics(obj.shapeHandle), delete(obj.shapeHandle); end
            if isgraphics(obj.titleHandle), delete(obj.titleHandle); end
            delete(obj.inputLines(isgraphics(obj.inputLines)));
            if isgraphics(obj.outputLine), delete(obj.outputLine); end
        end

    end

    methods (Static)
        function plot_hvh_line(start_point, end_point)
            % PLOT_HVH_LINE - plots a horizontal-vertical-horizontal line

            x = [start_point(1), ...
                 start_point(1) + (end_point(1)-start_point(1))/2, ...
                 start_point(1) + (end_point(1)-start_point(1))/2, ...
                 end_point(1)];

            y = [start_point(2), ...
                 start_point(2), ...
                 end_point(2), ...
                 end_point(2)];

            plot(x,y,'k-');
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
            h = w; % Height is same as width for a nice look
            pos = obj.position;

            if strcmpi(obj.shape, 'triangle')
                % Triangle vertices centered at obj.position
                vertices = [pos(1)-w/2, pos(2)-h/2; pos(1)-w/2, pos(2)+h/2; pos(1)+w/2, pos(2)];
                % Close the triangle
                plot_verts = [vertices; vertices(1,:)];
                obj.shapeHandle = plot(plot_verts(:,1), plot_verts(:,2), 'k', 'LineWidth', 2);
            else % Rectangle
                x = pos(1) - w/2;
                y = pos(2) - h/2;
                obj.shapeHandle = rectangle('Position', [x, y, w, h], 'LineWidth', 2);
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
                obj.inputLines(i) = plot([p_start(1) p_end(1)], [p_start(2) p_end(2)], 'k', 'LineWidth', 1);
            end

            p_start = obj.outputPort;
            p_end = [p_start(1)+w/4, p_start(2)];
            obj.outputLine = plot([p_start(1) p_end(1)], [p_start(2) p_end(2)], 'k', 'LineWidth', 1);

            % Title
            title_pos = [0, 0];
            switch obj.titleLocation
                case 'middle'
                    title_pos = [pos(1), pos(2)];
                case 'above'
                    title_pos = [pos(1), pos(2) + h/2 + h*0.2];
                case 'below'
                    title_pos = [pos(1), pos(2) - h/2 - h*0.2];
            end

            obj.titleHandle = text(title_pos(1), title_pos(2), obj.name, ...
                'FontName', obj.titleFontName, ...
                'FontSize', obj.titleFontSize, ...
                'HorizontalAlignment', 'center');

        end
    end
end