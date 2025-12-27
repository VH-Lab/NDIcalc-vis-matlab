classdef speed_tuning < ndi.calc.tuning_fit

    methods
        function speed_tuning_obj = speed_tuning(session)
            % SPEED_TUNING - a speed_tuning demonstration of an ndi.calculator object
            %
            % SPEED_TUNING_OBJ = SPEED_TUNING(SESSION)
            %
            % Creates a SPEED_TUNING ndi.calculator object
            %
            speed_tuning_obj = speed_tuning_obj@ndi.calc.tuning_fit(session, 'speedtuning_calc', ...
                'speedtuning_calc');
            speed_tuning_obj.defaultParametersCanFunction = true;

            speed_tuning_obj.numberOfSelfTests = 18;
        end % speed_tuning()

        function doc = calculate(ndi_calculator_obj, parameters)
            % CALCULATE - perform the calculator for ndi.calc.example.speed_tuning
            %
            % DOC = CALCULATE(NDI_CALCULATOR_OBJ, PARAMETERS)
            %
            % Creates a speed_tuning_calc document given input parameters.
            %
            % The document that is created speed_tuning
            % by the input parameters.
            arguments
                ndi_calculator_obj
                parameters (1,1) struct {ndi.validators.mustHaveFields(parameters,{'input_parameters','depends_on'})}
            end

            % Step 1: set up the output structure
            speed_tuning_calc = parameters;

            tuning_response_doc = ndi_calculator_obj.session.database_search(ndi.query('base.id', 'exact_string', ...
                did.db.struct_name_value_search(parameters.depends_on, 'stimulus_tuningcurve_id'), ''));
            if numel(tuning_response_doc) ~= 1
                error('Could not find stimulus tuning doc..');
            end
            tuning_response_doc = tuning_response_doc{1};

            % Step 2: perform the calculator, which here creates a speed_tuning doc
            doc = ndi_calculator_obj.calculate_speed_indexes(tuning_response_doc) + ...
                ndi_calculator_obj.newdocument();

            if ~isempty(doc)
                doc = ndi.document(ndi_calculator_obj.doc_document_types{1}, 'speedtuning_calc', speed_tuning_calc) + doc;
                doc = doc.set_dependency_value('stimulus_tuningcurve_id', tuning_response_doc.id());
                doc = doc.set_dependency_value('element_id', tuning_response_doc.dependency_value('element_id'));
            end
        end % calculate

        function parameters = default_search_for_input_parameters(obj)
            % DEFAULT_SEARCH_FOR_INPUT_PARAMETERS - default parameters for searching for inputs
            %
            % PARAMETERS = DEFAULT_SEARCH_FOR_INPUT_PARAMETERS(OBJ)
            %
            % Returns a list of the default search parameters for finding appropriate inputs
            % to the calculator. For speed_tuning_calc, there is no appropriate default parameters
            % so this search will yield empty.
            %
            parameters.input_parameters = struct([]);
            parameters.depends_on = did.datastructures.emptystruct('name', 'value');
            parameters.query = obj.default_parameters_query(parameters);

        end % default_search_for_input_parameters

        function query = default_parameters_query(obj, parameters_specification)
            % DEFAULT_PARAMETERS_QUERY - what queries should be used to search for input parameters if none are provided?
            %
            % QUERY = DEFAULT_PARAMETERS_QUERY(OBJ, PARAMETERS_SPECIFICATION)
            %
            % When one calls SEARCH_FOR_INPUT_PARAMETERS, it is possible to specify a 'query' structure to
            % select particular documents to be placed into the parameters 'depends_on' specification.
            % If one does not provide any 'query' structure, then the default values here are used.
            %
            % The function returns:
            % |-----------------------|----------------------------------------------|
            % | query                 | A structure with 'name' and 'query' fields   |
            % |                       |   that describes a search to be performed to |
            % |                       |   identify inputs for the 'depends_on' field |
            % |                       |   in the PARAMETERS output.                  |
            % |-----------------------|-----------------------------------------------
            %
            % For the ndi.calc.stimulus.speed_tuning_calc class, this looks for
            % documents of type 'stimulus_response_scalar.json' with 'response_type' fields
            % the contain 'mean' or 'F1'.
            %
            %
            q1 = ndi.query('', 'isa', 'stimulus_tuningcurve', '');
            q2 = ndi.query('stimulus_tuningcurve.independent_variable_label', 'contains_string', 'spatial_frequency', '');
            q3 = ndi.query('stimulus_tuningcurve.independent_variable_label', 'contains_string', 'temporal_frequency', '');
            %q4 = ndi.query('stimulus_tuningcurve.independent_variable_label','hassize',[2 1],'');
            q_total = q1 & q2 & q3; % & q4;

            query = struct('name', 'stimulus_tuningcurve_id', 'query', q_total);
        end % default_parameters_query()

        function b = is_valid_dependency_input(obj, name, value)
            % IS_VALID_DEPENDENCY_INPUT - is a potential dependency input actually valid for this calculator?
            %
            % B = IS_VALID_DEPENDENCY_INPUT(OBJ, NAME, VALUE)
            %
            % Tests whether a potential input to a calculator is valid.
            % The potential dependency name is provided in NAME and its base id is
            % provided in VALUE.
            %
            % The base class behavior of this function is simply to return true, but it
            % can be overriden if additional criteria beyond an ndi.query are needed to
            % assess if a document is an appropriate input for the calculator.
            %
            b = 1;
            return;

            % could also use the below, but will require an extra query operation
            % and updating for speed

            switch lower(name)
                case lower('stimulus_tuningcurve_id')
                    q = ndi.query('base.id', 'exact_string', value, '');
                    d = obj.S.database_search(q);
                    b = (numel(d.document_properties.independent_variable_label) == 2);
                case lower('element_id')
                    b = 1;
            end
        end % is_valid_dependency_input()

        function h = plot(obj, doc_or_parameters, varargin)
            % PLOT - provide a diagnostic plot to show the results of the calculator
            %
            % H=PLOT(OBJ, DOC_OR_PARAMETERS, ...)
            %
            % Produce a plot of the tuning curve.
            %
            % Handles to the figure, the axes, and any objects created are returned in H.
            %
            % This function takes additional input arguments as name/value pairs.
            % See ndi.calculator.plot_parameters for a description of those parameters.

            % call superclass plot method to set up axes
            h = plot@ndi.calculator(obj, doc_or_parameters, varargin{:});

            if isa(doc_or_parameters, 'ndi.document')
                doc = doc_or_parameters;
            else
                error('Do not know how to proceed without an ndi document for doc_or_parameters.');
            end

            sp = doc.document_properties.speed_tuning; % shorten our typing
            tc = sp.tuning_curve; % shorten our typing
            ft = sp.fit;

            % First plot fit
            hold on;

            %h_baseline = plot([min(tc.speed) max(tc.speed)],...
            %    [0 0],'k--','linewidth',1.0001);
            %h_baseline.Annotation.LegendInformation.IconDisplayStyle = 'off';

            % now call the plot routine

            [SF, TF, MNs] = vlt.math.vector2mesh(tc.spatial_frequency, tc.temporal_frequency, tc.mean);
            MNs_fit = vis.speed.tuningfunc(SF, TF, ft.Priebe_fit_parameters);

            significant = 0;
            linestyle = '--';
            if sp.significance.visual_response_anova_p < 0.05
                significant = 1;
                linestyle = '-';
            end
            vis.speed.plottuning(SF, TF, MNs_fit, 'marker', 'none', 'linestyle', linestyle);

            % now plot raw responses
            vis.speed.plottuning(SF, TF, MNs);

            ch = get(gcf, 'children');
            currentaxes = gca;
            axes(ch(1));
            title(['Speed tuning:' num2str(ft.Priebe_fit_parameters(3)) ', pref: ' num2str(ft.Priebe_fit_parameters(7) / ft.Priebe_fit_parameters(6))]);

            element_str = '';
            ele = obj.session.database_search(ndi.query('base.id', 'exact_string', doc.dependency_value('element_id')));
            epoch_id = '';
            tc_id = obj.session.database_search(ndi.query('base.id', 'exact_string', doc.dependency_value('stimulus_tuningcurve_id')));
            if ~isempty(tc_id)
                srs = obj.session.database_search(ndi.query('base.id', 'exact_string', tc_id{1}.dependency_value('stimulus_response_scalar_id')));
                if ~isempty(srs)
                    epoch_id = [srs{1}.document_properties.stimulus_response.element_epochid ' ' srs{1}.document_properties.stimulus_response_scalar.response_type];
                end
            end

            if ~isempty(ele)
                element = ndi.database.fun.ndi_document2ndi_object(ele{1}, obj.session);
                element_str = element.elementstring();
                axes(ch(2));
                title([element_str ' ' epoch_id], 'interp', 'none');
            end
            axes(currentaxes);

            if 0 % plot function already does this
                if ~h.params.suppress_x_label
                    h.xlabel = xlabel('Speed (deg/sec)');
                end
                if ~h.params.suppress_y_label
                    h.ylabel = ylabel(['Response (' sp.properties.response_type ', ' sp.properties.response_units ')']);
                end
                box off;
            end

        end % plot()

        function speed_props_doc = calculate_speed_indexes(obj, tuning_doc)
            % CALCULATE_SPEED_INDEXES - calculate speed index values from a tuning curve
            %
            % SPEED_PROPS_DOC = CALCULATE_SPEED_INDEXES(OBJ, TUNING_DOC)
            %
            % Given a 2-dimensional tuning curve document with measurements at many spatial and
            % and temporal frequencies, this function calculates speed response
            % parameters and stores them in SPEED_TUNING document SPEED_PROPS_DOC.
            %
            %
            properties.response_units = tuning_doc.document_properties.stimulus_tuningcurve.response_units;

            stim_response_doc = obj.session.database_search(ndi.query('base.id', ...
                'exact_string', tuning_doc.dependency_value('stimulus_response_scalar_id'), ''));
            if numel(stim_response_doc) ~= 1
                error('Could not find stimulus response scalar document.');
            end
            if iscell(stim_response_doc)
                stim_response_doc = stim_response_doc{1};
            end

            properties.response_type = stim_response_doc.document_properties.stimulus_response_scalar.response_type;

            sf_coord = 1;
            tf_coord = 2;
            if contains(tuning_doc.document_properties.stimulus_tuningcurve.independent_variable_label{1}, 'temporal', 'IgnoreCase', true)
                sf_coord = 2;
                tf_coord = 1;
            end

            resp = ndi.app.stimulus.tuning_response.tuningcurvedoc2vhlabrespstruct(tuning_doc);

            [anova_across_stims, anova_across_stims_blank] = neural_response_significance(resp);

            tuning_curve = struct(...
                'spatial_frequency', vlt.data.colvec(tuning_doc.document_properties.stimulus_tuningcurve.independent_variable_value(:, 1)), ...
                'temporal_frequency', vlt.data.colvec(tuning_doc.document_properties.stimulus_tuningcurve.independent_variable_value(:, 2)), ...
                'mean', vlt.data.colvec(resp.curve(2, :)), ...
                'stddev', vlt.data.colvec(resp.curve(3, :)), ...
                'stderr', vlt.data.colvec(resp.curve(4, :)), ...
                'individual', vlt.data.cellarray2mat(resp.ind), ...
                'control_stddev', resp.blankresp(2), ...
                'control_stderr', resp.blankresp(3));

            significance = struct('visual_response_anova_p', anova_across_stims_blank, ...
                'across_stimuli_anova_p', anova_across_stims);

            [f, sse_withspeed, r2_withspeed] = vis.speed.fit(tuning_curve.spatial_frequency(:), tuning_curve.temporal_frequency(:), tuning_curve.mean(:));
            sfs = logspace(0.01, 60, 200);
            tfs = logspace(0.01, 120, 200);
            [SFs, TFs] = meshgrid(sfs, tfs);
            fit_values = vis.speed.tuningfunc(SFs(:), TFs(:), f);

            fit.Priebe_fit_parameters = f;
            fit.Priebe_fit_spatial_frequencies = SFs(:);
            fit.Priebe_fit_temporal_frequencies = TFs(:);
            fit.Priebe_fit_values = fit_values;
            fit.Priebe_fit_speed_tuning_index = fit.Priebe_fit_parameters(3);
            fit.Priebe_fit_spatial_frequency_preference = fit.Priebe_fit_parameters(6);
            fit.Priebe_fit_temporal_frequency_preference = fit.Priebe_fit_parameters(7);
            fit.sse = sse_withspeed;
            fit.r_squared = r2_withspeed;

            %add fit with speed parameter set to 0
            [f_no_speed, sse_nospeed, r2_nospeed] = vis.speed.fit_nospeed(tuning_curve.spatial_frequency(:), tuning_curve.temporal_frequency(:), tuning_curve.mean(:));
            fit_no_speed_values = vis.speed.tuningfunc(SFs(:), TFs(:), f_no_speed);
            fit_no_speed.Priebe_fit_parameters = f_no_speed;
            fit_no_speed.Priebe_fit_spatial_frequencies = SFs(:);
            fit_no_speed.Priebe_fit_temporal_frequencies = TFs(:);
            fit_no_speed.Priebe_fit_values = fit_no_speed_values;
            fit_no_speed.Priebe_fit_speed_tuning_index = fit_no_speed.Priebe_fit_parameters(3);
            fit_no_speed.Priebe_fit_spatial_frequency_preference = fit_no_speed.Priebe_fit_parameters(6);
            fit_no_speed.Priebe_fit_temporal_frequency_preference = fit_no_speed.Priebe_fit_parameters(7);
            fit_no_speed.sse = sse_nospeed;
            fit_no_speed.r_squared = r2_nospeed;
            
            %add nested-F test
            num_responses = numel(fit_values);
            fit_no_speed.Priebe_fit_nested_F_test_p_value = vis.speed.speed_nested_f(num_responses, sse_withspeed, sse_nospeed);

            speed_tuning.properties = properties;
            speed_tuning.tuning_curve = tuning_curve;
            speed_tuning.significance = significance;
            speed_tuning.fit = fit;
            speed_tuning.fit_no_speed = fit_no_speed;

            speed_props_doc = ndi.document('speed_tuning', ...
                'speed_tuning', speed_tuning);
            speed_props_doc = speed_props_doc.set_dependency_value('element_id', ...
                tuning_doc.dependency_value('element_id'));
            speed_props_doc = speed_props_doc.set_dependency_value('stimulus_tuningcurve_id', tuning_doc.id());
        end % calculate_speed_indexes()

        % TESTING METHODS

        function [param_struct, independent_variable, x, r] = generate_mock_parameters(obj, scope, index)
            % generate_mock_parameters - generate mock parameters for testing ndi.calc.vis.oridir_tuning
            %
            % [PARAM_STRUCT, INDEPENDENT_VARIABLE, X, R] = GENERATE_MOCK_PARAMETERS(OBJ, SCOPE, INDEX)
            %
            % Generates a parameter set for generating a mock document with a given index value.
            % P will be a row vector of parameters [a1 b1 a2 b2].
            %
            % SCOPE can be 'standard', 'random_nonoise', or 'random_noisy'.
            % INDEX selects which parameters are used to generate a mock document (from 1..TOTAL, wrapped
            % using MOD).
            %
            %1st test taken from demo script
            %2nd test neuron not selective for speed
            %3rd test neuron maximally selective for speed
            %4th through 7th: tiling TF/SF space with TF-tuned
            %cells
            %8th through 11th: tiling TF/SF space with SF-tuned
            %cells
            %12th through 15th: tiling TF/SF space with speed-tuned
            %cells
            %16th: try a different peak neuron response?
            %17th and 18th: skew temporal freq tuning curve in
            %different directions?
            %add more tests by adding to end of each parameter's vector
            A = [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 1, 2, 2]; %Peak response of the neuron
            zeta = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, .5, -.5]; %Skew of the temporal freq tuning curve
            xi = [.2, 0.0001, 1, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 1, 1, 1, 1, .5, .5, .5]; %speed tuning index (0 to 1)
            sigma_sf = [1, 1, 1, 1.2, 1.2, 1.2, 1.2, .5, .5, .5, .5, 1, 1, 1, 1, 1, 1, 1]; % Tuning width of the neuron for spatial frequency (Cycles per degree)
            sigma_tf = [1, 1, 1, .5, .5, .5, .5, 1.2, 1.2, 1.2, 1.2, 1, 1, 1, 1, 1, 1, 1]; % Tuning width of the neuron for temporal frequency (Cycles per second)
            sf0 = [0.2, 0.2, 0.2, .2, .8, .2, .8, .2, .8, .2, .8, .2, .8, .2, .8, sqrt(2) / 5, sqrt(2) / 5, sqrt(2) / 5]; %preferred spatial frequency averaged across temporal frequencies
            tf0 = [2, 2, 2, 2, 2, 16, 16, 2, 2, 16, 16, 2, 2, 16, 16, 4, 4, 4]; %preferred temporal frequency averaged across spatial frequencies

            %for trying out new tests:
            % A = 20;
            % zeta = 0;
            % xi = .000009;
            % sigma_sf = 1;
            % sigma_tf = 1;
            % sf0 = sqrt(2)/5;
            % tf0 = 4;
            P_ = [A(:) zeta(:) xi(:) sigma_sf(:) sigma_tf(:) sf0(:) tf0(:)];

            % %1st test from spatial_frequency_tuning:
            % sf_params = [1 1 0 1];
            % %1st test from temporal_frequency_tuning:
            % tf_params = [1 1 0 1];
            % P_ = [sf_params tf_params];
            total = size(P_, 1);

            actual_index = 1 + mod(index - 1, total);

            % no dependence on scope for this stimulus type

            P = P_(actual_index, :);

            %taken from calculate_speed_indexes method:
            numsteps = 5;
            %sfs = logspace(log10(0.01),log10(60),numsteps);
            sfs = [0.05 0.08 0.1 0.2 0.4 0.8 1.2]; %taken from demo
            %tfs = logspace(log10(0.01),log10(120),numsteps);
            tfs = [0.5 1 2 4 8 16 32]; %taken from demo
            [SFs, TFs] = meshgrid(sfs, tfs);
            function_params = P;
            r_ = vlt.neuro.vision.speed.tuningfunc(SFs, TFs, function_params);

            %r = vlt.math.dog(speed_values,function_params);
            %param_struct = struct([]); %this didn't work - need at
            %least one parameter?
            param_struct = struct('contrast', .5); %random choice, can be anything between 0 and 1
            independent_variable = {'temporal_frequency', 'spatial_frequency'};
            x = [SFs(:), TFs(:)]; % columns
            r = r_(:); % columns
            %why do we have these? Set control (blank) stimulus to
            %firing rate = 0?
            %should we have just one blank stimulus row or
            %multiple? Add nan to end of x or end of sfs and tfs?
            x(end + 1, :) = NaN;
            r(end + 1, 1) = 0;

        end % generate_mock_parameters
    end % methods()
end % speed_tuning

