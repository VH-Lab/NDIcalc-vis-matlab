classdef speed_tuning_bootstrap < ndi.calc.tuning_fit
    % SPEED_TUNING_BOOTSTRAP - bootstrap version of the speed_tuning calculator
    %
    % =====================================================================
    % SCAFFOLD / DRAFT -- NOT YET IMPLEMENTED OR RUN.
    % This class is a starting point for discussion (see the design issue).
    % It has not been executed (no licensed MATLAB in the authoring env; per
    % AGENTS.md, changes here are verified by CI on push). Sections marked TODO
    % are intentionally left for the implementation session after the design is
    % agreed. Do not assume any method below is complete.
    % =====================================================================
    %
    % SPEED_TUNING_BOOTSTRAP computes the same spatiotemporal-coupling fits as
    % ndi.calc.vis.speed_tuning, but repeats the fit over S bootstrap resamples
    % of the individual trials, so that a confidence interval can be placed on
    % the coupling index xi and the other fit parameters. The motivation is the
    % Suarez-Casanova et al. review: cells that "pass" the xi==1 nested F-test
    % but have a fitted xi<0 are poorly-constrained fits, and a CI on xi
    % (rather than an accept-the-null test) separates genuinely speed-tuned
    % cells from unconstrained ones.
    %
    % Output document type: 'speedtuning_bootstrap_calc'
    % Result document type:  'speed_tuning_bootstrap'
    %
    % The per-sample fit parameters are stored as an S-by-P matrix (S bootstrap
    % samples in rows, P Priebe parameters in columns), so that
    % prctile(params,[2.5 97.5]) returns a 2-by-P per-parameter CI with no
    % transpose. (Orientation is a design decision -- see the issue.)
    %
    % See also: ndi.calc.vis.speed_tuning, vis.speed.fit, vis.speed.fit_nospeed,
    %           vis.speed.fit_fullspeed, vis.speed.speed_nested_f

    methods
        function obj = speed_tuning_bootstrap(session)
            % SPEED_TUNING_BOOTSTRAP - create a speed_tuning_bootstrap calculator
            %
            % OBJ = SPEED_TUNING_BOOTSTRAP(SESSION)
            %
            ndi.fun.checkCalcDirectory();

            obj = obj@ndi.calc.tuning_fit(session, 'speedtuning_bootstrap_calc', ...
                'speedtuning_bootstrap_calc');
            obj.defaultParametersCanFunction = true;

            % TODO: set once self-test mocks exist (see +vis/mock/speed_tuning_bootstrap/).
            obj.numberOfSelfTests = 0;
        end % speed_tuning_bootstrap()

        function doc = calculate(ndi_calculator_obj, parameters)
            % CALCULATE - run the bootstrap speed-tuning calculator
            %
            % DOC = CALCULATE(NDI_CALCULATOR_OBJ, PARAMETERS)
            %
            % Creates a speedtuning_bootstrap_calc document. Input parameters
            % mirror speed_tuning (min_xi, max_xi) plus:
            %   input_parameters.numBootstrap : number of resamples S (default 100)
            %   input_parameters.bootstrapSeed: optional; see AGENTS.md on seeding
            %
            arguments
                ndi_calculator_obj
                parameters (1,1) struct {ndi.validators.mustHaveFields(parameters,{'input_parameters','depends_on'})}
            end

            speedtuning_bootstrap_calc = parameters;

            tuning_response_doc = ndi_calculator_obj.session.database_search(ndi.query('base.id', 'exact_string', ...
                did.db.struct_name_value_search(parameters.depends_on, 'stimulus_tuningcurve_id'), ''));
            if numel(tuning_response_doc) ~= 1
                error('Could not find stimulus tuning doc..');
            end
            tuning_response_doc = tuning_response_doc{1};

            min_xi = 0; max_xi = 1; numBootstrap = 100;
            if isfield(parameters.input_parameters, 'min_xi'), min_xi = parameters.input_parameters.min_xi; end
            if isfield(parameters.input_parameters, 'max_xi'), max_xi = parameters.input_parameters.max_xi; end
            if isfield(parameters.input_parameters, 'numBootstrap'), numBootstrap = parameters.input_parameters.numBootstrap; end

            app_doc = ndi_calculator_obj.newdocument();
            doc = ndi_calculator_obj.calculate_speed_indexes_bootstrap(tuning_response_doc, ...
                'min_xi', min_xi, 'max_xi', max_xi, 'numBootstrap', numBootstrap) + app_doc;

            if ~isempty(doc)
                doc = ndi.document(ndi_calculator_obj.doc_document_types{1}, 'speedtuning_bootstrap_calc', speedtuning_bootstrap_calc) + doc;
                doc = doc.setproperties('app', app_doc.document_properties.app);
                doc = doc.set_dependency_value('stimulus_tuningcurve_id', tuning_response_doc.id());
                doc = doc.set_dependency_value('element_id', tuning_response_doc.dependency_value('element_id'));
            end
        end % calculate

        function parameters = default_search_for_input_parameters(obj)
            % DEFAULT_SEARCH_FOR_INPUT_PARAMETERS - default search parameters
            parameters.input_parameters = struct('min_xi', 0, 'max_xi', 1, 'numBootstrap', 100);
            parameters.depends_on = did.datastructures.emptystruct('name', 'value');
            parameters.query = obj.default_parameters_query(parameters);
        end % default_search_for_input_parameters

        function query = default_parameters_query(obj, parameters_specification)
            % DEFAULT_PARAMETERS_QUERY - queries to find inputs (same inputs as speed_tuning)
            q1 = ndi.query('', 'isa', 'stimulus_tuningcurve', '');
            q2 = ndi.query('stimulus_tuningcurve.independent_variable_label', 'contains_string', 'spatial_frequency', '');
            q3 = ndi.query('stimulus_tuningcurve.independent_variable_label', 'contains_string', 'temporal_frequency', '');
            q_total = q1 & q2 & q3;
            query = struct('name', 'stimulus_tuningcurve_id', 'query', q_total);
        end % default_parameters_query()

        function b = is_valid_dependency_input(obj, name, value)
            % IS_VALID_DEPENDENCY_INPUT - accept any input (as speed_tuning does)
            b = 1;
        end % is_valid_dependency_input()

        function speed_props_doc = calculate_speed_indexes_bootstrap(obj, tuning_doc, kwargs)
            % CALCULATE_SPEED_INDEXES_BOOTSTRAP - bootstrap the Priebe fit over trials
            %
            % SPEED_PROPS_DOC = CALCULATE_SPEED_INDEXES_BOOTSTRAP(OBJ, TUNING_DOC, ...)
            %
            % Extracts the tuning curve (with its individual trial responses) the
            % same way ndi.calc.vis.speed_tuning.calculate_speed_indexes does,
            % then, for each of S = numBootstrap resamples, resamples the trials
            % of each stimulus condition with replacement, recomputes the
            % per-condition mean, and refits the three Priebe models. The S fits
            % are stacked so that a confidence interval can be taken over the
            % rows.
            %
            % Name/value:
            %   min_xi (0), max_xi (1), numBootstrap (100)
            %
            arguments
                obj
                tuning_doc
                kwargs.min_xi (1,1) double = 0
                kwargs.max_xi (1,1) double = 1
                kwargs.numBootstrap (1,1) double = 100
            end

            % --- Extract the tuning curve exactly as speed_tuning does -------
            % (kept in sync with ndi.calc.vis.speed_tuning.calculate_speed_indexes;
            %  a future refactor could factor this extraction into a shared helper.)
            properties.response_units = tuning_doc.document_properties.stimulus_tuningcurve.response_units;
            stim_response_doc = obj.session.database_search(ndi.query('base.id', ...
                'exact_string', tuning_doc.dependency_value('stimulus_response_scalar_id'), ''));
            if numel(stim_response_doc) ~= 1
                error('Could not find stimulus response scalar document.');
            end
            if iscell(stim_response_doc), stim_response_doc = stim_response_doc{1}; end
            properties.response_type = stim_response_doc.document_properties.stimulus_response_scalar.response_type;

            resp = ndi.app.stimulus.tuning_response.tuningcurvedoc2vhlabrespstruct(tuning_doc);

            sf = vlt.data.colvec(tuning_doc.document_properties.stimulus_tuningcurve.independent_variable_value(:, 1));
            tf = vlt.data.colvec(tuning_doc.document_properties.stimulus_tuningcurve.independent_variable_value(:, 2));

            % resp.ind{k} holds the individual single-trial responses for
            % condition k; this is the resampling unit for the bootstrap.
            ind = resp.ind;                     % cell array, one entry per condition
            nCond = numel(ind);

            SFs_grid = logspace(log10(0.01), log10(60),  14);
            TFs_grid = logspace(log10(0.01), log10(120), 14);
            [SFg, TFg] = meshgrid(SFs_grid, TFs_grid);

            P = 7;                              % number of Priebe parameters
            S = kwargs.numBootstrap;

            % Pre-allocate S-by-P / S-by-1 outputs (S rows = bootstrap samples).
            fit_params            = nan(S, P);
            fit_no_speed_params   = nan(S, P);
            fit_fullspeed_params  = nan(S, P);
            speed_tuning_index    = nan(S, 1);
            sf_preference         = nan(S, 1);
            tf_preference         = nan(S, 1);
            r_squared             = nan(S, 1);
            nested_F_no_speed_p   = nan(S, 1);
            nested_F_fullspeed_p  = nan(S, 1);

            for s = 1:S
                % --- resample trials within each condition, recompute the mean ---
                meanResp = nan(nCond, 1);
                for k = 1:nCond
                    tk = ind{k};
                    tk = tk(:);
                    if isempty(tk)
                        meanResp(k) = NaN;      % TODO: decide handling of empty conditions
                    else
                        meanResp(k) = mean(tk(randi(numel(tk), numel(tk), 1)));
                    end
                end

                % --- refit the three Priebe models on the resampled means -------
                [f_ns, sse_ns, ~]  = vis.speed.fit_nospeed(sf, tf, meanResp);
                [f_fs, sse_fs, ~]  = vis.speed.fit_fullspeed(sf, tf, meanResp);
                [f, sse_ws, r2_ws] = vis.speed.fit(sf, tf, meanResp, ...
                    kwargs.min_xi, kwargs.max_xi, 'SpecificStartPoint', [f_ns f_fs]);

                fit_params(s, :)           = f(:).';
                fit_no_speed_params(s, :)  = f_ns(:).';
                fit_fullspeed_params(s, :) = f_fs(:).';
                speed_tuning_index(s)      = f(3);
                sf_preference(s)           = f(6);
                tf_preference(s)           = f(7);
                r_squared(s)               = r2_ws;
                nested_F_no_speed_p(s)     = vis.speed.speed_nested_f(nCond, sse_ws, sse_ns);
                nested_F_fullspeed_p(s)    = vis.speed.speed_nested_f(nCond, sse_ws, sse_fs);
            end

            % --- assemble the result document -------------------------------
            % TODO (design issue): finalize exactly which fields are stored and
            % their orientation. Draft below stores per-sample matrices/vectors
            % plus the observed-data point estimate (via the parent calculator)
            % and percentile CIs. Confidence level is a parameter to decide.
            speed_tuning_bootstrap.properties   = properties;
            speed_tuning_bootstrap.input        = struct('numBootstrap', S, ...
                'min_xi', kwargs.min_xi, 'max_xi', kwargs.max_xi);
            speed_tuning_bootstrap.bootstrap    = struct( ...
                'Priebe_fit_parameters',            fit_params, ...           % S x P
                'Priebe_fit_no_speed_parameters',   fit_no_speed_params, ...  % S x P
                'Priebe_fit_fullspeed_parameters',  fit_fullspeed_params, ... % S x P
                'Priebe_fit_speed_tuning_index',    speed_tuning_index, ...   % S x 1
                'Priebe_fit_spatial_frequency_preference', sf_preference, ... % S x 1
                'Priebe_fit_temporal_frequency_preference', tf_preference, ...% S x 1
                'r_squared',                        r_squared, ...            % S x 1
                'nested_F_no_speed_p_value',        nested_F_no_speed_p, ...  % S x 1
                'nested_F_fullspeed_p_value',       nested_F_fullspeed_p);    % S x 1

            % TODO: also store the point-estimate fit (call the sibling
            % speed_tuning calculator on the same tuning_doc) and percentile CIs
            % so downstream code does not have to recompute them.

            speed_props_doc = ndi.document('speed_tuning_bootstrap', ...
                'speed_tuning_bootstrap', speed_tuning_bootstrap);
            speed_props_doc = speed_props_doc.set_dependency_value('element_id', ...
                tuning_doc.dependency_value('element_id'));
            speed_props_doc = speed_props_doc.set_dependency_value('stimulus_tuningcurve_id', tuning_doc.id());
        end % calculate_speed_indexes_bootstrap()

        % TODO: plot() override to show the bootstrap distribution of xi / a CI.
        % TODO: generate_mock_parameters() for the self-test framework (mirror
        %       speed_tuning; store expectations on the CIs, not point digits).

    end % methods()
end % speed_tuning_bootstrap
