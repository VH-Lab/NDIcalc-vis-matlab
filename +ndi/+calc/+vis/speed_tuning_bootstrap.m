classdef speed_tuning_bootstrap < ndi.calc.tuning_fit
    % SPEED_TUNING_BOOTSTRAP - bootstrap version of the speed_tuning calculator
    %
    % SPEED_TUNING_BOOTSTRAP computes the same spatiotemporal-coupling fits as
    % ndi.calc.vis.speed_tuning, but repeats the fit over S bootstrap resamples
    % of the individual trials, so a confidence interval can be placed on the
    % coupling index xi and the other Priebe fit parameters.
    %
    % Motivation (Suarez-Casanova et al. review, R1.4/R1.5): the nested F test
    % calls a site "speed tuned" whenever it *fails to reject* xi==1, which
    % conflates "the data support xi==1" with "the data do not constrain xi at
    % all" -- some such cells have a fitted xi<0. A confidence interval on xi
    % (an equivalence-style classification) separates genuinely speed-tuned
    % cells from unconstrained ones. The tuning-curve document already stores
    % the individual trial responses (resp.ind), so a nonparametric bootstrap
    % over trials is feasible.
    %
    % Output document type: 'speedtuning_bootstrap_calc'
    % Result document type: 'speed_tuning_bootstrap'
    %
    % The per-sample free-fit parameters are stored as an S-by-P matrix (S
    % bootstrap samples in rows, P=7 Priebe parameters in columns), so that
    % prctile(params,[2.5 97.5]) returns a 2-by-P per-parameter confidence
    % interval with no transpose.
    %
    % See also: ndi.calc.vis.speed_tuning, vis.speed.extract_tuning_curve,
    %           vis.speed.fit_priebe_triplet, vis.speed.fit

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

            % The field-by-field self-tests (generate_mock_parameters + stored
            % expected documents) require the expected mock documents to be
            % generated in a licensed-MATLAB session, as for every other
            % calculator here. Until those are generated and committed, this
            % stays 0. See generate_mock_parameters below.
            obj.numberOfSelfTests = 0;
        end % speed_tuning_bootstrap()

        function doc = calculate(ndi_calculator_obj, parameters)
            % CALCULATE - run the bootstrap speed-tuning calculator
            %
            % DOC = CALCULATE(NDI_CALCULATOR_OBJ, PARAMETERS)
            %
            % Creates a speedtuning_bootstrap_calc document. Input parameters
            % mirror speed_tuning (min_xi, max_xi) plus:
            %   input_parameters.numBootstrap    : number of resamples S (default 200)
            %   input_parameters.confidenceLevel : percent for the CIs (default 95)
            %   input_parameters.bootstrapSeed   : optional; [] draws fresh (see
            %                                      AGENTS.md on seeding)
            %   input_parameters.useParallel     : use an already-open parallel
            %                                      pool if one exists (default true;
            %                                      never opens a pool)
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

            min_xi = 0; max_xi = 1; numBootstrap = 200; confidenceLevel = 95;
            bootstrapSeed = []; useParallel = true;
            ip = parameters.input_parameters;
            if isfield(ip, 'min_xi'),          min_xi = ip.min_xi;                   end
            if isfield(ip, 'max_xi'),          max_xi = ip.max_xi;                   end
            if isfield(ip, 'numBootstrap'),    numBootstrap = ip.numBootstrap;       end
            if isfield(ip, 'confidenceLevel'), confidenceLevel = ip.confidenceLevel; end
            if isfield(ip, 'bootstrapSeed'),   bootstrapSeed = ip.bootstrapSeed;     end
            if isfield(ip, 'useParallel'),     useParallel = ip.useParallel;         end

            app_doc = ndi_calculator_obj.newdocument();
            doc = ndi_calculator_obj.calculate_speed_indexes_bootstrap(tuning_response_doc, ...
                'min_xi', min_xi, 'max_xi', max_xi, 'numBootstrap', numBootstrap, ...
                'confidenceLevel', confidenceLevel, 'bootstrapSeed', bootstrapSeed, ...
                'useParallel', useParallel) + app_doc;

            if ~isempty(doc)
                doc = ndi.document(ndi_calculator_obj.doc_document_types{1}, 'speedtuning_bootstrap_calc', speedtuning_bootstrap_calc) + doc;
                doc = doc.setproperties('app', app_doc.document_properties.app);
                doc = doc.set_dependency_value('stimulus_tuningcurve_id', tuning_response_doc.id());
                doc = doc.set_dependency_value('element_id', tuning_response_doc.dependency_value('element_id'));
            end
        end % calculate

        function parameters = default_search_for_input_parameters(obj)
            % DEFAULT_SEARCH_FOR_INPUT_PARAMETERS - default search parameters
            parameters.input_parameters = struct('min_xi', 0, 'max_xi', 1, ...
                'numBootstrap', 200, 'confidenceLevel', 95);
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

        function h = plot(obj, doc_or_parameters, varargin)
            % PLOT - show the bootstrap distribution of the speed index xi
            %
            % H = PLOT(OBJ, DOC_OR_PARAMETERS, ...)
            %
            % Plots a histogram of the S bootstrap estimates of xi, with the
            % point estimate and the confidence interval marked, and reference
            % lines at xi==0 (no speed tuning) and xi==1 (full speed tuning) so
            % the equivalence-style classification reads off the figure.
            %
            % Handles are returned in H.

            h = plot@ndi.calculator(obj, doc_or_parameters, varargin{:});

            if isa(doc_or_parameters, 'ndi.document')
                doc = doc_or_parameters;
            else
                error('Do not know how to proceed without an ndi document for doc_or_parameters.');
            end

            sb = doc.document_properties.speed_tuning_bootstrap;
            xi = sb.bootstrap.Priebe_fit_parameters(:, 3);
            ci = sb.confidence_interval;
            xi_hat = sb.point_estimate.fit.Priebe_fit_speed_tuning_index;

            hold on;
            h.histogram = histogram(xi, 'Normalization', 'probability');
            yl = ylim;

            h.ci_low  = plot([ci.speed_tuning_index_ci(1) ci.speed_tuning_index_ci(1)], yl, 'r--', 'linewidth', 1.5);
            h.ci_high = plot([ci.speed_tuning_index_ci(2) ci.speed_tuning_index_ci(2)], yl, 'r--', 'linewidth', 1.5);
            h.point_estimate = plot([xi_hat xi_hat], yl, 'k-', 'linewidth', 2);
            h.ref0 = plot([0 0], yl, 'b:');
            h.ref1 = plot([1 1], yl, 'b:');

            if ~h.params.suppress_x_label
                h.xlabel = xlabel('\xi (speed tuning index)');
            end
            if ~h.params.suppress_y_label
                h.ylabel = ylabel('bootstrap probability');
            end
            if ~h.params.suppress_title
                h.title = title(['\xi = ' num2str(xi_hat, 3) ' [' ...
                    num2str(ci.speed_tuning_index_ci(1), 3) ', ' ...
                    num2str(ci.speed_tuning_index_ci(2), 3) '] (' ...
                    num2str(ci.confidenceLevel) '% CI)']);
            end
            box off;

        end % plot()

        function speed_props_doc = calculate_speed_indexes_bootstrap(obj, tuning_doc, kwargs)
            % CALCULATE_SPEED_INDEXES_BOOTSTRAP - bootstrap the Priebe fit over trials
            %
            % SPEED_PROPS_DOC = CALCULATE_SPEED_INDEXES_BOOTSTRAP(OBJ, TUNING_DOC, ...)
            %
            % Extracts the tuning curve (with its individual trial responses)
            % via vis.speed.extract_tuning_curve, then, for each of S =
            % numBootstrap resamples, resamples the trials of each stimulus
            % condition with replacement, recomputes the per-condition mean, and
            % refits the three Priebe models via vis.speed.fit_priebe_triplet.
            % The S free-fit parameter vectors are stacked S-by-P so a
            % confidence interval can be taken over the rows.
            %
            % Name/value:
            %   min_xi (0), max_xi (1), numBootstrap (200), confidenceLevel (95),
            %   bootstrapSeed ([] = fresh draw), useParallel (true)
            %
            arguments
                obj
                tuning_doc
                kwargs.min_xi (1,1) double = 0
                kwargs.max_xi (1,1) double = 1
                kwargs.numBootstrap (1,1) double = 200
                kwargs.confidenceLevel (1,1) double = 95
                kwargs.bootstrapSeed = []
                kwargs.useParallel (1,1) logical = true
            end

            % --- Extract the tuning curve (shared with speed_tuning) ---------
            tc = vis.speed.extract_tuning_curve(obj.session, tuning_doc);
            properties = tc.properties;
            sf  = tc.spatial_frequency;
            tf  = tc.temporal_frequency;
            ind = tc.resp.ind;                       % cell, one entry per condition
            nCond = numel(ind);
            observed_mean = vlt.data.colvec(tc.resp.curve(2, :));   % per-condition mean

            % --- Point estimate: the observed-data fit, computed by the sibling
            %     speed_tuning calculator on the same document, so the point
            %     estimate cannot drift from speed_tuning. -------------------
            st = ndi.calc.vis.speed_tuning(obj.session);
            pe_doc = st.calculate_speed_indexes(tuning_doc, 'min_xi', kwargs.min_xi, 'max_xi', kwargs.max_xi);
            point_estimate = pe_doc.document_properties.speed_tuning;

            % --- Private stream for the trial resampling (see AGENTS.md) -----
            if isempty(kwargs.bootstrapSeed)
                rs = vis.randomstream('shuffle');
            else
                rs = vis.randomstream(kwargs.bootstrapSeed);
            end

            S = kwargs.numBootstrap;
            P = 7;                                   % number of Priebe parameters
            min_xi = kwargs.min_xi;
            max_xi = kwargs.max_xi;

            % --- Serial pre-pass: draw every resampled per-condition mean up
            %     front, from the private stream, so the result does not depend
            %     on whether the fit loop below runs serially or in parallel.
            %     A condition with no stored trials keeps its observed mean
            %     (rather than injecting NaN, which would break the fit). -----
            meanRespAll = repmat(observed_mean(:).', S, 1);   % S x nCond
            for s = 1:S
                for k = 1:nCond
                    tk = ind{k}(:);
                    nk = numel(tk);
                    if nk > 0
                        idx = randi(rs, nk, nk, 1);
                        meanRespAll(s, k) = mean(tk(idx));
                    end
                end
            end

            fit_params = nan(S, P);
            r_squared  = nan(S, 1);

            % Use an already-open parallel pool if the caller asked for parallel
            % and one exists; never open a pool here (gcp('nocreate') returns []
            % when there is no pool and does not create one).
            runParallel = false;
            if kwargs.useParallel && exist('gcp', 'file')
                runParallel = ~isempty(gcp('nocreate'));
            end

            if runParallel
                parfor s = 1:S
                    [f, ~, ~, stats] = vis.speed.fit_priebe_triplet(sf, tf, meanRespAll(s, :).', min_xi, max_xi);
                    fit_params(s, :) = f(:).';
                    r_squared(s)     = stats.r_squared;
                end
            else
                for s = 1:S
                    [f, ~, ~, stats] = vis.speed.fit_priebe_triplet(sf, tf, meanRespAll(s, :).', min_xi, max_xi);
                    fit_params(s, :) = f(:).';
                    r_squared(s)     = stats.r_squared;
                end
            end

            % --- Percentile confidence intervals over the S rows ------------
            alphaPct  = (100 - kwargs.confidenceLevel) / 2;
            ci        = prctile(fit_params, [alphaPct, 100 - alphaPct], 1);   % 2 x P
            ci_median = prctile(fit_params, 50, 1);                           % 1 x P

            confidence_interval = struct( ...
                'confidenceLevel',                  kwargs.confidenceLevel, ...
                'parameter_ci_low',                 ci(1, :), ...
                'parameter_ci_high',                ci(2, :), ...
                'parameter_median',                 ci_median, ...
                'speed_tuning_index_ci',            [ci(1, 3) ci(2, 3)], ...
                'speed_tuning_index_median',        ci_median(3), ...
                'spatial_frequency_preference_ci',  [ci(1, 6) ci(2, 6)], ...
                'temporal_frequency_preference_ci', [ci(1, 7) ci(2, 7)]);

            % --- Assemble the result document -------------------------------
            speed_tuning_bootstrap.properties = properties;
            speed_tuning_bootstrap.input = struct('numBootstrap', S, ...
                'min_xi', min_xi, 'max_xi', max_xi, 'confidenceLevel', kwargs.confidenceLevel);
            speed_tuning_bootstrap.point_estimate = point_estimate;
            speed_tuning_bootstrap.bootstrap = struct( ...
                'Priebe_fit_parameters', fit_params, ...   % S x P
                'r_squared',             r_squared);       % S x 1
            speed_tuning_bootstrap.confidence_interval = confidence_interval;

            speed_props_doc = ndi.document('speed_tuning_bootstrap', ...
                'speed_tuning_bootstrap', speed_tuning_bootstrap);
            speed_props_doc = speed_props_doc.set_dependency_value('element_id', ...
                tuning_doc.dependency_value('element_id'));
            speed_props_doc = speed_props_doc.set_dependency_value('stimulus_tuningcurve_id', tuning_doc.id());
        end % calculate_speed_indexes_bootstrap()

        % TESTING METHODS

        function [param_struct, independent_variable, x, r] = generate_mock_parameters(obj, scope, index)
            % GENERATE_MOCK_PARAMETERS - generate mock parameters for testing
            %
            % [PARAM_STRUCT, INDEPENDENT_VARIABLE, X, R] = GENERATE_MOCK_PARAMETERS(OBJ, SCOPE, INDEX)
            %
            % Generates a parameter set for a mock speed-tuning document, mirroring
            % ndi.calc.vis.speed_tuning.generate_mock_parameters. Three cells are
            % defined so the bootstrap self-test can check the confidence interval
            % behaves as expected across the regimes the paper cares about:
            %   1  speed-tuned, well constrained (xi = 1)
            %   2  not speed tuned              (xi ~ 0)
            %   3  weak response / poorly constrained (small A, xi = 0.5)
            %
            % SCOPE can be 'standard', 'random_nonoise', or 'random_noisy'; the
            % framework adds trial-to-trial noise and replicate trials (the
            % resampling unit for the bootstrap). INDEX selects the cell (1..TOTAL,
            % wrapped with MOD).
            %
            % The field-by-field self-tests are enabled by generating the expected
            % mock documents in a licensed-MATLAB session and raising
            % numberOfSelfTests; see the class constructor.

            %          cell:  1(tuned) 2(not)   3(weak)
            A        = [   5,        5,      1 ];   % peak response
            zeta     = [   0,        0,      0 ];   % temporal-frequency skew
            xi       = [   1,   0.0001,    0.5 ];   % speed tuning index (0..1)
            sigma_sf = [   1,        1,      1 ];   % spatial-frequency tuning width
            sigma_tf = [   1,        1,      1 ];   % temporal-frequency tuning width
            sf0      = [ 0.2,      0.2, sqrt(2)/5 ];% preferred spatial frequency
            tf0      = [   2,        2,      4 ];   % preferred temporal frequency

            P_ = [A(:) zeta(:) xi(:) sigma_sf(:) sigma_tf(:) sf0(:) tf0(:)];
            total = size(P_, 1);

            actual_index = 1 + mod(index - 1, total);

            % no dependence on scope for this stimulus type
            P = P_(actual_index, :);

            % grid of stimulus conditions (taken from speed_tuning's mock / the demo)
            sfs = [0.05 0.08 0.1 0.2 0.4 0.8 1.2];
            tfs = [0.5 1 2 4 8 16 32];
            [SFs, TFs] = meshgrid(sfs, tfs);
            function_params = P;
            r_ = vlt.neuro.vision.speed.tuningfunc(SFs, TFs, function_params);

            param_struct = struct('contrast', .5);
            independent_variable = {'temporal_frequency', 'spatial_frequency'};
            x = [SFs(:), TFs(:)];
            r = r_(:);

            % blank (control) stimulus with firing rate 0
            x(end + 1, :) = NaN;
            r(end + 1, 1) = 0;

        end % generate_mock_parameters

    end % methods()
end % speed_tuning_bootstrap
