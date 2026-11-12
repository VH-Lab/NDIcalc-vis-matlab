classdef speed_tuning_bootstrap < ndi.calc.tuning_fit
    % SPEED_TUNING_BOOTSTRAP - bootstrap the speed-tuning fit and derive coupling indexes
    %
    % SPEED_TUNING_BOOTSTRAP fits the Priebe et al. (2006) speed-tuning model
    % (vis.speed.tuningfunc) to a spatial x temporal frequency tuning curve, then
    % repeats the same free fit (vis.speed.fit, coupling index xi free) on
    % numBootstrap resamples of the single trials (vis.speed.bootstrap.fit).
    %
    % Every bootstrap fit is stored, in a binary ngrid file attached to the
    % output document (numBootstrap x 9: the 7 Priebe parameters, SSE, r^2), so
    % that any index can be recomputed later from the stored fits without
    % refitting. The document's JSON holds the point estimate and a summary of
    % indexes computed from the fits by vis.speed.bootstrap.indexes:
    %   - xi with its bootstrap confidence interval and p-values for xi == 0
    %     and xi == 1 (percentile bootstrap by default);
    %   - spatiotemporal inseparability by three rules (differs from 0,
    %     minimum effect |xi| > deltaInseparable, preferred-TF shift across the
    %     cell's SF passband > tfShiftOctaves);
    %   - speed tuning by three rules (Priebe et al. 2006: interval excludes 0
    %     and includes 1; equivalence within deltaSpeedTuned of 1; preferred
    %     speed changes by < speedConstancyOctaves across the passband);
    %   - preferred speed, SF passband, TF shift and speed change across the
    %     passband, each with intervals, and fit-quality flags.
    % The thresholds are input parameters and may be vectors, giving one
    % classification per value.
    %
    % Why a bootstrap: in simulations with the noise of recorded ferret V1
    % responses, the standard nested F test (as in ndi.calc.vis.speed_tuning)
    % rejected xi == 0 for ~15% of truly separable cells, while the percentile
    % bootstrap test from these fits rejected ~4% (nominal 5%).
    %
    % Output document type: 'speedtuning_bootstrap_calc' (a calculator document
    % that is also a 'speed_tuning_bootstrap' and an 'ngrid' document), with the
    % attached binary file 'speed_tuning_bootstrap_fits.ngrid'.
    %
    % See the documentation files for details:
    %   ndi.calculator.docfiletext('ndi.calc.vis.speed_tuning_bootstrap', 'general')
    %   ndi.calculator.docfiletext('ndi.calc.vis.speed_tuning_bootstrap', 'searching for inputs')
    %   ndi.calculator.docfiletext('ndi.calc.vis.speed_tuning_bootstrap', 'output')
    %
    % See also: ndi.calc.vis.speed_tuning, vis.speed.bootstrap.fit,
    %           vis.speed.bootstrap.indexes, vis.speed.bootstrap.xi_test

    properties (Constant)
        % Name of the attached binary file holding the bootstrap fits
        fitsFileName = 'speed_tuning_bootstrap_fits.ngrid';
        % Column names of the bootstrap fit matrix
        fitColumns = {'A', 'zeta', 'xi', 'sigma_sf', 'sigma_tf', 'sf0', 'tf0', 'sse', 'r_squared'};
    end

    methods
        function obj = speed_tuning_bootstrap(session)
            % SPEED_TUNING_BOOTSTRAP - create a speed_tuning_bootstrap calculator
            %
            % OBJ = ndi.calc.vis.speed_tuning_bootstrap(SESSION)
            %
            ndi.fun.checkCalcDirectory();

            obj = obj@ndi.calc.tuning_fit(session, 'speedtuning_bootstrap_calc', ...
                'speedtuning_bootstrap_calc');
            obj.defaultParametersCanFunction = true;
            obj.numberOfSelfTests = 4;
        end % speed_tuning_bootstrap()

        function doc = calculate(obj, parameters)
            % CALCULATE - run the bootstrap speed-tuning calculator on one input
            %
            % DOC = CALCULATE(OBJ, PARAMETERS)
            %
            % Creates a speedtuning_bootstrap_calc document for the
            % stimulus_tuningcurve named in PARAMETERS.depends_on. Missing
            % fields of PARAMETERS.input_parameters take the defaults of
            % ndi.calc.vis.speed_tuning_bootstrap.default_input_parameters.
            %
            arguments
                obj
                parameters (1,1) struct {ndi.validators.mustHaveFields(parameters,{'input_parameters','depends_on'})}
            end

            tuning_doc = obj.session.database_search(ndi.query('base.id', 'exact_string', ...
                did.db.struct_name_value_search(parameters.depends_on, 'stimulus_tuningcurve_id'), ''));
            if numel(tuning_doc) ~= 1
                error('ndi:calc:vis:speed_tuning_bootstrap:noInput', ...
                    'Could not find the stimulus_tuningcurve document.');
            end
            tuning_doc = tuning_doc{1};

            ip = ndi.calc.vis.speed_tuning_bootstrap.complete_input_parameters(parameters.input_parameters);
            parameters.input_parameters = ip;

            [result, B] = obj.calculate_bootstrap(tuning_doc, ip);

            ngridp = struct('data_size', 8, 'data_type', 'double', ...
                'data_dim', size(B), 'coordinates', []);

            app_doc = obj.newdocument();
            doc = ndi.document(obj.doc_document_types{1}, ...
                'speedtuning_bootstrap_calc', parameters, ...
                'speed_tuning_bootstrap', result, 'ngrid', ngridp) + app_doc;
            doc = doc.setproperties('app', app_doc.document_properties.app);
            doc = doc.set_dependency_value('stimulus_tuningcurve_id', tuning_doc.id());
            doc = doc.set_dependency_value('element_id', tuning_doc.dependency_value('element_id'));

            % Write the fits to a binary file and attach it. As in
            % ndi.calc.vis.hartley, the file is staged in a folder of the
            % session until the document is added to the database.
            stageDir = fullfile(obj.session.path, 'speed_tuning_bootstrap');
            if ~isfolder(stageDir)
                mkdir(stageDir);
            end
            fitsFile = fullfile(stageDir, [doc.id() '_' obj.fitsFileName]);
            ndi.fun.data.writengrid(B, fitsFile, 'double');
            doc = doc.add_file(obj.fitsFileName, fitsFile);
        end % calculate()

        function [result, B] = calculate_bootstrap(obj, tuning_doc, ip)
            % CALCULATE_BOOTSTRAP - point estimate, bootstrap fits and indexes for one tuning curve
            %
            % [RESULT, B] = CALCULATE_BOOTSTRAP(OBJ, TUNING_DOC, IP)
            %
            % RESULT is the speed_tuning_bootstrap property structure of the
            % output document; B is the numBootstrap x 9 matrix of bootstrap
            % fits. IP is a complete input_parameters structure (see
            % complete_input_parameters).
            %
            tc = vis.speed.extract_tuning_curve(obj.session, tuning_doc);
            sf = tc.spatial_frequency;
            tf = tc.temporal_frequency;
            trials = ndi.calc.vis.speed_tuning_bootstrap.trials_matrix(tc.resp.ind);

            if isempty(ip.bootstrapSeed)
                seed = 'shuffle';
            else
                seed = ip.bootstrapSeed;
            end
            % Use an already-open parallel pool if asked; never open one.
            runParallel = false;
            if ip.useParallel && exist('gcp', 'file')
                runParallel = ~isempty(gcp('nocreate'));
            end

            [B, info] = vis.speed.bootstrap.fit(sf, tf, trials, ip.numBootstrap, ...
                'min_xi', ip.min_xi, 'max_xi', ip.max_xi, ...
                'pointEstimateStartPoints', ip.pointEstimateStartPoints, ...
                'numberStartPoints', ip.numberStartPoints, ...
                'scaleResiduals', logical(ip.scaleResiduals), ...
                'RandomSeed', seed, 'useParallel', runParallel);

            pe = info.pointEstimate;
            ix = vis.speed.bootstrap.indexes(B, pe, 'alpha', ip.alpha, 'method', ip.method, ...
                'deltaInseparable', ip.deltaInseparable, 'deltaSpeedTuned', ip.deltaSpeedTuned, ...
                'tfShiftOctaves', ip.tfShiftOctaves, ...
                'speedConstancyOctaves', ip.speedConstancyOctaves, ...
                'min_xi', ip.min_xi, 'max_xi', ip.max_xi, ...
                'sfRange', [min(sf) max(sf)]);

            result = struct();
            result.properties = tc.properties;
            result.significance = tc.significance;
            result.point_estimate = struct( ...
                'Priebe_fit_parameters', pe(:).', ...
                'Priebe_fit_speed_tuning_index', pe(3), ...
                'Priebe_fit_spatial_frequency_preference', pe(6), ...
                'Priebe_fit_temporal_frequency_preference', pe(7), ...
                'sse', info.pointEstimateSSE, ...
                'r_squared', info.pointEstimateR2);
            result.bootstrap = struct( ...
                'numBootstrap', ip.numBootstrap, ...
                'numValid', ix.quality.numBootstrap, ...
                'columns', {obj.fitColumns}, ...
                'trials_per_condition_median', median(sum(~isnan(trials), 2)), ...
                'stimulus_spatial_frequency_range', [min(sf) max(sf)], ...
                'stimulus_temporal_frequency_range', [min(tf) max(tf)]);
            result.indexes = ndi.calc.vis.speed_tuning_bootstrap.indexes_for_document(ix);
        end % calculate_bootstrap()

        function parameters = default_search_for_input_parameters(obj)
            % DEFAULT_SEARCH_FOR_INPUT_PARAMETERS - default search and input parameters
            %
            % PARAMETERS = DEFAULT_SEARCH_FOR_INPUT_PARAMETERS(OBJ)
            %
            % PARAMETERS.input_parameters holds every input parameter at its
            % default (see default_input_parameters); PARAMETERS.query finds
            % stimulus_tuningcurve documents that vary spatial and temporal
            % frequency, as for ndi.calc.vis.speed_tuning.
            %
            parameters.input_parameters = ndi.calc.vis.speed_tuning_bootstrap.default_input_parameters();
            parameters.depends_on = did.datastructures.emptystruct('name', 'value');
            parameters.query = obj.default_parameters_query(parameters);
        end % default_search_for_input_parameters()

        function query = default_parameters_query(obj, parameters_specification)
            % DEFAULT_PARAMETERS_QUERY - query for inputs (same inputs as speed_tuning)
            %
            % QUERY = DEFAULT_PARAMETERS_QUERY(OBJ, PARAMETERS_SPECIFICATION)
            %
            % Finds stimulus_tuningcurve documents whose independent variables
            % include spatial_frequency and temporal_frequency.
            %
            q1 = ndi.query('', 'isa', 'stimulus_tuningcurve', '');
            q2 = ndi.query('stimulus_tuningcurve.independent_variable_label', 'contains_string', 'spatial_frequency', '');
            q3 = ndi.query('stimulus_tuningcurve.independent_variable_label', 'contains_string', 'temporal_frequency', '');
            query = struct('name', 'stimulus_tuningcurve_id', 'query', q1 & q2 & q3);
        end % default_parameters_query()

        function b = is_valid_dependency_input(obj, name, value)
            % IS_VALID_DEPENDENCY_INPUT - accept any input found by the query (as speed_tuning does)
            b = 1;
        end % is_valid_dependency_input()

        function B = read_bootstrap_fits(obj, doc)
            % READ_BOOTSTRAP_FITS - read the stored bootstrap fits of an output document
            %
            % B = READ_BOOTSTRAP_FITS(OBJ, DOC)
            %
            % Returns the numBootstrap x 9 matrix stored in the document's
            % attached ngrid file; columns are OBJ.fitColumns
            % (A zeta xi sigma_sf sigma_tf sf0 tf0 sse r_squared). DOC must
            % be in the session's database.
            %
            dims = doc.document_properties.ngrid.data_dim(:).';
            dtype = doc.document_properties.ngrid.data_type;
            % database_openbinarydoc returns an already-open, read-only file
            % object; read it directly (as ndi.calc.vis.hartley does).
            % ndi.fun.data.readngrid cannot be used here: it reopens the
            % object with mode 'r', which a read-only file object rejects
            % (it requires 'rb').
            f = obj.session.database_openbinarydoc(doc, obj.fitsFileName);
            try
                B = fread(f, prod(dims), dtype);
            catch err
                obj.session.database_closebinarydoc(f);
                rethrow(err);
            end
            obj.session.database_closebinarydoc(f);
            if numel(B) ~= prod(dims)
                error('ndi:calc:vis:speed_tuning_bootstrap:shortFile', ...
                    'Read %d values from %s; expected %d.', numel(B), obj.fitsFileName, prod(dims));
            end
            B = reshape(B, dims);
        end % read_bootstrap_fits()

        function ix = recompute_indexes(obj, doc, options)
            % RECOMPUTE_INDEXES - recompute the indexes from stored fits with new thresholds
            %
            % IX = RECOMPUTE_INDEXES(OBJ, DOC, 'Name', Value, ...)
            %
            % Reads the stored bootstrap fits of DOC and recomputes the indexes
            % with vis.speed.bootstrap.indexes, without refitting. Thresholds
            % not given take the values stored in DOC's input parameters.
            % Name/value pairs: alpha, method, deltaInseparable,
            % deltaSpeedTuned, tfShiftOctaves, speedConstancyOctaves (vectors
            % allowed for the four thresholds).
            %
            % Example:
            %   ix = c.recompute_indexes(doc, 'deltaSpeedTuned', [0.25 0.5 1]);
            %   ix.speedTuned.equivalence   % one value per margin
            %
            arguments
                obj
                doc (1,1) ndi.document
                options.alpha = []
                options.method = []
                options.deltaInseparable = []
                options.deltaSpeedTuned = []
                options.tfShiftOctaves = []
                options.speedConstancyOctaves = []
            end
            ip = ndi.calc.vis.speed_tuning_bootstrap.complete_input_parameters( ...
                doc.document_properties.speedtuning_bootstrap_calc.input_parameters);
            names = fieldnames(options);
            for k = 1:numel(names)
                if ~isempty(options.(names{k}))
                    ip.(names{k}) = options.(names{k});
                end
            end
            B = obj.read_bootstrap_fits(doc);
            pe = doc.document_properties.speed_tuning_bootstrap.point_estimate.Priebe_fit_parameters;
            sfRange = doc.document_properties.speed_tuning_bootstrap.bootstrap.stimulus_spatial_frequency_range;
            ix = vis.speed.bootstrap.indexes(B, pe, 'alpha', ip.alpha, 'method', ip.method, ...
                'deltaInseparable', ip.deltaInseparable, 'deltaSpeedTuned', ip.deltaSpeedTuned, ...
                'tfShiftOctaves', ip.tfShiftOctaves, 'speedConstancyOctaves', ip.speedConstancyOctaves, ...
                'min_xi', ip.min_xi, 'max_xi', ip.max_xi, 'sfRange', sfRange(:).');
        end % recompute_indexes()

        function h = plot(obj, doc_or_parameters, varargin)
            % PLOT - bootstrap distribution of the coupling index xi
            %
            % H = PLOT(OBJ, DOC, ...)
            %
            % Histogram of the bootstrap estimates of xi, with the point
            % estimate (black), the confidence interval (red dashed), and
            % reference lines at xi = 0 (separable) and xi = 1 (speed tuned).
            % The title gives xi, its interval, and the inseparable and
            % speed-tuned classifications. DOC must be in the database.
            % Handles are returned in H.
            %
            h = plot@ndi.calculator(obj, doc_or_parameters, varargin{:});
            if ~isa(doc_or_parameters, 'ndi.document')
                error('ndi:calc:vis:speed_tuning_bootstrap:plotNeedsDoc', ...
                    'PLOT needs an ndi.document output of this calculator.');
            end
            doc = doc_or_parameters;
            sb = doc.document_properties.speed_tuning_bootstrap;
            B = obj.read_bootstrap_fits(doc);
            xi = B(:, 3);
            xi = xi(~isnan(xi));
            ci = sb.indexes.xi.ci;
            xiHat = sb.indexes.xi.estimate;

            hold on;
            h.histogram = histogram(xi, 'Normalization', 'probability');
            yl = ylim;
            h.ci_low  = plot([ci(1) ci(1)], yl, 'r--', 'linewidth', 1.5);
            h.ci_high = plot([ci(2) ci(2)], yl, 'r--', 'linewidth', 1.5);
            h.point_estimate = plot([xiHat xiHat], yl, 'k-', 'linewidth', 2);
            h.ref0 = plot([0 0], yl, 'b:');
            h.ref1 = plot([1 1], yl, 'b:');

            if ~h.params.suppress_x_label
                h.xlabel = xlabel('\xi (spatiotemporal coupling index)');
            end
            if ~h.params.suppress_y_label
                h.ylabel = ylabel('bootstrap probability');
            end
            if ~h.params.suppress_title
                h.title = title(sprintf(['\\xi = %.3g [%.3g, %.3g] (%g%% CI); ' ...
                    'inseparable: %d; speed tuned (Priebe): %d'], xiHat, ci(1), ci(2), ...
                    100 * (1 - sb.indexes.parameters.alpha), ...
                    sb.indexes.inseparable.differentFrom0, sb.indexes.speedTuned.priebe));
            end
            box off;
        end % plot()

        % TESTING METHODS

        function [param_struct, independent_variable, x, r] = generate_mock_parameters(obj, scope, index)
            % GENERATE_MOCK_PARAMETERS - mock tuning curves for the self-tests
            %
            % [PARAM_STRUCT, INDEPENDENT_VARIABLE, X, R] = GENERATE_MOCK_PARAMETERS(OBJ, SCOPE, INDEX)
            %
            % Four model cells (vis.speed.tuningfunc) on the spatial x temporal
            % frequency grid used by ndi.calc.vis.speed_tuning's mocks:
            %   1  separable                         (xi = 0)
            %   2  speed tuned                       (xi = 1)
            %   3  negative ("anti-speed") coupling  (xi = -0.5)
            %   4  preferred SF below the lowest stimulus SF (sf0 = 0.03 vs
            %      0.05), xi = 0.5: the low-SF edge case, where the fitted sf0
            %      sits at the edge of the sampled range
            % The framework (ndi.calc.tuning_fit.generate_mock_docs) adds
            % trial-to-trial noise and replicate trials, which are the
            % resampling unit of the bootstrap. INDEX selects the cell
            % (wrapped with MOD); SCOPE does not change the cells.
            %
            %          cell:  1      2      3      4
            A        = [   10,    10,    10,     5 ];
            zeta     = [    0,     0,     0,     0 ];
            xi       = [    0,     1,  -0.5,   0.5 ];
            sigma_sf = [  0.3,   0.3,   0.3,   0.3 ];   % log10 units (SF passband ~2.3 octaves)
            sigma_tf = [  0.4,   0.4,   0.4,   0.4 ];   % log10 units
            sf0      = [  0.2,   0.2,   0.2,  0.03 ];
            tf0      = [    4,     4,     4,     4 ];

            P_ = [A(:) zeta(:) xi(:) sigma_sf(:) sigma_tf(:) sf0(:) tf0(:)];
            P = P_(1 + mod(index - 1, size(P_, 1)), :);

            sfs = [0.05 0.08 0.1 0.2 0.4 0.8 1.2];
            tfs = [0.5 1 2 4 8 16 32];
            [SFs, TFs] = meshgrid(sfs, tfs);
            r_ = vis.speed.tuningfunc(SFs, TFs, P);

            param_struct = struct('contrast', .5);
            independent_variable = {'temporal_frequency', 'spatial_frequency'};
            x = [SFs(:), TFs(:)];
            r = r_(:);

            % blank (control) stimulus with response 0
            x(end + 1, :) = NaN;
            r(end + 1, 1) = 0;
        end % generate_mock_parameters()

    end % methods

    methods (Static)

        function ip = default_input_parameters()
            % DEFAULT_INPUT_PARAMETERS - every input parameter at its default value
            %
            % IP = ndi.calc.vis.speed_tuning_bootstrap.default_input_parameters()
            %
            % See the 'searching for inputs' documentation for the meaning of
            % each field.
            %
            ip = struct( ...
                'min_xi', -4, ...
                'max_xi', 4, ...
                'numBootstrap', 200, ...
                'pointEstimateStartPoints', 40, ...
                'numberStartPoints', 5, ...
                'scaleResiduals', 1, ...
                'bootstrapSeed', [], ...
                'useParallel', 1, ...
                'alpha', 0.05, ...
                'method', 'percentile', ...
                'deltaInseparable', 0.2, ...
                'deltaSpeedTuned', 0.5, ...
                'tfShiftOctaves', 1, ...
                'speedConstancyOctaves', 1);
        end % default_input_parameters()

        function ip = complete_input_parameters(ip_in)
            % COMPLETE_INPUT_PARAMETERS - fill in missing input parameters with defaults
            %
            % IP = ndi.calc.vis.speed_tuning_bootstrap.complete_input_parameters(IP_IN)
            %
            % Returns the default input parameters with every field present in
            % IP_IN taking IP_IN's value. Vector thresholds are made rows.
            %
            ip = ndi.calc.vis.speed_tuning_bootstrap.default_input_parameters();
            if isempty(ip_in)
                return;
            end
            names = fieldnames(ip_in);
            for k = 1:numel(names)
                ip.(names{k}) = ip_in.(names{k});
            end
            for f = {'deltaInseparable', 'deltaSpeedTuned', 'tfShiftOctaves', 'speedConstancyOctaves'}
                ip.(f{1}) = ip.(f{1})(:).';
            end
            ip.method = char(ip.method);
            if ~(ip.min_xi < 0 && ip.max_xi > 1)
                warning('ndi:calc:vis:speed_tuning_bootstrap:xiRange', ...
                    ['min_xi/max_xi = [%g %g] do not contain 0 and 1 in their interior; ' ...
                    'the tests of xi == 0 and xi == 1 are not meaningful at a bound.'], ...
                    ip.min_xi, ip.max_xi);
            end
        end % complete_input_parameters()

        function trials = trials_matrix(ind)
            % TRIALS_MATRIX - single-trial responses as an nCond x maxTrials matrix
            %
            % TRIALS = ndi.calc.vis.speed_tuning_bootstrap.trials_matrix(IND)
            %
            % IND is the cell array of single-trial responses (one cell per
            % condition, resp.ind of the vhlab response structure). Rows are
            % NaN-padded to the largest number of trials.
            %
            nCond = numel(ind);
            nMax = max([cellfun(@numel, ind(:)); 0]);
            trials = nan(nCond, nMax);
            for k = 1:nCond
                v = ind{k}(:).';
                trials(k, 1:numel(v)) = v;
            end
        end % trials_matrix()

        function s = indexes_for_document(ix)
            % INDEXES_FOR_DOCUMENT - indexes struct with JSON-friendly types
            %
            % S = ndi.calc.vis.speed_tuning_bootstrap.indexes_for_document(IX)
            %
            % Converts the output of vis.speed.bootstrap.indexes for storage:
            % logical classifications become double 0/1 (so that they can be
            % compared numerically in the self-tests), and the method string
            % and sfRange are kept as stored.
            %
            s = ix;
            groups = {'inseparable', 'speedTuned'};
            for g = 1:numel(groups)
                names = fieldnames(s.(groups{g}));
                for k = 1:numel(names)
                    v = s.(groups{g}).(names{k});
                    if islogical(v)
                        s.(groups{g}).(names{k}) = double(v);
                    end
                end
            end
        end % indexes_for_document()

    end % methods (Static)
end % speed_tuning_bootstrap
