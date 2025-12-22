classdef contrast_sensitivity < ndi.calculator

    methods
        function contrast_sensitivity_obj = contrast_sensitivity(session)
            % CONTRAST_TUNING - a contrast_sensitivity demonstration of an ndi.calculator object
            %
            % CONTRAST_TUNING_OBJ = CONTRAST_TUNING(SESSION)
            %
            % Creates a CONTRAST_TUNING ndi.calculator object
            %
                contrast_sensitivity_obj = contrast_sensitivity_obj@ndi.calculator(session,'contrastsensitivity_calc',...
                    'contrastsensitivity_calc');
                contrast_sensitivity_obj.numberOfSelfTests = 1;
        end % contrast_sensitivity()

        function [docs, doc_output, doc_expected_output] = generate_mock_docs(obj, scope, number_of_tests, kwargs)
            % GENERATE_MOCK_DOCS - generate mock documents for testing
            %
            % [DOCS, DOC_OUTPUT, DOC_EXPECTED_OUTPUT] = GENERATE_MOCK_DOCS(OBJ, SCOPE, NUMBER_OF_TESTS, ...)
            %
            % Generates mock documents and runs the calculator for testing purposes.
            %
            % Inputs:
            %   OBJ - the calculator object
            %   SCOPE - 'highSNR' or 'lowSNR'
            %   NUMBER_OF_TESTS - number of tests to run (only 1 supported for now)
            %   KWARGS - keyword arguments:
            %     'generate_expected_docs' (boolean, default false) - if true, save expected output
            %     'specific_test_inds' (vector, default []) - run only specific test indices
            %

            arguments
                obj
                scope (1,:) char {mustBeMember(scope,{'highSNR','lowSNR'})}
                number_of_tests
                kwargs.generate_expected_docs (1,1) logical = false
                kwargs.specific_test_inds (1,:) double = []
            end

            docs = cell(1,number_of_tests);
            doc_output = cell(1,number_of_tests);
            doc_expected_output = cell(1,number_of_tests);

            % Define Mock Directory
            p = fileparts(mfilename('fullpath'));
            mock_dir = fullfile(p, 'mock', 'contrast_sensitivity');
            if ~isfolder(mock_dir)
                mkdir(mock_dir);
            end

            % Define Noise based on Scope
            if strcmpi(scope, 'lowSNR')
                noise_level = 1.0;
            else
                noise_level = 0.1;
            end

            for i = 1:number_of_tests
                if ~isempty(kwargs.specific_test_inds) && ~ismember(i, kwargs.specific_test_inds)
                    continue;
                end

                % Use Object Session
                S = obj.session;

                % Mock Data Generation
                % 1. Element
                mock_data = ndi.mock.fun.subject_stimulator_neuron(S);
                nde = mock_data.mock_spikes;
                nde_stim = mock_data.mock_stimulator;
                refNum = nde.reference;
                nde_control = ndi.element.timeseries(S, 'mock control', refNum, 'stimulus_control', [], 0, mock_data.mock_subject.id());

                current_docs = {mock_data.mock_subject, nde_stim.load_element_doc(), nde.load_element_doc(), nde_control.load_element_doc()};

                % Parameters
                sFrequencies = [0.05, 0.1, 0.2];
                contrasts = 0:0.1:1;
                Rm_base = 10;
                b = 2;
                c50 = 0.3;

                % Adjust parameters based on SF
                % SF 0.05: High response
                % SF 0.1: Lower response
                % SF 0.2: Very low response

                scale_factors = [1, 0.5, 0.1];

                for j = 1:numel(sFrequencies)
                    sf = sFrequencies(j);
                    Rm = Rm_base * scale_factors(j);

                    % Naka Rushton Response
                    resp_mean = Rm * (contrasts.^b) ./ (c50^b + contrasts.^b);

                    % 2. Stimulus Presentation
                    stim_params = struct();
                    stim_params.stimuli(1).parameters.sFrequency = sf;
                    stim_params.stimuli(1).parameters.contrast = contrasts; % Optional but good for completeness
                    stim_params.stimuli(2).parameters.is_blank = 1;
                    stim_params.presentation_order = 1;

                    stim_pres_doc = ndi.document('stimulus_presentation', 'stimulus_presentation', stim_params) + ...
                        obj.session.newdocument();
                    stim_pres_doc = stim_pres_doc.set_dependency_value('stimulus_element_id', nde_stim.id());
                    S.database_add(stim_pres_doc);
                    current_docs{end+1} = stim_pres_doc;

                    % 3. Stimulus Response Scalar
                    stim_resp_param_struct = struct();
                    stim_resp_param_doc = ndi.document('stimulus_response_scalar_parameters', 'stimulus_response_scalar_parameters', stim_resp_param_struct) + ...
                        obj.session.newdocument();
                    S.database_add(stim_resp_param_doc);
                    current_docs{end+1} = stim_resp_param_doc;

                    stim_resp_struct.response_type = 'mean';
                    stim_resp_struct.responses = [];
                    stim_resp_scalar_doc = ndi.document('stimulus_response_scalar', 'stimulus_response_scalar', stim_resp_struct) + ...
                        obj.session.newdocument();
                    stim_resp_scalar_doc = stim_resp_scalar_doc.set_dependency_value('element_id', nde.id());
                    stim_resp_scalar_doc = stim_resp_scalar_doc.set_dependency_value('stimulus_presentation_id', stim_pres_doc.id());
                    stim_resp_scalar_doc = stim_resp_scalar_doc.set_dependency_value('stimulator_id', nde_stim.id());
                    stim_resp_scalar_doc = stim_resp_scalar_doc.set_dependency_value('stimulus_control_id', nde_control.id());
                    stim_resp_scalar_doc = stim_resp_scalar_doc.set_dependency_value('stimulus_response_scalar_parameters_id', stim_resp_param_doc.id());
                    S.database_add(stim_resp_scalar_doc);
                    current_docs{end+1} = stim_resp_scalar_doc;

                    % 4. Stimulus Tuning Curve
                    n_trials = 5;
                    tuning_struct.independent_variable_label = {'Contrast'};
                    tuning_struct.independent_variable_value = contrasts;
                    tuning_struct.stimid = ones(size(contrasts));
                    tuning_struct.stimulus_presentation_number = 1;
                    tuning_struct.response_units = 'Hz';
                    tuning_struct.response_mean = resp_mean;
                    tuning_struct.response_stddev = resp_mean * noise_level; % Use noise level
                    tuning_struct.response_stderr = resp_mean * (noise_level * 0.5);

                    tuning_struct.individual_responses_real = repmat(resp_mean, n_trials, 1) + randn(n_trials, numel(resp_mean)) * noise_level;
                    tuning_struct.individual_responses_imaginary = 0 * tuning_struct.individual_responses_real;

                    tuning_struct.control_stimid = 2;
                    tuning_struct.control_response_mean = zeros(size(resp_mean));
                    tuning_struct.control_response_stddev = zeros(size(resp_mean));
                    tuning_struct.control_response_stderr = zeros(size(resp_mean));
                    tuning_struct.control_individual_responses_real = zeros(size(tuning_struct.individual_responses_real));
                    tuning_struct.control_individual_responses_imaginary = zeros(size(tuning_struct.individual_responses_real));

                    stim_tuning_doc = ndi.document('stimulus_tuningcurve', 'stimulus_tuningcurve', tuning_struct) + ...
                        obj.session.newdocument();
                    stim_tuning_doc = stim_tuning_doc.set_dependency_value('element_id', nde.id());
                    stim_tuning_doc = stim_tuning_doc.set_dependency_value('stimulus_response_scalar_id', stim_resp_scalar_doc.id());
                    S.database_add(stim_tuning_doc);
                    current_docs{end+1} = stim_tuning_doc;
                end

                docs{i} = current_docs;

                % Run Contrast Tuning Calculator
                ct = ndi.calc.vis.contrast_tuning(S);
                ct.run('Replace');

                % Run Contrast Sensitivity Calculator (Target)
                ccc = ndi.calc.vis.contrast_sensitivity(S);

                parameters = ccc.default_search_for_input_parameters();
                parameters.query.query = ndi.query('element.type','exact_string','spikes','');

                ccc.run('Replace', parameters);

                % Get Output
                q = ndi.query('','isa','contrastsensitivity_calc','');
                out_docs = S.database_search(q);

                if ~isempty(out_docs)
                    doc_output{i} = out_docs{1};
                else
                    error('No output document generated.');
                end

                % Handle Expected Docs
                 if kwargs.generate_expected_docs
                    doc_expected_output{i} = doc_output{i};
                    vlt.file.str2text(fullfile(mock_dir, ['mock.' int2str(i) '.json']), jsonencode(doc_expected_output{i}.document_properties));
                else
                    expected_file = fullfile(mock_dir, ['mock.' int2str(i) '.json']);
                    if isfile(expected_file)
                         doc_expected_output{i} = ndi.document(jsondecode(fileread(expected_file)));
                    else
                         doc_expected_output{i} = [];
                    end
                end
            end
        end

        function doc = calculate(ndi_calculator_obj, parameters)
            % CALCULATE - perform the calculator for ndi.calc.example.contrast_sensitivity
            %
            % DOC = CALCULATE(NDI_CALCULATOR_OBJ, PARAMETERS)
            %
            % Creates a contrast_sensitivity_calc document given input parameters.
            %
            % The document that is created contrast_sensitivity
            % by the input parameters.
                arguments
                    ndi_calculator_obj
                    parameters (1,1) struct {ndi.validators.mustHaveFields(parameters,{'input_parameters','depends_on'})}
                end

                % Step 1: set up the output structure
                contrastsensitivity_calc = parameters;

                element_doc = ndi_calculator_obj.session.database_search(ndi.query('base.id','exact_string',...
                    did.db.struct_name_value_search(parameters.depends_on,'element_id'),''));
                if numel(element_doc)~=1
                    error('Could not find element doc..');
                end
                element_doc = element_doc{1};

                % Step 2: search for possible stimulus_response_scalar docs

                q1a = ndi.query('','depends_on','element_id',element_doc.id());
                q1b = ndi.query('','isa','stimulus_response_scalar','');
                stim_resp_scalar = ndi_calculator_obj.session.database_search(q1a&q1b);

                tuning_curve_app = ndi.calc.stimulus.tuningcurve(ndi_calculator_obj.session);

                % Step 3: find out how many stimulus presentations we have here; could be several
                stim_pres_ids_all = {};

                for i=1:numel(stim_resp_scalar)
                    stim_pres_ids_all{i} = stim_resp_scalar{i}.dependency_value('stimulus_presentation_id');
                end

                stim_pres_id = unique(stim_pres_ids_all);

                % Now, if possible, construct a constrast sensitivity curve for each stimulus presentation

                doc = {};

                for i=1:numel(stim_pres_id)
                    % now see if the stimulus presentations vary in contrast and spatial frequency
                    q2 = ndi.query('base.id','exact_string',stim_pres_id{i},'');
                    stim_pres_doc = ndi_calculator_obj.session.database_search(q2);
                    if numel(stim_pres_doc) ~=1
                        error(['Missing stimulus presentation document for ' stim_pres_id{i} '. (Should not happen).']);
                    end
                    stim_pres_doc = stim_pres_doc{1};
                    good = 0;
                    stim_resp_doc_indexes = find(strcmp(stim_pres_id{i},stim_pres_ids_all)); % find all the responses associated with this stimulus presentation
                    v1 = tuning_curve_app.property_value_array(stim_resp_scalar{stim_resp_doc_indexes(1)},'contrast');
                    if numel(v1)>2
                        good = 1;
                    end
                    if good
                        v2 = tuning_curve_app.property_value_array(stim_resp_scalar{stim_resp_doc_indexes(1)},'sFrequency');
                    else
                        v2 = {};
                    end
                    if numel(v2)<=2
                        good = 0;
                    end
                    if good
                        if numel(stim_resp_doc_indexes)==1% we're okay, just use the one choice
                            stim_resp_index_value = 1;
                        else
                            [b,r,dummy,dummy,mean_i,mod_i] = ndi.app.stimulus.tuning_response.modulated_or_mean(stim_resp_scalar(stim_resp_doc_indexes));
                            if b==-1
                                warning(['Skipping responses to stimulus presentation ' stim_pres_id{i} ' because we found more than one response and could not identify mean and modulated response.']);
                                stim_resp_index_value = [];
                            elseif b==0
                                stim_resp_index_value = mean_i;
                                response_type_here = 'mean';
                            elseif b==1
                                stim_resp_index_value = mod_i;
                                response_type_here = 'F1';
                            end
                        end

                        if ~isempty(stim_resp_index_value)
                            % Step 3: Search for contrast tuning curve objects that depend on this stimulus response document
                            q3 = ndi.query('stimulus_tuningcurve.independent_variable_label','contains_string','Contrast','');
                            q4 = ndi.query('','depends_on','stimulus_response_scalar_id',stim_resp_scalar{stim_resp_index_value}.id());
                            q5 = ndi.query('','isa','stimulus_tuningcurve','');
                            tuning_curves = ndi_calculator_obj.session.database_search(q3 & q4 & q5);

                            spatial_frequencies = [];
                            sensitivity_RB = [];
                            sensitivity_RBN = [];
                            sensitivity_RBNS = [];
                            response_type = stim_resp_scalar{stim_resp_index_value}.document_properties.stimulus_response_scalar.response_type;

                            relative_max_gain_RB = [];
                            relative_max_gain_RBN = [];
                            relative_max_gain_RBNS = [];

                            empirical_c50_RB = [];
                            empirical_c50_RBN = [];
                            empirical_c50_RBNS = [];

                            saturation_index_RB = [];
                            saturation_index_RBN = [];
                            saturation_index_RBNS = [];

                            parameters_RB = [];
                            parameters_RBN = [];
                            parameters_RBNS = [];

                            fitless_interpolated_c50 = [];

                            visual_response_p = [];
                            across_stims_p = [];

                            all_contrast_tuning_curves_ids = {};

                            for k=1:numel(tuning_curves)
                                q6 = ndi.query('','isa','contrast_tuning','');
                                q7 = ndi.query('','depends_on','stimulus_tuningcurve_id',tuning_curves{k}.id());

                                contrast_tuning_props = ndi_calculator_obj.session.database_search(q6&q7);

                                if numel(contrast_tuning_props)>1
                                    error('Found multiple contrast tuning property records for a single tuning curve.');
                                elseif numel(contrast_tuning_props)==0
                                    error(['Found contrast tuning curve but no contrast tuning curve properties for element with id ' element_doc.id()]);
                                else
                                    contrast_tuning_props = contrast_tuning_props{1};
                                end
                                all_contrast_tuning_curves_ids{end+1} = contrast_tuning_props.id();
                                
                                stimid = tuning_curves{k}.document_properties.stimulus_tuningcurve.stimid(1);
                                params_here = stim_pres_doc.document_properties.stimulus_presentation.stimuli(stimid).parameters;
                                if isfield(params_here,'sFrequency')
                                    spatial_frequencies(end+1) = getfield(params_here,'sFrequency');
                                else
                                    error('Expected spatial frequency information.'); % should this be an error or just a skip?
                                end

                                ctp = contrast_tuning_props.document_properties.contrast_tuning;

                                sensitivity_RB = [ sensitivity_RB vlt.data.colvec(ctp.fit.naka_rushton_RB_sensitivity) ];
                                sensitivity_RBN = [ sensitivity_RBN vlt.data.colvec(ctp.fit.naka_rushton_RBN_sensitivity) ];
                                sensitivity_RBNS = [ sensitivity_RBNS vlt.data.colvec(ctp.fit.naka_rushton_RBNS_sensitivity) ];

                                relative_max_gain_RB = [ relative_max_gain_RB vlt.data.colvec(ctp.fit.naka_rushton_RB_relative_max_gain)];
                                relative_max_gain_RBN = [ relative_max_gain_RBN vlt.data.colvec(ctp.fit.naka_rushton_RBN_relative_max_gain)];
                                relative_max_gain_RBNS = [ relative_max_gain_RBNS vlt.data.colvec(ctp.fit.naka_rushton_RBNS_relative_max_gain)];

                                empirical_c50_RB = [ empirical_c50_RB vlt.data.colvec(ctp.fit.naka_rushton_RB_empirical_c50)];
                                empirical_c50_RBN = [ empirical_c50_RBN vlt.data.colvec(ctp.fit.naka_rushton_RBN_empirical_c50)];
                                empirical_c50_RBNS = [ empirical_c50_RBNS vlt.data.colvec(ctp.fit.naka_rushton_RBNS_empirical_c50)];

                                saturation_index_RB = [ saturation_index_RB vlt.data.colvec(ctp.fit.naka_rushton_RB_empirical_c50)];
                                saturation_index_RBN = [ saturation_index_RBN vlt.data.colvec(ctp.fit.naka_rushton_RBN_empirical_c50)];
                                saturation_index_RBNS = [ saturation_index_RBNS vlt.data.colvec(ctp.fit.naka_rushton_RBNS_empirical_c50)];

                                fitless_interpolated_c50 = [ fitless_interpolated_c50 vlt.data.colvec(ctp.fitless.interpolated_c50)];

                                parameters_RB = [ parameters_RB vlt.data.colvec(ctp.fit.naka_rushton_RB_parameters)];
                                parameters_RBN = [ parameters_RBN vlt.data.colvec(ctp.fit.naka_rushton_RBN_parameters)];
                                parameters_RBNS = [ parameters_RBNS vlt.data.colvec(ctp.fit.naka_rushton_RBNS_parameters)];

                                visual_response_p(end+1) = contrast_tuning_props.document_properties.contrast_tuning.significance.visual_response_anova_p;
                                across_stims_p(end+1) = contrast_tuning_props.document_properties.contrast_tuning.significance.across_stimuli_anova_p;
                            end

                            [spatial_frequencies,order] = sort(spatial_frequencies);
                            sensitivity_RB = sensitivity_RB(:,order);
                            sensitivity_RBN = sensitivity_RBN(:,order);
                            sensitivity_RBNS = sensitivity_RBNS(:,order);
                            relative_max_gain_RB = relative_max_gain_RB(:,order);
                            relative_max_gain_RBN = relative_max_gain_RBN(:,order);
                            relative_max_gain_RBNS = relative_max_gain_RBNS(:,order);
                            empirical_c50_RB = empirical_c50_RB(:,order);
                            empirical_c50_RBN = empirical_c50_RBN(:,order);
                            empirical_c50_RBNS = empirical_c50_RBNS(:,order);
                            saturation_index_RB = saturation_index_RB(:,order);
                            saturation_index_RBN = saturation_index_RBN(:,order);
                            saturation_index_RBNS = saturation_index_RBNS(:,order);
                            fitless_interpolated_c50 = fitless_interpolated_c50(:,order);
                            parameters_RB = parameters_RB(:,order);
                            parameters_RBN = parameters_RBN(:,order);
                            parameters_RBNS = parameters_RBNS(:,order);

                            % make the doc

                            parameters_here = contrastsensitivity_calc;
                            parameters_here.spatial_frequencies = vlt.data.rowvec(spatial_frequencies);
                            parameters_here.sensitivity_RB = sensitivity_RB;
                            parameters_here.sensitivity_RBN = sensitivity_RBN;
                            parameters_here.sensitivity_RBNS = sensitivity_RBNS;
                            parameters_here.relative_max_gain_RB = relative_max_gain_RB;
                            parameters_here.relative_max_gain_RBN = relative_max_gain_RBN;
                            parameters_here.relative_max_gain_RBNS = relative_max_gain_RBNS;
                            parameters_here.empirical_c50_RB = empirical_c50_RB;
                            parameters_here.empirical_c50_RBN = empirical_c50_RBN;
                            parameters_here.empirical_c50_RBNS = empirical_c50_RBNS;
                            parameters_here.saturation_index_RB = saturation_index_RB;
                            parameters_here.saturation_index_RBN = saturation_index_RBN;
                            parameters_here.saturation_index_RBNS = saturation_index_RBNS;
                            parameters_here.fitless_interpolated_c50 = fitless_interpolated_c50;
                            parameters_here.parameters_RB = parameters_RB;
                            parameters_here.parameters_RBN = parameters_RBN;
                            parameters_here.parameters_RBNS = parameters_RBNS;
                            parameters_here.is_modulated_response = b;
                            % could actually do 2-factor ANOVA on responses; would be better
                            parameters_here.visual_response_p_bonferroni = nanmin(visual_response_p)*numel(visual_response_p);
                            parameters_here.response_varies_p_bonferroni = nanmin(across_stims_p)*numel(across_stims_p);
                            parameters_here.depends_on = did.datastructures.emptystruct('name','value');
                            parameters_here.response_type = response_type_here;

                            if numel(tuning_curves)>0
                                doc{end+1} = ndi.document(ndi_calculator_obj.doc_document_types{1},...
                                    'contrastsensitivity_calc',parameters_here) + ndi_calculator_obj.newdocument();
                                doc{end} = doc{end}.set_dependency_value('element_id',element_doc.id());
                                doc{end} = doc{end}.set_dependency_value('stimulus_presentation_id', stim_pres_id{i});
                                doc{end} = doc{end}.set_dependency_value('stimulus_response_scalar_id',...
                                    stim_resp_scalar{stim_resp_index_value}.id());
                                for q_=1:numel(all_contrast_tuning_curves_ids)
                                    doc{end} = doc{end}.add_dependency_value_n('contrasttuning_id',all_contrast_tuning_curves_ids{q_});
                                end
                            end
                        end % if stim_resp_index_value is not empty
                    end % if good
                end

        end % calculate

        function parameters = default_search_for_input_parameters(ndi_calculator_obj)
            % DEFAULT_SEARCH_FOR_INPUT_PARAMETERS - default parameters for searching for inputs
            %
            % PARAMETERS = DEFAULT_SEARCH_FOR_INPUT_PARAMETERS(NDI_CALCULATOR_OBJ)
            %
            % Returns a list of the default search parameters for finding appropriate inputs
            % to the calculator. For contrast_sensitivity_calc, there is no appropriate default parameters
            % so this search will yield empty.
            %
                parameters.input_parameters = struct([]);
                parameters.depends_on = did.datastructures.emptystruct('name','value');
                parameters.query = ndi_calculator_obj.default_parameters_query(parameters);

        end % default_search_for_input_parameters

                function query = default_parameters_query(ndi_calculator_obj, parameters_specification)
            % DEFAULT_PARAMETERS_QUERY - what queries should be used to search for input parameters if none are provided?
            %
            % QUERY = DEFAULT_PARAMETERS_QUERY(NDI_CALCULATOR_OBJ, PARAMETERS_SPECIFICATION)
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
            % For the ndi.calc.stimulus.contrast_sensitivity_calc class, this looks for
            % documents of type 'stimulus_response_scalar' with 'response_type' fields
            % the contain 'mean' or 'F1'.
            %
            %
                q_total = ndi.query('','isa','base','');

                query = struct('name','element_id','query',q_total);
        end % default_parameters_query()

        function h=plot(ndi_calculator_obj, doc_or_parameters, varargin)
                        % PLOT - provide a diagnostic plot to show the results of the calculator
                        %
                        % H=PLOT(NDI_CALCULATOR_OBJ, DOC_OR_PARAMETERS, ...)
                        %
                        % Produce a plot of the tuning curve.
            %
                        % Handles to the figure, the axes, and any objects created are returned in H.
                        %
                        % This function takes additional input arguments as name/value pairs.
                        % See ndi.calculator.plot_parameters for a description of those parameters.

                % call superclass plot method to set up axes
                h=plot@ndi.calculator(ndi_calculator_obj, doc_or_parameters, varargin{:});

                if isa(doc_or_parameters,'ndi.document')
                    doc = doc_or_parameters;
                else
                    error('Do not know how to proceed without an ndi document for doc_or_parameters.');
                end

                cs = doc.document_properties.contrastsensitivity_calc; % shorten our typing

                noise_threshold_indexes = [5];

                % First plot responses

                line_type = '-';
                if ~(cs.visual_response_p_bonferroni<0.05)
                    line_type = '--';
                end

                for i=1:numel(noise_threshold_indexes)
                    hold on;
                    h_baseline = plot([min(cs.spatial_frequencies) max(cs.spatial_frequencies)],...
                        [0 0],'k--','linewidth',1.0001);
                    h_baseline.Annotation.LegendInformation.IconDisplayStyle = 'off';
                    h.objects(end+1) = h_baseline;

                    % plot all fits

                    h.objects(end+1) = plot(cs.spatial_frequencies, cs.sensitivity_RB(noise_threshold_indexes(i),:), ...
                        ['o' line_type], 'color', (1/(max(1,numel(noise_threshold_indexes)-1))) * [0 1 0],'linewidth',1.5);
                    h.objects(end+1) = plot(cs.spatial_frequencies, cs.sensitivity_RBN(noise_threshold_indexes(i),:), ...
                        ['d' line_type], 'color', (1/(max(1,numel(noise_threshold_indexes)-1))) * [0 0 1],'linewidth',1.5);
                    h.objects(end+1) = plot(cs.spatial_frequencies, cs.sensitivity_RBNS(noise_threshold_indexes(i),:), ...
                        ['s' line_type], 'color', (1/(max(1,numel(noise_threshold_indexes)-1))) * [1 0 1],'linewidth',1.5);
                end

                if ~h.params.suppress_x_label
                    h.xlabel = xlabel('Spatial frequency');
                end
                if ~h.params.suppress_y_label
                    h.ylabel = ylabel(['Sensitivity']);
                end

                if 0% when database is faster :-/
                    if ~h.params.suppress_title
                        element = ndi.database.fun.ndi_document2ndi_object(doc.dependency_value('element_id'),ndi_calculator_obj.session);
                        h.title = title(element.elementstring(), 'interp','none');
                    end
                end
                set(gca,'xscale','log');
                box off;

        end % plot()

    end % methods()
end % contrast_sensitivity
