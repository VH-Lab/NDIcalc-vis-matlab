classdef speed_tuning < ndi.calculator

	methods
		function speed_tuning_obj = speed_tuning(session)
			% SPEED_TUNING - a speed_tuning demonstration of an ndi.calculator object
			%
			% SPEED_TUNING_OBJ = SPEED_TUNING(SESSION)
			%
			% Creates a SPEED_TUNING ndi.calculator object
			%
				ndi.globals;
				w = which('ndi.calc.vis.contrast_tuning');
				parparparpar = fileparts(fileparts(fileparts(fileparts(w))));                
				speed_tuning_obj = speed_tuning_obj@ndi.calculator(session,'speedtuning_calc',...
					fullfile(parparparpar,'ndi_common','database_documents','calc','speedtuning_calc.json'));
		end; % speed_tuning()

		function doc = calculate(ndi_calculator_obj, parameters)
			% CALCULATE - perform the calculator for ndi.calc.example.speed_tuning
			%
			% DOC = CALCULATE(NDI_CALCULATOR_OBJ, PARAMETERS)
			%
			% Creates a speed_tuning_calc document given input parameters.
			%
			% The document that is created speed_tuning
			% by the input parameters.
				% check inputs
				if ~isfield(parameters,'input_parameters'), error(['parameters structure lacks ''input_parameters''.']); end;
				if ~isfield(parameters,'depends_on'), error(['parameters structure lacks ''depends_on''.']); end;
				
				% Step 1: set up the output structure
				speed_tuning_calc = parameters;

				tuning_response_doc = ndi_calculator_obj.session.database_search(ndi.query('base.id','exact_string',...
					vlt.db.struct_name_value_search(parameters.depends_on,'stimulus_tuningcurve_id'),''));
				if numel(tuning_response_doc)~=1, 
					error(['Could not find stimulus tuning doc..']);
				end;
				tuning_response_doc = tuning_response_doc{1};

				% Step 2: perform the calculator, which here creates a speed_tuning doc
				doc = ndi_calculator_obj.calculate_speed_indexes(tuning_response_doc) + ...
					ndi_calculator_obj.newdocument();
				
				if ~isempty(doc), 
					doc = ndi.document(ndi_calculator_obj.doc_document_types{1},'speedtuning_calc',speed_tuning_calc) + doc;
					doc = doc.set_dependency_value('stimulus_tuningcurve_id',tuning_response_doc.id());
					doc = doc.set_dependency_value('element_id',tuning_response_doc.dependency_value('element_id'));
				end;
		end; % calculate

		function parameters = default_search_for_input_parameters(ndi_calculator_obj)
			% DEFAULT_SEARCH_FOR_INPUT_PARAMETERS - default parameters for searching for inputs
			%
			% PARAMETERS = DEFAULT_SEARCH_FOR_INPUT_PARAMETERS(NDI_CALCULATOR_OBJ)
			%
			% Returns a list of the default search parameters for finding appropriate inputs
			% to the calculator. For speed_tuning_calc, there is no appropriate default parameters
			% so this search will yield empty.
			%
				parameters.input_parameters = struct([]);
				parameters.depends_on = vlt.data.emptystruct('name','value');
				parameters.query = ndi_calculator_obj.default_parameters_query(parameters);
					
		end; % default_search_for_input_parameters

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
			% For the ndi.calc.stimulus.speed_tuning_calc class, this looks for 
			% documents of type 'stimulus_response_scalar.json' with 'response_type' fields
			% the contain 'mean' or 'F1'.
			%
			%
				q1 = ndi.query('','isa','stimulus_tuningcurve','');
				q2 = ndi.query('stimulus_tuningcurve.independent_variable_label','contains_string','spatial_frequency','');
				q3 = ndi.query('stimulus_tuningcurve.independent_variable_label','contains_string','temporal_frequency','');
				%q4 = ndi.query('stimulus_tuningcurve.independent_variable_label','hassize',[2 1],'');
				q_total = q1 & q2 & q3; % & q4;

				query = struct('name','stimulus_tuningcurve_id','query',q_total);
		end; % default_parameters_query()

		function b = is_valid_dependency_input(ndi_calculator_obj, name, value)
			% IS_VALID_DEPENDENCY_INPUT - is a potential dependency input actually valid for this calculator?
			%
			% B = IS_VALID_DEPENDENCY_INPUT(NDI_CALCULATOR_OBJ, NAME, VALUE)
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
	
				switch lower(name),
					case lower('stimulus_tuningcurve_id'),
						q = ndi.query('base.id','exact_string',value,'');
						d = ndi_calculator_obj.S.database_search(q);
						b = (numel(d.document_properties.independent_variable_label) ==2);
					case lower('element_id'),
						b = 1;
				end;
		end; % is_valid_dependency_input()

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

				if isa(doc_or_parameters,'ndi.document'),
					doc = doc_or_parameters;
				else,
					error(['Do not know how to proceed without an ndi document for doc_or_parameters.']);
				end;

				sp = doc.document_properties.speed_tuning; % shorten our typing
				tc = sp.tuning_curve; % shorten our typing
				ft = sp.fit;

				% First plot fit
				hold on;

				%h_baseline = plot([min(tc.speed) max(tc.speed)],...
				%	[0 0],'k--','linewidth',1.0001);
				%h_baseline.Annotation.LegendInformation.IconDisplayStyle = 'off';

				% now call the plot routine

				[SF,TF,MNs] = vlt.math.vector2mesh(tc.spatial_frequency,tc.temporal_frequency,tc.mean);
				MNs_fit = vlt.neuro.vision.speed.tuningfunc(SF,TF,ft.Priebe_fit_parameters);

				significant = 0;
				linestyle = '--';
				if sp.significance.visual_response_anova_p<0.05,
					significant = 1;
					linestyle = '-';
				end;
				vlt.neuro.vision.speed.plottuning(SF,TF,MNs_fit,'marker','none','linestyle',linestyle);

				% now plot raw responses
				vlt.neuro.vision.speed.plottuning(SF,TF,MNs);

                ch = get(gcf,'children');
                currentaxes = gca;
                axes(ch(1));
                title(['Speed tuning:' num2str(ft.Priebe_fit_parameters(3)) ', pref: ' num2str(ft.Priebe_fit_parameters(7)/ft.Priebe_fit_parameters(6))]);

                element_str = '';
                ele = ndi_calculator_obj.session.database_search(ndi.query('base.id','exact_string',doc.dependency_value('element_id')));
                epoch_id = '';
                tc_id = ndi_calculator_obj.session.database_search(ndi.query('base.id','exact_string',doc.dependency_value('stimulus_tuningcurve_id')));
                if ~isempty(tc_id),
                    srs = ndi_calculator_obj.session.database_search(ndi.query('base.id','exact_string',tc_id{1}.dependency_value('stimulus_response_scalar_id')));
                    if ~isempty(srs),
                        epoch_id = [srs{1}.document_properties.stimulus_response.element_epochid ' ' srs{1}.document_properties.stimulus_response_scalar.response_type];
                    end;
                end;

                if ~isempty(ele),
                    element = ndi.database.fun.ndi_document2ndi_object(ele{1},ndi_calculator_obj.session);
                    element_str = element.elementstring();
                    axes(ch(2));
                    title([element_str ' ' epoch_id],'interp','none');
                end;
                axes(currentaxes);

				if 0, % plot function already does this
				if ~h.params.suppress_x_label,
					h.xlabel = xlabel('Speed (deg/sec)');
				end;
				if ~h.params.suppress_y_label,
					h.ylabel = ylabel(['Response (' sp.properties.response_type ', ' sp.properties.response_units ')']);
				end;
				box off;
				end;


		end; % plot()

		function speed_props_doc = calculate_speed_indexes(ndi_calculator_obj, tuning_doc)
			% CALCULATE_SPEED_INDEXES - calculate speed index values from a tuning curve
			%
			% SPEED_PROPS_DOC = CALCULATE_SPEED_INDEXES(NDI_SPEED_TUNING_CALC_OBJ, TUNING_DOC)
			%
			% Given a 2-dimensional tuning curve document with measurements at many spatial and
			% and temporal frequencies, this function calculates speed response
			% parameters and stores them in SPEED_TUNING document SPEED_PROPS_DOC.
			%
			%
				properties.response_units = tuning_doc.document_properties.stimulus_tuningcurve.response_units;
				
				stim_response_doc = ndi_calculator_obj.session.database_search(ndi.query('base.id',...
					'exact_string',tuning_doc.dependency_value('stimulus_response_scalar_id'),''));
				if numel(stim_response_doc)~=1,
					error(['Could not find stimulus response scalar document.']);
				end;
				if iscell(stim_response_doc),
					stim_response_doc = stim_response_doc{1};
				end;

				properties.response_type = stim_response_doc.document_properties.stimulus_response_scalar.response_type;

				sf_coord = 1;
				tf_coord = 2;
				if contains(tuning_doc.document_properties.stimulus_tuningcurve.independent_variable_label{1},'temporal','IgnoreCase',true),
					sf_coord = 2;
					tf_coord = 1;
				end;

				resp = ndi.app.stimulus.tuning_response.tuningcurvedoc2vhlabrespstruct(tuning_doc);

				[anova_across_stims, anova_across_stims_blank] = neural_response_significance(resp);

				tuning_curve = struct(...
					'spatial_frequency', ...
						vlt.data.rowvec(tuning_doc.document_properties.stimulus_tuningcurve.independent_variable_value(:,1)), ...
					'temporal_frequency', ...
						vlt.data.rowvec(tuning_doc.document_properties.stimulus_tuningcurve.independent_variable_value(:,2)), ...
					'mean', resp.curve(2,:), ...
					'stddev', resp.curve(3,:), ...
					'stderr', resp.curve(4,:), ...
					'individual', {resp.ind}, ...
					'control_stddev', resp.blankresp(2),...
					'control_stderr', resp.blankresp(3));

				significance = struct('visual_response_anova_p',anova_across_stims_blank,...
					'across_stimuli_anova_p', anova_across_stims);

				f = vlt.neuro.vision.speed.fit(tuning_curve.spatial_frequency(:),tuning_curve.temporal_frequency(:),tuning_curve.mean(:));
				sfs = logspace(0.01,60,200);
				tfs = logspace(0.01,120,200);
				[SFs,TFs] = meshgrid(sfs,tfs);
				fit_values = vlt.neuro.vision.speed.tuningfunc(SFs(:),TFs(:),f);

				fit.Priebe_fit_parameters = f;
				fit.Priebe_fit_spatial_frequencies = SFs(:);
				fit.Priebe_fit_temporal_frequencies = TFs(:);
				fit.Priebe_fit_values = fit_values;
				fit.Priebe_fit_speed_tuning_index = fit.Priebe_fit_parameters(3);
				fit.Priebe_fit_spatial_frequency_preference = fit.Priebe_fit_parameters(6);
				fit.Priebe_fit_temporal_frequency_preference = fit.Priebe_fit_parameters(7);

				speed_tuning.properties = properties;
				speed_tuning.tuning_curve = tuning_curve;
				speed_tuning.significance = significance;
				speed_tuning.fit = fit;

				speed_props_doc = ndi.document('speed_tuning',...
					'speed_tuning',speed_tuning);
				speed_props_doc = speed_props_doc.set_dependency_value('element_id', ...
					tuning_doc.dependency_value('element_id'));
				speed_props_doc = speed_props_doc.set_dependency_value('stimulus_tuningcurve_id',tuning_doc.id());
		end; % calculate_speed_indexes()

        % TESTING METHODS

        function [docs, doc_output, doc_expected_output] = generate_mock_docs(speed_calc_obj, scope, number_of_tests, varargin)
			% GENERATE_MOCK_DOCS - generate mock documents and expected answers for tests
			%
			% [DOCS, DOC_OUTPUT, DOC_EXPECTED_OUTPUT] = GENERATE_MOCK_DOCS(SPEED_CALC_OBJ, ...
			%    SCOPE, NUMBER_OF_TESTS, ...)
			%
			% Creates a set of documents to test ndi.calc.vis.speed_tuning.
			%
			% SCOPE is the scope to be tested: 'standard', 'low_noise', 'high_noise'
			% NUMBER_OF_TESTS indicates the number of tests to be performed.
			%
			% DOCS{i} is the set of helper documents that may have been created
			%   in generating the ith test.
			% DOC_OUTPUT{i} is the actual output of the calculator when operating on
			%   DOCS{i} (the ith test).
			% DOC_EXPECTED_OUTPUT{i} is what the output of the calculator should be, if there
			%   were no noise.
			%
			% The quality of these outputs are evaluted using the function COMPARE_MOCK_DOCS
			% as part of the TEST function for ndi.calculator objects.
			%
			% This function's behavior can be modified by name/value pairs.
			% --------------------------------------------------------------------------------
			% | Parameter (default):     | Description:                                      |
			% |--------------------------|---------------------------------------------------|
			% | generate_expected_docs(0)| Should we generate the expected docs? (That is,   |
			% |                          |   generate the "right answer"?) Use carefully.    |
			% |--------------------------|---------------------------------------------------|
			%

				generate_expected_docs = 0;
				vlt.data.assign(varargin{:});

				docs = {};
				doc_output = {};
				doc_expected_output = {};

				for i=1:number_of_tests,
					docs{i} = {};
                    S = speed_calc_obj.session;
                    
                    %taken from calculate_speed_indexes method:
                    numsteps = 5;
                    sfs = logspace(log10(0.01),log10(60),numsteps);
				    tfs = logspace(log10(0.01),log10(120),numsteps);
				    [SFs,TFs] = meshgrid(sfs,tfs);
                    function_params = speed_calc_obj.generate_mock_parameters(scope, i);
                    r_ = vlt.neuro.vision.speed.tuningfunc(SFs,TFs,function_params);
                    
                    %r = vlt.math.dog(speed_values,function_params);
                    %param_struct = struct([]); %this didn't work - need at
                    %least one parameter?
                    param_struct = struct('contrast',.5); %random choice, can be anything between 0 and 1
					independent_variable = {'temporal_frequency','spatial_frequency'};
                    x = [SFs(:),TFs(:)]; % columns
					r = r_(:); % columns
					%why do we have these? Set control (blank) stimulus to
                    %firing rate = 0?
                    %should we have just one blank stimulus row or
                    %multiple? Add nan to end of x or end of sfs and tfs?
                    x(end+1,:) = NaN;
                    r(end+1,1) = 0;
					
					switch (scope),
						case 'standard',
							reps = 5; % need reps to test significance measures
							noise = 0;
						case 'low_noise',
							reps = 10;
							noise = 0.1;
						case 'high_noise',
							reps = 10;
							noise = 1;
						otherwise,
							error(['Unknown scope ' scope '.']);
					end; % switch

					docs{i} = ndi.mock.fun.stimulus_response(S,...
						param_struct, independent_variable, x, r, noise, reps);

                    calcparameters = speed_calc_obj.default_search_for_input_parameters();
                    calcparameters.query.query = ndi.query('stimulus_tuningcurve.independent_variable_label','contains_string','temporal_frequency','');
					calcparameters.query.query = calcparameters.query.query & ...
                        ndi.query('stimulus_tuningcurve.independent_variable_label','contains_string','spatial_frequency','');
                    calcparameters.query.query = calcparameters.query.query & ...
						ndi.query('','depends_on','element_id',docs{i}{3}.id());
                    %I = speed_calc_obj.search_for_input_parameters(calcparameters);
                    doc_output{i} = speed_calc_obj.run('Replace',calcparameters);
					if numel(doc_output{i})>1,
						error(['Generated more than one output doc when one was expected.']);
                    elseif numel(doc_output{i})==0,
						error(['Generated no output docs when one was expected.']);
					end;
					doc_output{i} = doc_output{i}{1};

					if generate_expected_docs,
						speed_calc_obj.write_mock_expected_output(i,doc_output{i});
					end;

					doc_expected_output{i} = speed_calc_obj.load_mock_expected_output(i);

				end; % for
		end; % generate_mock_docs()

		function [b,errormsg] = compare_mock_docs(speed_calc_obj, expected_doc, actual_doc, scope)
			% COMPARE_MOCK_DOCS - compare an expected calculation answer with an actual answer
			%
			% [B, ERRORMSG] = COMPARE_MOCK_DOCS(CTEST_OBJ, EXPECTED_DOC, ACTUAL_DOC, SCOPE)
			%
			% Given an NDI document with the expected answer to a calculation (EXPECTED_DOC),
			% the ACTUAL_DOC computed, and the SCOPE (a string: 'standard', 'low_noise','high_noise'),
			% this function computes whether the ACTUAL_DOC is within an allowed tolerance of
			% EXPECTED_DOC.
			%
			% B is 1 if the differences in the documents are within the tolerance of the class.
			% Otherwise, B is 0.
			% If B is 0, ERRORMSG is a string that indicates where the ACTUAL_DOC is out of tolerance.
			%

				[b_,errormsg] = ndi.calc.vis.test.speed_tuning_compare_docs(expected_doc,actual_doc,scope);
        		b = ~isempty(find(b_, 1)); %b is 1 if b_ has no 0s, i.e. there are no errors

		end;

        function [P, total] = generate_mock_parameters(speed_calc_obj, scope, index)
			% generate_mock_parameters - generate mock parameters for testing ndi.calc.vis.oridir_tuning
			%
			% [P, TOTAL] = ndi.calc.vis.generate_mock_parameters(scope, index)
			%
			% Generates a parameter set for generating a mock document with a given index value.
			% P will be a row vector of parameters [a1 b1 a2 b2].
			% TOTAL is the total number of mock stimuli that are available to be generated.
			% 
			% SCOPE can be 'standard', 'random_nonoise', or 'random_noisy'.
			% INDEX selects which parameters are used to generate a mock document (from 1..TOTAL, wrapped
			% using MOD).
			% 
                %1st test taken from example in vlt.neuro.vision.speed.tuningfunc:
                A = [1,1,1]; %add more tests by adding to end of each parameter's vector
                zeta = [0,0,0];
                xi = [0,0,0];
                sigma_sf = [0.2,0.2,2]; % Cycles per degree
                sigma_tf = [4,4,2]; % Cycles per second; this is the fall off
                sf0 = [0.1,1,1];
                tf0 = [4,1,1];
                P_ = [ A(:) zeta(:) xi(:) sigma_sf(:) sigma_tf(:) sf0(:) tf0(:)] ; 

                % %1st test from spatial_frequency_tuning:
                % sf_params = [1 1 0 1]; 
                % %1st test from temporal_frequency_tuning:
                % tf_params = [1 1 0 1];
                % P_ = [sf_params tf_params];
				total = size(P_,1);

				actual_index = 1+mod(index-1,total);

				% no dependence on scope for this stimulus type

				P = P_(actual_index,:);
		end; % generate_mock_parameters
	end; % methods()
end % speed_tuning
