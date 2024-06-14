classdef temporal_frequency_tuning < ndi.calculator

	methods
		function temporal_frequency_tuning_obj = temporal_frequency_tuning(session)
			% TEMPORAL_FREQUENCY_TUNING - a temporal_frequency_tuning demonstration of an ndi.calculator object
			%
			% TEMPORAL_FREQUENCY_TUNING_OBJ = TEMPORAL_FREQUENCY_TUNING(SESSION)
			%
			% Creates a TEMPORAL_FREQUENCY_TUNING ndi.calculator object
			%
				ndi.globals;
				w = which('ndi.calc.vis.temporal_frequency_tuning');
				parparparpar = fileparts(fileparts(fileparts(fileparts(w))));
				temporal_frequency_tuning_obj = temporal_frequency_tuning_obj@ndi.calculator(session,'temporal_frequency_tuning_calc',...
					fullfile(parparparpar,'ndi_common','database_documents','calc','temporal_frequency_tuning_calc.json'));
		end; % temporal_frequency_tuning()

		function doc = calculate(ndi_calculator_obj, parameters)
			% CALCULATE - perform the calculator for ndi.calc.example.temporal_frequency_tuning
			%
			% DOC = CALCULATE(NDI_CALCULATOR_OBJ, PARAMETERS)
			%
			% Creates a temporal_frequency_tuning_calc document given input parameters.
			%
			% The document that is created temporal_frequency_tuning
			% by the input parameters.
				% check inputs
				if ~isfield(parameters,'input_parameters'), error(['parameters structure lacks ''input_parameters''.']); end;
				if ~isfield(parameters,'depends_on'), error(['parameters structure lacks ''depends_on''.']); end;
				
				% Step 1: set up the output structure
				temporal_frequency_tuning_calc = parameters;

				tuning_response_doc = ndi_calculator_obj.session.database_search(ndi.query('base.id','exact_string',...
					did.db.struct_name_value_search(parameters.depends_on,'stimulus_tuningcurve_id'),''));
				if numel(tuning_response_doc)~=1, 
					error(['Could not find stimulus tuning doc..']);
				end;
				tuning_response_doc = tuning_response_doc{1};

				% Step 2: perform the calculator, which here creates a temporal_frequency_tuning doc
				doc = ndi_calculator_obj.calculate_temporal_frequency_indexes(tuning_response_doc) + ...
					ndi_calculator_obj.newdocument();
				
				if ~isempty(doc), 
					doc = ndi.document(ndi_calculator_obj.doc_document_types{1},'temporal_frequency_tuning_calc',temporal_frequency_tuning_calc) + doc;
				end;
		end; % calculate

		function parameters = default_search_for_input_parameters(ndi_calculator_obj)
			% DEFAULT_SEARCH_FOR_INPUT_PARAMETERS - default parameters for searching for inputs
			%
			% PARAMETERS = DEFAULT_SEARCH_FOR_INPUT_PARAMETERS(NDI_CALCULATOR_OBJ)
			%
			% Returns a list of the default search parameters for finding appropriate inputs
			% to the calculator. For temporal_frequency_tuning_calc, there is no appropriate default parameters
			% so this search will yield empty.
			%
				parameters.input_parameters = struct([]);
				parameters.depends_on = did.datastructures.emptystruct('name','value');
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
			% For the ndi.calc.stimulus.temporal_frequency_tuning_calc class, this looks for 
			% documents of type 'stimulus_response_scalar' with 'response_type' fields
			% the contain 'mean' or 'F1'.
			%
			%
				q1 = ndi.query('','isa','stimulus_tuningcurve','');
				q2 = ndi.query('stimulus_tuningcurve.independent_variable_label','contains_string','Temporal_Frequency','');
				q3 = ndi.query('stimulus_tuningcurve.independent_variable_label','contains_string','temporal_frequency','');
				q4 = ndi.query('stimulus_tuningcurve.independent_variable_label','contains_string','TEMPORAL_FREQUENCY','');
				q22 = ndi.query('stimulus_tuningcurve.independent_variable_label','contains_string','Temporal Frequency','');
				q32 = ndi.query('stimulus_tuningcurve.independent_variable_label','contains_string','temporal frequency','');
				q42 = ndi.query('stimulus_tuningcurve.independent_variable_label','contains_string','TEMPORAL FREQUENCY','');
				q234 = q2 | q3 | q4 | q22 | q32 | q42;
				q_total = q1 & q234;

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
	
				switch lower(name),
					case lower('stimulus_tuningcurve_id'),
						q = ndi.query('base.id','exact_string',value,'');
						d = ndi_calculator_obj.S.database_search(q);
						b = (numel(d.document_properties.independent_variable_label) ==1);
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
                        % See ndi.calculator.plot_parameters for a description of many parameters.
			% Also takes:
			% |----------------------|-------------------------------------------------------|
			% | Parameter (default)  | Description                                           |
			% |----------------------|-------------------------------------------------------|
			% | useAbsolute (0)      | Plot the absolute value of the responses and fits of  |
			% |                      |   absolute value.                                     |
			% |----------------------|-------------------------------------------------------|
			% 

				x_axis = [0.01 120];
				useAbsolute = 0;
				did.datastructures.assign(varargin{:});

				% call superclass plot method to set up axes
				h=plot@ndi.calculator(ndi_calculator_obj, doc_or_parameters, varargin{:});

				if isa(doc_or_parameters,'ndi.document'),
					doc = doc_or_parameters;
				else,
					error(['Do not know how to proceed without an ndi document for doc_or_parameters.']);
				end;

				tft = doc.document_properties.temporal_frequency_tuning; % shorten our typing
				tc = tft.tuning_curve; % shorten our typing

				% First plot responses
				hold on;
				h_baseline = plot([min(tft.fit_spline.values) max(tft.fit_spline.values)],...
					[0 0],'k--','linewidth',1.0001);
				h_baseline.Annotation.LegendInformation.IconDisplayStyle = 'off';
				h.objects(end+1) = h_baseline;
				[v,sortorder] = sort(tc.temporal_frequency);

				myfun = @(x) x;
				if useAbsolute,
					myfun = @(x) abs(x);
				end

				h_errorbar = errorbar(tc.temporal_frequency(sortorder(:)),...
					myfun(tc.mean(sortorder(:))),...
					tc.stderr(sortorder(:)),tc.stderr(sortorder(:)));
				set(h_errorbar,'color',[0 0 0],'linewidth',1,'linestyle','none');
				set(gca,'xscale','log');
				h.objects = cat(2,h.objects,h_errorbar);
				
				% Second plot all fits

				linestyle = '--';
				if tft.significance.visual_response_anova_p<0.05,
					linestyle = '-';
				end;
				tft_o = tft;

				if useAbsolute,
					tft = doc.document_properties.temporal_frequency_tuning.abs;
				end

					% drop spline, gausslog because not good
				%h_fit = plot(tft.fit_spline.values,tft.fit_spline.fit,['k' linestyle] );
				%h.objects = cat(2,h.objects,h_fit);
				h_fit = plot(tft.fit_dog.values,tft.fit_dog.fit,['m' linestyle]);
				h.objects = cat(2,h.objects,h_fit);
				h_fit = plot(tft.fit_movshon.values,tft.fit_movshon.fit,['b' linestyle],'linewidth',2);
				h.objects = cat(2,h.objects,h_fit);
				h_fit = plot(tft.fit_movshon_c.values,tft.fit_movshon_c.fit,['g' linestyle],'linewidth',2);
				h.objects = cat(2,h.objects,h_fit);
				%h_fit = plot(tft.fit_gausslog.values,tft.fit_gausslog.fit,['g' linestyle]);
				%h.objects = cat(2,h.objects,h_fit);

				if ~h.params.suppress_x_label,
					h.xlabel = xlabel('Temporal frequency');
				end;
				if ~h.params.suppress_y_label,
					h.ylabel = ylabel(['Response (' tft_o.properties.response_type ', ' tft_o.properties.response_units ')']);
				end;

				set(gca,'xlim',x_axis);

				if 0, % when database is faster :-/
					if ~h.params.suppress_title,
						element = ndi.database.fun.ndi_document2ndi_object(doc.dependency_value('element_id'),ndi_calculator_obj.session);
						h.title = title(element.elementstring(), 'interp','none');
					end;
				end;
				box off;

		end; % plot()

		function temporal_frequency_props_doc = calculate_temporal_frequency_indexes(ndi_calculator_obj, tuning_doc)
			% CALCULATE_TEMPORAL_FREQUENCY_INDEXES - calculate contrast index values from a tuning curve
			%
			% TEMPORAL_FREQUENCY_PROPS_DOC = CALCULATE_TEMPORAL_FREQUENCY_INDEXES(NDI_TEMPORAL_FREQUENCY_TUNING_CALC_OBJ, TUNING_DOC)
			%
			% Given a 1-dimensional tuning curve document, this function calculates contrast response
			% parameters and stores them in TEMPORAL_FREQUENCY_TUNING document TEMPORAL_FREQUENCY_PROPS_DOC.
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

				resp = ndi.app.stimulus.tuning_response.tuningcurvedoc2vhlabrespstruct(tuning_doc);

				stc = tuning_doc.document_properties.stimulus_tuningcurve;

				tuning_curve = struct(...
					'temporal_frequency', ...
						vlt.data.rowvec(tuning_doc.document_properties.stimulus_tuningcurve.independent_variable_value), ...
					'mean', resp.curve(2,:), ...
					'stddev', resp.curve(3,:), ...
					'stderr', resp.curve(4,:), ...
					'individual', {resp.ind}, ...
					'control_mean', stc.control_response_mean,...
					'control_stddev', stc.control_response_stddev,...
					'control_stderr', stc.control_response_stderr,...
					'control_mean_stddev', resp.blankresp(2),...
					'control_mean_stderr', resp.blankresp(3));

				[anova_across_stims, anova_across_stims_blank] = neural_response_significance(resp);

				significance = struct('visual_response_anova_p',anova_across_stims_blank,...
					'across_stimuli_anova_p', anova_across_stims);

				tf_props = vis.temporal_frequency_analysis(resp);

				temporal_frequency_tuning.properties = properties;
				temporal_frequency_tuning.tuning_curve = tuning_curve;
				temporal_frequency_tuning.significance = significance;
				temporal_frequency_tuning.fitless = tf_props.fitless;
				temporal_frequency_tuning.fit_dog = tf_props.fit_dog;
				temporal_frequency_tuning.fit_movshon = tf_props.fit_movshon;
				temporal_frequency_tuning.fit_movshon_c = tf_props.fit_movshon_c;
				temporal_frequency_tuning.fit_spline = tf_props.fit_spline;
				temporal_frequency_tuning.fit_gausslog = tf_props.fit_gausslog;

				resp_abs = resp;
				resp_abs.curve(2,:) = abs(resp_abs.curve(2,:));
				abs_tf_props = vis.temporal_frequency_analysis(resp_abs);

				temporal_frequency_tuning.abs.fitless = abs_tf_props.fitless;
				temporal_frequency_tuning.abs.fit_dog = abs_tf_props.fit_dog;
				temporal_frequency_tuning.abs.fit_movshon = abs_tf_props.fit_movshon;
				temporal_frequency_tuning.abs.fit_movshon_c = abs_tf_props.fit_movshon_c;
				temporal_frequency_tuning.abs.fit_spline = abs_tf_props.fit_spline;
				temporal_frequency_tuning.abs.fit_gausslog = abs_tf_props.fit_gausslog;

				temporal_frequency_props_doc = ndi.document('temporal_frequency_tuning',...
					'temporal_frequency_tuning',temporal_frequency_tuning);
				temporal_frequency_props_doc = temporal_frequency_props_doc.set_dependency_value('element_id', ...
					tuning_doc.dependency_value('element_id'));
				temporal_frequency_props_doc = temporal_frequency_props_doc.set_dependency_value('stimulus_tuningcurve_id',tuning_doc.id());

		end; % calculate_temporal_frequency_indexes()
        
        % TESTING METHODS

        function [docs, doc_output, doc_expected_output] = generate_mock_docs(temporal_freq_calc_obj, scope, number_of_tests, varargin)
			% GENERATE_MOCK_DOCS - generate mock documents and expected answers for tests
			%
			% [DOCS, DOC_OUTPUT, DOC_EXPECTED_OUTPUT] = GENERATE_MOCK_DOCS(TEMPORAL_FREQ_CALC_OBJ, ...
			%    SCOPE, NUMBER_OF_TESTS, ...)
			%
			% Creates a set of documents to test ndi.calc.vis.temporal_frequency_tuning.
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
                    S = temporal_freq_calc_obj.session;
                    function_params = temporal_freq_calc_obj.generate_mock_parameters(scope, i);
                    numsteps = 100; %sets size of x
					%spatial_freq_values = [.05, .1, .15, .2, .3, .5, .8]; %spatial frequency values commonly used in experiments, in units of cpd (cycles per degree)
                    temporal_freq_values = logspace(-2,log10(60),numsteps);
                    r = vlt.math.dog(temporal_freq_values,function_params);
                    %param_struct = struct([]); %this didn't work - need at
                    %least one parameter?
                    param_struct = struct('spatial_frequency',.5);
					independent_variable = {'temporal_frequency'}; % is the underscore required?
					x = temporal_freq_values(:); % column
					r = r(:); % column
					%why do we have these?
                    x(numsteps+1,1) = NaN;
					r(numsteps+1,1) = 0;
					
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

                    calcparameters = temporal_freq_calc_obj.default_search_for_input_parameters();
                    calcparameters.query.query = ndi.query('stimulus_tuningcurve.independent_variable_label','contains_string','temporal_frequency','');
					calcparameters.query.query = calcparameters.query.query & ...
						ndi.query('','depends_on','element_id',docs{i}{3}.id());
                    I = temporal_freq_calc_obj.search_for_input_parameters(calcparameters);
                    doc_output{i} = temporal_freq_calc_obj.run('Replace',calcparameters);
					if numel(doc_output{i})>1,
						error(['Generated more than one output doc when one was expected.']);
                    elseif numel(doc_output{i})==0,
						error(['Generated no output docs when one was expected.']);
					end;
					doc_output{i} = doc_output{i}{1};

					if generate_expected_docs,
						temporal_freq_calc_obj.write_mock_expected_output(i,doc_output{i});
					end;

					doc_expected_output{i} = temporal_freq_calc_obj.load_mock_expected_output(i);

				end; % for
		end; % generate_mock_docs()

		function [b,errormsg] = compare_mock_docs(temporal_freq_calc_obj, expected_doc, actual_doc, scope)
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

				[b_,errormsg] = ndi.calc.vis.test.temporal_frequency_tuning_compare_docs(expected_doc,actual_doc,scope);	%need to implement
        		b = ~isempty(find(b_, 1)); %b is 1 if b_ has no 0s, i.e. there are no errors

		end;

        function [P, total] = generate_mock_parameters(temporal_freq_calc_obj, scope, index)
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
				P_(1,:) = [ 1 1 0 1 ] ; %regular gaussian with peak 1 and width parameter set to 1
				total = size(P_,1);

				actual_index = 1+mod(index-1,total);

				% no dependence on scope for this stimulus type

				P = P_(actual_index,:);
		end; % generate_mock_parameters
	end; % methods()
end % temporal_frequency_tuning
