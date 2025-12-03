classdef spatial_frequency_tuning < ndi.calculator

	methods
		function spatial_frequency_tuning_obj = spatial_frequency_tuning(session)
			% SPATIAL_FREQUENCY_TUNING - a spatial_frequency_tuning demonstration of an ndi.calculator object
			%
			% SPATIAL_FREQUENCY_TUNING_OBJ = SPATIAL_FREQUENCY_TUNING(SESSION)
			%
			% Creates a SPATIAL_FREQUENCY_TUNING ndi.calculator object
			%
				w = which('ndi.calc.vis.spatial_frequency_tuning');
				parparparpar = fileparts(fileparts(fileparts(fileparts(w))));
				spatial_frequency_tuning_obj = spatial_frequency_tuning_obj@ndi.calculator(session,'spatial_frequency_tuning_calc',...
					fullfile(parparparpar,'ndi_common','database_documents','calc','spatial_frequency_tuning_calc.json'));
				spatial_frequency_tuning_obj.numberOfSelfTests = 23;
		end; % spatial_frequency_tuning()

		function doc = calculate(ndi_calculator_obj, parameters)
			% CALCULATE - perform the calculator for ndi.calc.example.spatial_frequency_tuning
			%
			% DOC = CALCULATE(NDI_CALCULATOR_OBJ, PARAMETERS)
			%
			% Creates a spatial_frequency_tuning_calc document given input parameters.
			%
			% The document that is created spatial_frequency_tuning
			% by the input parameters.
				% check inputs
				if ~isfield(parameters,'input_parameters'), error(['parameters structure lacks ''input_parameters''.']); end;
				if ~isfield(parameters,'depends_on'), error(['parameters structure lacks ''depends_on''.']); end;
				
				% Step 1: set up the output structure
				spatial_frequency_tuning_calc = parameters;

				tuning_response_doc = ndi_calculator_obj.session.database_search(ndi.query('base.id','exact_string',...
					did.db.struct_name_value_search(parameters.depends_on,'stimulus_tuningcurve_id'),''));
				if numel(tuning_response_doc)~=1, 
					error(['Could not find stimulus tuning doc..']);
				end;
				tuning_response_doc = tuning_response_doc{1};

				% Step 2: perform the calculator, which here creates a spatial_frequency_tuning doc
				doc = ndi_calculator_obj.calculate_spatial_frequency_indexes(tuning_response_doc) + ...
					ndi_calculator_obj.newdocument();
				
				if ~isempty(doc), 
					doc = ndi.document(ndi_calculator_obj.doc_document_types{1},'spatial_frequency_tuning_calc',spatial_frequency_tuning_calc) + doc;
				end;
		end; % calculate

		function parameters = default_search_for_input_parameters(ndi_calculator_obj)
			% DEFAULT_SEARCH_FOR_INPUT_PARAMETERS - default parameters for searching for inputs
			%
			% PARAMETERS = DEFAULT_SEARCH_FOR_INPUT_PARAMETERS(NDI_CALCULATOR_OBJ)
			%
			% Returns a list of the default search parameters for finding appropriate inputs
			% to the calculator. For spatial_frequency_tuning_calc, there is no appropriate default parameters
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
			% For the ndi.calc.stimulus.spatial_frequency_tuning_calc class, this looks for 
			% documents of type 'stimulus_response_scalar' with 'response_type' fields
			% the contain 'mean' or 'F1'.
			%
			%
				q1 = ndi.query('','isa','stimulus_tuningcurve','');
				q2 = ndi.query('stimulus_tuningcurve.independent_variable_label','contains_string','Spatial_Frequency','');
				q3 = ndi.query('stimulus_tuningcurve.independent_variable_label','contains_string','spatial_frequency','');
				q4 = ndi.query('stimulus_tuningcurve.independent_variable_label','contains_string','SPATIAL_FREQUENCY','');
				q22 = ndi.query('stimulus_tuningcurve.independent_variable_label','contains_string','Spatial Frequency','');
				q32 = ndi.query('stimulus_tuningcurve.independent_variable_label','contains_string','spatial frequency','');
				q42 = ndi.query('stimulus_tuningcurve.independent_variable_label','contains_string','SPATIAL FREQUENCY','');
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

				useAbsolute = 0;
				vlt.data.assign(varargin{:});
				% call superclass plot method to set up axes
				h=plot@ndi.calculator(ndi_calculator_obj, doc_or_parameters, varargin{:});

				if isa(doc_or_parameters,'ndi.document'),
					doc = doc_or_parameters;
				else,
					error(['Do not know how to proceed without an ndi document for doc_or_parameters.']);
				end;

				sft = doc.document_properties.spatial_frequency_tuning; % shorten our typing
				tc = sft.tuning_curve; % shorten our typing

				% First plot responses
				hold on;
				h_baseline = plot([min(sft.fit_spline.values) max(sft.fit_spline.values)],...
					[0 0],'k--','linewidth',1.0001);
				h_baseline.Annotation.LegendInformation.IconDisplayStyle = 'off';
				h.objects(end+1) = h_baseline;
				[v,sortorder] = sort(tc.spatial_frequency);

				myfun = @(x) x;
				if useAbsolute,
					myfun = @(x) abs(x);
				end

				h_errorbar = errorbar(tc.spatial_frequency(sortorder(:)),...
					myfun(tc.mean(sortorder(:))),...
					tc.stderr(sortorder(:)),tc.stderr(sortorder(:)));
				set(h_errorbar,'color',[0 0 0],'linewidth',1,'linestyle','none');
				set(gca,'xscale','log');
				h.objects = cat(2,h.objects,h_errorbar);
				
				% Second plot all fits

				linestyle = '--';
				if sft.significance.visual_response_anova_p<0.05,
					linestyle = '-';
				end;
				sft_o = sft;

				if useAbsolute,
					sft = doc.document_properties.spatial_frequency_tuning.abs;
				end

					% the spline fits are terrible
				%h_fit = plot(sft.fit_spline.values,sft.fit_spline.fit,['k' linestyle] );
				%h.objects = cat(2,h.objects,h_fit);
				h_fit = plot(sft.fit_dog.values,sft.fit_dog.fit,['m' linestyle]);
				h.objects = cat(2,h.objects,h_fit);
				h_fit = plot(sft.fit_movshon.values,sft.fit_movshon.fit,['b' linestyle],'linewidth',2);
				h.objects = cat(2,h.objects,h_fit);
				h_fit = plot(sft.fit_movshon_c.values,sft.fit_movshon_c.fit,['g' linestyle],'linewidth',1.5);
				h.objects = cat(2,h.objects,h_fit);
					% the gauss log fits are terrible
				%h_fit = plot(sft.fit_gausslog.values,sft.fit_gausslog.fit,['g' linestyle]);
				%h.objects = cat(2,h.objects,h_fit);

				if ~h.params.suppress_x_label,
					h.xlabel = xlabel('Spatial frequency');
				end;
				if ~h.params.suppress_y_label,
					h.ylabel = ylabel(['Response (' sft_o.properties.response_type ', ' sft_o.properties.response_units ')']);
				end;

				if 0, % when database is faster :-/
					if ~h.params.suppress_title,
						element = ndi.database.fun.ndi_document2ndi_object(doc.dependency_value('element_id'),ndi_calculator_obj.session);
						h.title = title(element.elementstring(), 'interp','none');
					end;
				end;
				box off;

		end; % plot()

		function spatial_frequency_props_doc = calculate_spatial_frequency_indexes(ndi_calculator_obj, tuning_doc)
			% CALCULATE_SPATIAL_FREQUENCY_INDEXES - calculate contrast index values from a tuning curve
			%
			% SPATIAL_FREQUENCY_PROPS_DOC = CALCULATE_SPATIAL_FREQUENCY_INDEXES(NDI_SPATIAL_FREQUENCY_TUNING_CALC_OBJ, TUNING_DOC)
			%
			% Given a 1-dimensional tuning curve document, this function calculates contrast response
			% parameters and stores them in SPATIAL_FREQUENCY_TUNING document SPATIAL_FREQUENCY_PROPS_DOC.
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
					'spatial_frequency', ...
						vlt.data.colvec(tuning_doc.document_properties.stimulus_tuningcurve.independent_variable_value), ...
					'mean', vlt.data.colvec(resp.curve(2,:)), ...
					'stddev', vlt.data.colvec(resp.curve(3,:)), ...
					'stderr', vlt.data.colvec(resp.curve(4,:)), ...
					'individual', vlt.data.cellarray2mat(resp.ind), ...
					'control_mean', stc.control_response_mean,...
					'control_stddev', stc.control_response_stddev,...
					'control_stderr', stc.control_response_stderr,...
					'control_mean_stddev', resp.blankresp(2),...
					'control_mean_stderr', resp.blankresp(3));

				[anova_across_stims, anova_across_stims_blank] = neural_response_significance(resp);

				significance = struct('visual_response_anova_p',anova_across_stims_blank,...
					'across_stimuli_anova_p', anova_across_stims);

				sf_props = vis.frequency.spatial_frequency_analysis(resp);

				spatial_frequency_tuning.properties = properties;
				spatial_frequency_tuning.tuning_curve = tuning_curve;
				spatial_frequency_tuning.significance = significance;
				spatial_frequency_tuning.fitless = sf_props.fitless;
				spatial_frequency_tuning.fit_dog = sf_props.fit_dog;
				spatial_frequency_tuning.fit_movshon = sf_props.fit_movshon;
				spatial_frequency_tuning.fit_movshon_c = sf_props.fit_movshon_c;
				spatial_frequency_tuning.fit_spline = sf_props.fit_spline;
				spatial_frequency_tuning.fit_gausslog = sf_props.fit_gausslog;

				resp_abs = resp;
				resp_abs.curve(2,:) = abs(resp_abs.curve(2,:));
				abs_sf_props = vis.frequency.spatial_frequency_analysis(resp_abs);

				spatial_frequency_tuning.abs.fitless = abs_sf_props.fitless;
				spatial_frequency_tuning.abs.fit_dog = abs_sf_props.fit_dog;
				spatial_frequency_tuning.abs.fit_movshon = abs_sf_props.fit_movshon;
				spatial_frequency_tuning.abs.fit_movshon_c = abs_sf_props.fit_movshon_c;
				spatial_frequency_tuning.abs.fit_spline = abs_sf_props.fit_spline;
				spatial_frequency_tuning.abs.fit_gausslog = abs_sf_props.fit_gausslog;

				spatial_frequency_props_doc = ndi.document('spatial_frequency_tuning',...
					'spatial_frequency_tuning',spatial_frequency_tuning);
				spatial_frequency_props_doc = spatial_frequency_props_doc.set_dependency_value('element_id', ...
					tuning_doc.dependency_value('element_id'));
				spatial_frequency_props_doc = spatial_frequency_props_doc.set_dependency_value('stimulus_tuningcurve_id',tuning_doc.id());

		end; % calculate_spatial_frequency_indexes()
        % TESTING METHODS

        function [docs, doc_output, doc_expected_output] = generate_mock_docs(spatial_freq_calc_obj, scope, number_of_tests, kwargs)
			% GENERATE_MOCK_DOCS - generate mock documents and expected answers for tests
			%
			% [DOCS, DOC_OUTPUT, DOC_EXPECTED_OUTPUT] = GENERATE_MOCK_DOCS(SPATIAL_FREQ_CALC_OBJ, ...
			%    SCOPE, NUMBER_OF_TESTS, ...)
			%
			% Creates a set of documents to test ndi.calc.vis.spatial_frequency_tuning.
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

				arguments
					spatial_freq_calc_obj
					scope
					number_of_tests
					kwargs.generate_expected_docs (1,1) logical = false
					kwargs.specific_test_inds double = []
				end
				specific_test_inds = kwargs.specific_test_inds;
				generate_expected_docs = kwargs.generate_expected_docs;

				docs = {};
				doc_output = {};
				doc_expected_output = {};

				if numel(specific_test_inds) == 0
					specific_test_inds = 1:number_of_tests;
				end

				for i=specific_test_inds,
					if i > numel(docs)
						docs{i} = {};
					else
						docs{i} = {};
					end
                    S = spatial_freq_calc_obj.session;
                    [function_params, function_choice] = spatial_freq_calc_obj.generate_mock_parameters(scope, i);
                    numsteps = 100; %sets size of x
					%spatial_freq_values = [.05, .1, .15, .2, .3, .5, .8]; %spatial frequency values commonly used in experiments, in units of cpd (cycles per degree)
                    spatial_freq_values = logspace(-2,log10(60),numsteps);
                    switch (function_choice)
                        case 'dog'
                            r = vlt.math.dog(spatial_freq_values,function_params);
                        case 'movshon' %function_params should have 4 values 
                            r = vis.frequency.movshon2005_func(spatial_freq_values,function_params);
                        case 'movshon_c' %function_params should have 5 params 
                            r = vis.frequency.movshon2005_func(spatial_freq_values,function_params);
                        otherwise
                            error(['Unkown function ' function_choice '.']);
                    end
                    %param_struct = struct([]); %this didn't work - need at
                    %least one parameter?
                    param_struct = struct('temporal_frequency',2);
					independent_variable = {'spatial_frequency'}; % is the underscore required?
					x = spatial_freq_values(:); % column
					r = r(:); % column
					%why do we have these?
                    x(end+1,1) = NaN;
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

                    calcparameters = spatial_freq_calc_obj.default_search_for_input_parameters();
                    calcparameters.query.query = ndi.query('stimulus_tuningcurve.independent_variable_label','contains_string','spatial_frequency','');
					calcparameters.query.query = calcparameters.query.query & ...
						ndi.query('','depends_on','element_id',docs{i}{3}.id());
                    %I = spatial_freq_calc_obj.search_for_input_parameters(calcparameters);
                    %%doesn't seem to be used
                    doc_output{i} = spatial_freq_calc_obj.run('Replace',calcparameters);
					if numel(doc_output{i})>1,
						error(['Generated more than one output doc when one was expected.']);
                    elseif numel(doc_output{i})==0,
						error(['Generated no output docs when one was expected.']);
					end;
					doc_output{i} = doc_output{i}{1};

					if generate_expected_docs,
						spatial_freq_calc_obj.write_mock_expected_output(i,doc_output{i});
					end;

					doc_expected_output{i} = spatial_freq_calc_obj.load_mock_expected_output(i);

				end; % for
		end; % generate_mock_docs()

		function [b,errormsg] = compare_mock_docs(spatial_freq_calc_obj, expected_doc, actual_doc, scope)
			% COMPARE_MOCK_DOCS - compare an expected calculation answer with an actual answer
			%
			% [B, ERRORMSG] = COMPARE_MOCK_DOCS(TEMPORAL_FREQ_CALC_OBJ, EXPECTED_DOC, ACTUAL_DOC, SCOPE)
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

				[b_,errormsg] = ndi.calc.vis.test.spatial_frequency_tuning_compare_docs(expected_doc,actual_doc,scope);	
        		b = ~isempty(find(b_, 1)); %b is 1 if b_ has no 0s, i.e. there are no errors

		end;

        function [P, function_choice, total] = generate_mock_parameters(spatial_freq_calc_obj, scope, index)
			% generate_mock_parameters - generate mock parameters for testing ndi.calc.vis.spatial_frequency_tuning
			%
			% [P, TOTAL] = ndi.calc.vis.generate_mock_parameters(scope, index)
			%
			% Generates a parameter set for generating a mock document with a given index value.
			% P will be a row vector of parameters 
			% TOTAL is the total number of mock stimuli that are available to be generated.
			% 
			% SCOPE can be 'standard', 'random_nonoise', or 'random_noisy'.
			% INDEX selects which parameters are used to generate a mock document (from 1..TOTAL, wrapped
			% using MOD).
			% 
                %plot is using dog, movshon, and movshon_c fits - let's use
                %those here
                function_choice_{1} = 'dog';
				P_{1} = [ 1 1 0 1 ] ; %regular gaussian with peak 1 and width parameter set to 1
                function_choice_{2} = 'movshon';
                P_{2} = [ 10 5 12 1]; %example from vis.frequency.movshon2005_func
                function_choice_{3} = 'movshon_c';
                P_{3} = [ 10 5 12 1 -1]; 
                %tests 4 through 10 made using difference of gaussians
                for i = 4:10 
                    function_choice_{i} = 'dog';
                end
                P_{4} = [1 1 .5 1];
                P_{5} = [1 1 .5 2 ];
                P_{6} = [1 1 2 .5] ;
                P_{7} = [1 1 1 .5];
                P_{8} = [2 1 2 .5];
                P_{9} = [1 2 1 1];
                P_{10} = [1 2 -1 1];
                for i = 11:23
                    function_choice_{i} = 'movshon_c';
                end
                P_{11} = [10 10 12 1 -1]; %double the characteristic spatial frequency
                P_{12} = [5 10 12 1 -1]; %flip the scaling factor and characteristic spatial frequency
                P_{13} = [10 10 18 1 -1]; %increase the corner spatial frequency on the lower frequency side
                P_{14} = [10 10 6 1 -1]; %halve the corner spatial frequency on the lower frequency side
                P_{15} = [10 10 6 2 -1]; %double the slope of the lower frequency side
                P_{16} = [10 10 6 5 -1]; %further increase the slope of the lower frequency side
                P_{17} = [10 10 6 .5 -1]; %halve the slope of the lower frequency side
                P_{18} = [10 10 6 .5 1]; %raise the constant to above baseline
                P_{19} = [10 10 6 .5 0]; %keep the constant at 0
                P_{20} = [20 10 6 .5 -.5]; %double the scaling factor
                P_{21} = [20 10 0 .5 -.5]; %set the corner frequency of the lower frequency side to lowest possible (0)
                P_{22} = [20 10 6 0 -5]; %another boundary condition: set slope of lower frequency side to 0
                P_{23} = [0 10 6 .5 -.5]; %set scaling factor to 0 as an example of no tuning
                total = size(P_,2); % P_ is a 1xN cell array, this sets total equal to N 
				actual_index = 1+mod(index-1,total);

				% no dependence on scope for this stimulus type

				P = P_{actual_index};
                function_choice = function_choice_{actual_index};
		end; % generate_mock_parameters
	end; % methods()
end % spatial_frequency_tuning
