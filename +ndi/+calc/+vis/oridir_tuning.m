classdef oridir_tuning < ndi.calculator

	methods
		function oridir_tuning_obj = oridir_tuning(session)
			% oridir_tuning - ndi.calculator object that
			% calculates orientation and direction tuning curves from spike
			% elements
			%
			% ORIDIRTUNING_OBJ = ORIDIRTUNING(SESSION)
			%
			% Creates a oridir_tuning ndi.calculator object
			%
				w = which('ndi.calc.vis.contrast_tuning');
				parparparpar = fileparts(fileparts(fileparts(fileparts(w))));                
				oridir_tuning_obj = oridir_tuning_obj@ndi.calculator(session,'oridir_tuning',...
					fullfile(parparparpar,'ndi_common','database_documents','calc','oridirtuning_calc.json'));
				oridir_tuning_obj.numberOfSelfTests = 9;
		end; % oridir_tuning() creator

		function doc = calculate(ndi_calculator_obj, parameters)
			% CALCULATE - perform the calculator for
			% ndi.calc.oridir_tuning
			%
			% DOC = CALCULATE(NDI_CALCULATION_OBJ, PARAMETERS)
			%
			% Creates a oridir_tuning_direction_tuning_calc document given input parameters.
			
				% Step 1. Check inputs
				if ~isfield(parameters,'input_parameters'),
					error(['parameters structure lacks ''input_parameters.''']);
				end;
				if ~isfield(parameters,'depends_on'),
					error(['parameters structure lacks ''depends_on.''']);
				end;
						
				% Step 2. Set up output structure
				oridir_tuning = parameters;
						
				tuning_doc = ndi_calculator_obj.session.database_search(ndi.query('base.id',...
					'exact_string',...
					did.db.struct_name_value_search(parameters.depends_on,'stimulus_tuningcurve_id'),''));
				if numel(tuning_doc)~=1, 
					error(['Could not find stimulus tuning curve doc..']);
				end;
				tuning_doc = tuning_doc{1};
						
				% Step 3. Calculate oridir_tuning and direction indexes from
				% stimulus responses and write output into an oridir_tuning document

				oriapp = ndi.app.oridirtuning(ndi_calculator_obj.session);
				doc = ndi_calculator_obj.calculate_oridir_indexes(tuning_doc) + ...
					ndi_calculator_obj.newdocument();

				% Step 4. Check if doc exists
				if ~isempty(doc), 
					doc = ndi.document(ndi_calculator_obj.doc_document_types{1},...
						'oridirtuning_calc',oridir_tuning) + doc;
					doc = doc.set_dependency_value('stimulus_tuningcurve_id',tuning_doc.id());
				end;
		end; % calculate

		function parameters = default_search_for_input_parameters(ndi_calculator_obj, varargin)
			% DEFAULT_SEARCH_FOR_INPUT_PARAMETERS - default parameters for searching for inputs
			%
			% PARAMETERS = DEFAULT_SEARCH_FOR_INPUT_PARAMETERS(NDI_CALCULATION_OBJ)
			%
			% Returns a list of the default search parameters for finding appropriate inputs
			% to the calculator.
			%
				% search for stimulus_tuningcurve_id
				parameters.input_parameters = struct([]);
				parameters.depends_on = did.datastructures.emptystruct('name','value');
				parameters.query = ndi_calculator_obj.default_parameters_query(parameters);
		end; % default_search_for_input_parameters

		function query = default_parameters_query(ndi_calculator_obj, parameters_specification)
			% DEFAULT_PARAMETERS_QUERY - what queries should be used to
			% search for input parameters
			%
			% QUERY = DEFAULT_PARAMETERS_QUERY(NDI_CALCULATION_OBJ,
			% PARAMETERS_SPECIFICATION)
			%
			% Calling SEARCH_FOR_INPUT_PARAMETERS allows for users to
			% specify a 'query' structure to select particular documents to
			% be placed into the 'depends_on' parameter specification.
			% If a 'query' structure is not provided, the default will be
			% used.
			%
			% The function returns: 
			% |-----------|--------------------------------------------|
			% | query     | A structure with 'name' and 'query' fields |
			% |           | that describes a search to be performed to |
			% |           | identify inputs for the 'depends_on' field |
			% |           | in the PARAMETERS output.                  |
			% |-----------|--------------------------------------------|
			%
			% For the ndi.calc.vision.tuning_curve class, this looks for
			% documents of type 'stimulus_tuningcurve.json' with
			% 'response_type' fields that contain 'mean' or 'F1'.
			%
			%         
				q1 = ndi.query('','isa','stimulus_tuningcurve.json','');
			
				q2 = ndi.query('tuning_curve.independent_variable_label','exact_string_anycase','Orientation','');
				q3 = ndi.query('tuning_curve.independent_variable_label','exact_string_anycase','Direction','');
				q4 = ndi.query('tuning_curve.independent_variable_label','exact_string_anycase','angle','');
				q234 = q2 | q3 | q4;
				q_total = q1 & q234;
			
				query = struct('name','stimulus_tuningcurve_id','query',q_total);
		
		end; % default_parameters_query()
		
		function b = is_valid_dependency_input(ndi_calculator_obj, name, value)
			% IS_VALID_DEPENDENCY_INPUT - checks if a potential dependency input
			% actually valid for this calculator
			% 
			% B = IS_VALID_DEPENDENCY_INPUT(NDI_CALCULATION_OBJ, NAME,
			% VALUE)
			%
			% Tests whether a potential input to a calculator is valid.
			% NAME - potential dependency name
			% VALUE - base id of the potential dependency name
			%
			% The base class behavior of this function will return true.
			% This is overridden if additional criteria beyond an ndi.query
			% are needed to assess if a document is an appropriate input
			% for the calculator.
				b = 1;
				return;

				% could also use the below, but will require an extra query operation
				% and updating for speed

				switch lower(name),
					case lower('stimulus_tuningcurve_id'),
						q = ndi.query('base.id','exact_string',value,'');
						d = ndi_calculator_obj.S.database_search(q);
						b = (numel(d.document_properties.independent_variable_label) ==2);
					end;
		end; % is_valid_dependency_input()

		function oriprops_doc = calculate_oridir_indexes(ndi_calculator_obj, tuning_doc)
			% CALCULATE_ORIDIR_INDEXES - calculate orientation and direction index values from a tuning curve
			%
			% ORIDIR_DOC = CALCULATE_ORIDIR_INDEXES(NDI_ORIDIRTUNING_CALC_OBJ, TUNING_DOC)
			%
			% Given a 2-dimensional tuning curve document with measurements
			% at orientation and direction frequencies, this function calculates oridir_tuning
			% parameters and stores them in ORIDIRTUNING document ORIDIR_DOC.
			%
			%
				ind = {};
				ind_real = {};
				control_ind = {};
				control_ind_real = {};
				response_ind = {};
				response_mean = [];
				response_stddev = [];
				response_stderr = [];
                
				% stim_response_doc
				stim_response_doc = ndi_calculator_obj.session.database_search(ndi.query('base.id',...
					'exact_string',tuning_doc.dependency_value('stimulus_response_scalar_id'),''));
				if numel(stim_response_doc)~=1,
					error(['Could not find stimulus response scalar document.']);
				end;
				if iscell(stim_response_doc),
					stim_response_doc = stim_response_doc{1};
				end;
                
				tuning_doc = ndi.app.stimulus.tuning_response.tuningdoc_fixcellarrays_static(tuning_doc);

				for i=1:numel(tuning_doc.document_properties.stimulus_tuningcurve.individual_responses_real),
					ind{i} = tuning_doc.document_properties.stimulus_tuningcurve.individual_responses_real{i} + ...
						sqrt(-1)*tuning_doc.document_properties.stimulus_tuningcurve.individual_responses_imaginary{i};
					ind_real{i} = ind{i};
					if any(~isreal(ind_real{i})), 
						ind_real{i} = abs(ind_real{i}); 
					end;
					control_ind{i} = tuning_doc.document_properties.stimulus_tuningcurve.control_individual_responses_real{i} + ...
						sqrt(-1)*tuning_doc.document_properties.stimulus_tuningcurve.control_individual_responses_imaginary{i};
					control_ind_real{i} = control_ind{i};
					if any(~isreal(control_ind_real{i})), 
						control_ind_real{i} = abs(control_ind_real{i}); 
					end;
					response_ind{i} = ind{i} - control_ind{i};
					response_mean(i) = nanmean(response_ind{i});
					if ~isreal(response_mean(i)), 
						response_mean(i) = abs(response_mean(i)); 
					end;
					response_stddev(i) = nanstd(response_ind{i});
					response_stderr(i) = vlt.data.nanstderr(response_ind{i});
					if any(~isreal(response_ind{i})),
						response_ind{i} = abs(response_ind{i});
					end;
				end;
                   
				properties.coordinates = 'compass';
				properties.response_units = tuning_doc.document_properties.stimulus_tuningcurve.response_units;
				properties.response_type = stim_response_doc.document_properties.stimulus_response_scalar.response_type;

				response.curve = ...
					[ tuning_doc.document_properties.stimulus_tuningcurve.independent_variable_value(:)' ; ...
						response_mean ; ...
						response_stddev ; ...
						response_stderr; ];
				response.ind = response_ind;

 				vi = vis.oridir.index.oridir_vectorindexes(response);
 				fi = vis.oridir.index.oridir_fitindexes(response);
             
				resp.ind = ind_real;
				resp.blankind = control_ind_real{1};
				resp = ndi.app.stimulus.tuning_response.tuningcurvedoc2vhlabrespstruct(tuning_doc);
				[anova_across_stims, anova_across_stims_blank] = neural_response_significance(resp);

				tuning_curve = struct(...
					'direction', vlt.data.colvec(tuning_doc.document_properties.stimulus_tuningcurve.independent_variable_value), ...
					'mean', response_mean(:), ...
					'stddev', response_stddev(:), ...
					'stderr', response_stderr(:), ...
					'individual', vlt.data.cellarray2mat(response_ind), ...
					'raw_individual', vlt.data.cellarray2mat(ind_real), ...
					'control_individual', vlt.data.cellarray2mat(control_ind_real));

				significance = struct('visual_response_anova_p',anova_across_stims_blank,...
					'across_stimuli_anova_p', anova_across_stims);

				vector = struct('circular_variance', vi.ot_circularvariance, ...
					'direction_circular_variance', vi.dir_circularvariance', ...
					'Hotelling2Test', vi.ot_HotellingT2_p, ...
					'orientation_preference', vi.ot_pref, ...
					'direction_preference', vi.dir_pref, ...
					'direction_hotelling2test', vi.dir_HotellingT2_p, ...
					'dot_direction_significance', vi.dir_dotproduct_sig_p);

				fit = struct('double_gaussian_parameters', fi.fit_parameters,...
					'double_gaussian_fit_angles', vlt.data.colvec(fi.fit(1,:)), ...
					'double_gaussian_fit_values', vlt.data.colvec(fi.fit(2,:)), ...
					'orientation_preferred_orthogonal_ratio', fi.ot_index, ...
					'direction_preferred_null_ratio', fi.dir_index, ...
					'orientation_preferred_orthogonal_ratio_rectified', fi.ot_index_rectified', ...
					'direction_preferred_null_ratio_rectified', fi.dir_index_rectified, ...
					'orientation_angle_preference', mod(fi.dirpref,180), ...
					'direction_angle_preference', fi.dirpref, ...
					'hwhh', fi.tuning_width);

				% create document and store in oridir_tuning
				oriprops_doc = ndi.document('orientation_direction_tuning',...
					'orientation_direction_tuning',vlt.data.var2struct('properties', 'tuning_curve', 'significance', 'vector', 'fit'));
                                oriprops_doc = oriprops_doc.set_dependency_value('element_id', stim_response_doc.dependency_value('element_id'));
				oriprops_doc = oriprops_doc.set_dependency_value('stimulus_tuningcurve_id',tuning_doc.id());
		end; %calculate_oridir_indexes()
    
		function h=plot(ndi_calculator_obj, doc_or_parameters, varargin)
			% PLOT - provide a diagnostic plot to show the results of the calculator
			%
			% H=PLOT(NDI_CALCULATION_OBJ, DOC_OR_PARAMETERS, ...)
			%
			% Produce a plot of the tuning curve.
			%
			% Handles to the figure, the axes, and any objects created are returned in H.
			%
			% This function takes additional input arguments as name/value pairs.
			% See ndi.calculator.plot_parameters for a description of those parameters.

				plot_tuning_curve_log = 0;

				vlt.data.assign(varargin{:});

				% call superclass plot method to set up axes
				h=plot@ndi.calculator(ndi_calculator_obj, doc_or_parameters, varargin{:});
				
				% Check doc parameters
				if isa(doc_or_parameters,'ndi.document'),
					doc = doc_or_parameters;
				else,
					error(['Do not know how to proceed without an ndi document for doc_or_parameters.']);
				end;
           
				ot = doc.document_properties.orientation_direction_tuning;  % set variable for less typing
            
				% Set up plot
				ha = vlt.plot.myerrorbar(ot.tuning_curve.direction, ...
					ot.tuning_curve.mean, ...
					ot.tuning_curve.stderr, ...
					ot.tuning_curve.stderr);
                
				delete(ha(2));
				set(ha(1), 'color', [0 0 0]);
				h.objects(end+1) = ha(1);
                
				% Plot responses
				hold on;
				h_baseline = plot([0 360], [0 0], 'k--');
				h_fitline = plot(ot.fit.double_gaussian_fit_angles,...
					ot.fit.double_gaussian_fit_values,'k-');
				h.objects(end+1) = h_baseline;
				h.objects(end+1) = h_fitline;
			
				% Set labels
				if ~h.params.suppress_x_label,
					h.xlabel = xlabel('Direction (\circ)');
				end;
				if ~h.params.suppress_y_label,
					h.ylabel = ylabel(ot.properties.response_units);
				end;

				if 1, % when database is faster :-/ it is faster
					if ~h.params.suppress_title,
						element = ndi.database.fun.ndi_document2ndi_object(doc.dependency_value('element_id'),ndi_calculator_obj.session);
						title_str = [element.elementstring() '.' element.type '; ' ot.properties.response_type];

						if plot_tuning_curve_log,
							log_str = ndi.fun.calc.stimulus_tuningcurve_log(ndi_calculator_obj.session,doc);
							title_str = {title_str, log_str};
						end;
						h.title = title(title_str,'interp','none');
					end;
				end;
				box off;
		end; % plot()

		% TESTING METHODS

		function [docs, doc_output, doc_expected_output] = generate_mock_docs(oridir_calc_obj, scope, number_of_tests, kwargs)
			% GENERATE_MOCK_DOCS - generate mock documents and expected answers for tests
			%
			% [DOCS, DOC_OUTPUT, DOC_EXPECTED_OUTPUT] = GENERATE_MOCK_DOCS(ORIDIR_CALC_OBJ, ...
			%    SCOPE, NUMBER_OF_TESTS, ...)
			%
			% Creates a set of documents to test ndi.calc.vis.oridir_tuning.
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
			% | specific_test_inds([])     | Should we specify which tests to run?             |
			% |--------------------------|---------------------------------------------------|
			%
				arguments
					oridir_calc_obj
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
				%if not specifying the number of tests, just use the number
				%given by number_of_tests; otherwise use the test indices
				%specified by specific_test_inds
				if numel(specific_test_inds) == 0
					specific_test_inds = 1:number_of_tests;
				end

				for i=specific_test_inds,
					docs{end+1} = {};

					parameters = oridir_calc_obj.generate_mock_parameters(scope, i);

					angles = 0:30:360-30; % use these angles
					r = vis.oridir.doublegaussianfunc(angles,parameters);

					param_struct = struct('sFrequency',0.5,'tFrequency',2);
					independent_variable = {'angle'};
					x = angles(:); % column
					r = r(:); % column
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
					docs{end} = ndi.mock.fun.stimulus_response(oridir_calc_obj.session,...
						param_struct, independent_variable, x, r, noise, reps);

					calcparameters = oridir_calc_obj.default_search_for_input_parameters();
					calcparameters.query.query = ndi.query('stimulus_tuningcurve.independent_variable_label','contains_string','angle','');
					calcparameters.query.query = calcparameters.query.query & ...
						ndi.query('','depends_on','element_id',docs{end}{3}.id());
					doc_output{end+1} = oridir_calc_obj.run('Replace',calcparameters);
					if numel(doc_output{end})>1,
						error(['Generated more than one output doc when one was expected.']);
					elseif numel(doc_output{end})==0,
						error(['Generated no output docs when one was expected.']);
					end;
					doc_output{end} = doc_output{end}{1}; %what's the point of this?

					if generate_expected_docs,
						oridir_calc_obj.write_mock_expected_output(i,doc_output{end});
					end;

					doc_expected_output{end+1} = oridir_calc_obj.load_mock_expected_output(i);

				end; % for
		end; % generate_mock_docs()

		function [b,errormsg] = compare_mock_docs(oridir_calc_obj, expected_doc, actual_doc, scope)
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
			% Uses the function ndi.calc.vis.test.oridir_compare_docs().
			%
				[b,errormsg] = ndi.calc.vis.test.oridir_compare_docs(expected_doc,actual_doc,scope);
				errormsg = cat(2,errormsg{:});
				b = all(b);
		end;

		function [P, total] = generate_mock_parameters(oridir_calc_obj, scope, index)
			% generate_mock_parameters - generate mock parameters for testing ndi.calc.vis.oridir_tuning
			%
			% [P, TOTAL] = ndi.calc.vis.generate_mock_parameters(scope, index)
			%
			% Generates a parameter set for generating a mock document with a given index value.
			% P will be a row vector of parameters [Rsp Rp Rn theta sigma].
			% TOTAL is the total number of mock stimuli that are available to be generated.
			% 
			% SCOPE can be 'standard', 'random_nonoise', or 'random_noisy'.
			% INDEX selects which parameters are used to generate a mock document (from 1..TOTAL, wrapped
			% using MOD).
			% 

				P_(1,:) = [ 0 20 10 45 30] ; % response of 20 in preferred direction of 45 degrees, 10 opposite
				P_(2,:) = [ 0 20 10 45 45] ; % broader tuning
				P_(3,:) = [ 0 20 10 45 90] ; % really broad tuning 
				P_(4,:) = [ 0 20 10 45 90] ; % really broad tuning 
				P_(5,:) = [ 10 20 10 45 30] ; % large offset
				P_(6,:) = [ 10 20 19 45 30] ; % really low direction index offset
				P_(7,:) = [0 20 10 20 45] ; % Narrower tuning
				P_(8,:) = [0 20 10 10 45] ; %Extremely narrow tuning
				P_(9,:) = [0 20 20 45 45] ; %Equal response Rp and Rn

					% we should add more

				total = size(P_,1);

				actual_index = 1+mod(index-1,total);

				% no dependence on scope for this stimulus type

				P = P_(actual_index,:);

		end; % generate_mock_parameters

	end; % methods()
			
end % oridir_tuning
