classdef oridir_bayes < ndi.calculator

	methods
		function oridir_bayes_obj = oridir_bayes(session)
			% oridir_bayes - ndi.calculator object that
			% calculates orientation and direction tuning curves from spike
			% elements
			%
			% ORIDIRTUNING_OBJ = ORIDIRTUNING(SESSION)
			%
			% Creates a oridir_bayes ndi.calculator object
			%
				ndi.globals;
				w = which('ndi.calc.vis.contrast_tuning');
				parparparpar = fileparts(fileparts(fileparts(fileparts(w))));                
				oridir_bayes_obj = oridir_bayes_obj@ndi.calculator(session,'oridir_bayes',...
					fullfile(parparparpar,'ndi_common','database_documents','calc','oridir_bayes_calc.json'));
		end; % oridir_bayes() creator

		function doc = calculate(ndi_calculator_obj, parameters)
			% CALCULATE - perform the calculator for
			% ndi.calc.oridir_bayes
			%
			% DOC = CALCULATE(NDI_CALCULATION_OBJ, PARAMETERS)
			%
			% Creates a oridir_bayes_direction_tuning_calc document given input parameters.
			
				% Step 1. Check inputs
				if ~isfield(parameters,'input_parameters'),
					error(['parameters structure lacks ''input_parameters.''']);
				end;
				if ~isfield(parameters,'depends_on'),
					error(['parameters structure lacks ''depends_on.''']);
				end;
						
				% Step 2. Set up output structure
						
				tuning_doc = ndi_calculator_obj.session.database_search(ndi.query('base.id',...
					'exact_string', vlt.db.struct_name_value_search(parameters.depends_on,'stimulus_tuningcurve_id'),''));
				if numel(tuning_doc)~=1, 
					error(['Could not find stimulus tuning curve doc..']);
				end;
				tuning_doc = tuning_doc{1};

				tapp = ndi.app.stimulus.tuning_response(ndi_calculator_obj.session);
				tuning_doc = tapp.tuningdoc_fixcellarrays(tuning_doc);
				element_id = tuning_doc.dependency_value('element_id');
				num_trials = [];


				for i=1:numel(tuning_doc.document_properties.stimulus_tuningcurve.individual_responses_real),
					ind{i} = tuning_doc.document_properties.stimulus_tuningcurve.individual_responses_real{i} + ...
						sqrt(-1)*tuning_doc.document_properties.stimulus_tuningcurve.individual_responses_imaginary{i};
					ind_real{i} = ind{i};
					if any(~isreal(ind_real{i})), ind_real{i} = abs(ind_real{i}); end;
					control_ind{i} = tuning_doc.document_properties.stimulus_tuningcurve.control_individual_responses_real{i} + ...
						sqrt(-1)*tuning_doc.document_properties.stimulus_tuningcurve.control_individual_responses_imaginary{i};
					control_ind_real{i} = control_ind{i};
					if any(~isreal(control_ind_real{i})), control_ind_real{i} = abs(control_ind_real{i}); end;
					response_ind{i} = ind{i} - control_ind{i};
					response_mean(i) = nanmean(response_ind{i});
					if ~isreal(response_mean(i)), response_mean(i) = abs(response_mean(i)); end;
					response_stddev(i) = nanstd(response_ind{i});
					response_stderr(i) = vlt.data.nanstderr(response_ind{i});
					if any(~isreal(response_ind{i})),
						response_ind{i} = abs(response_ind{i});
					end;
					num_trials(i) = numel(ind{i});
				end;

				data.num_trials = num_trials(:);
				data.mean_responses = response_mean(:);
				data.angles = tuning_doc.document_properties.stimulus_tuningcurve.independent_variable_value(:);
						
				% Step 3. Calculate oridir_bayes

				[output_struct] = vis.bayes.double_gaussian.grid_proportional_noise(parameters.input_parameters.grid,...
					data,parameters.input_parameters.noise_model);

				doc = ndi.document('orientation_direction_bayes','orientation_direction_bayes',output_struct);

				% Step 4. Check if doc exists
				if ~isempty(doc), 
					doc = ndi.document(ndi_calculator_obj.doc_document_types{1},...
						'oridir_bayes_calc',parameters) + doc;
					doc = doc.set_dependency_value('stimulus_tuningcurve_id',tuning_doc.id());
					doc = doc.set_dependency_value('element_id',element_id);
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
				grid = struct('Rp',logspace(log10(.5),log10(150),30), ...
					'Op',0:5:359, ...
					'Alpha',linspace(0,1,20), ...
					'Sig',logspace(log10(5),log10(90),10), ...
					'Rsp',logspace(log10(0.01),log10(30),10),...
					'di_bins', 0:0.05:1.0, ...
					'oi_bins', 0:0.05:1.0, ...
					'dir_cv_bins', 0:0.05:1.0,...
					'cv_bins', 0:0.05:1.0);
				parameters.input_parameters.grid = grid;
				parameters.input_parameters.noise_model = [ 0.25 0.73 ];
				parameters.depends_on = vlt.data.emptystruct('name','value');
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
				q1 = ndi.query('','isa','stimulus_tuningcurve','');
				q2 = ndi.query('stimulus_tuningcurve.independent_variable_label','contains_string','irection');
				q3 = ndi.query('stimulus_tuningcurve.independent_variable_label','contains_string','eal');
				q_total = q1 & q2 & q3;
			
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

	end; % methods()
			
end % oridir_bayes
