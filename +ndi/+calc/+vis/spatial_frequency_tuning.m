classdef spatial_frequency_tuning < ndi.calculator

	methods
		function spatial_frequency_tuning_obj = spatial_frequency_tuning(session)
			% SPATIAL_FREQUENCY_TUNING - a spatial_frequency_tuning demonstration of an ndi.calculator object
			%
			% SPATIAL_FREQUENCY_TUNING_OBJ = SPATIAL_FREQUENCY_TUNING(SESSION)
			%
			% Creates a SPATIAL_FREQUENCY_TUNING ndi.calculator object
			%
				ndi.globals;
				w = which('ndi.calc.vis.spatial_frequency_tuning');
				parparparpar = fileparts(fileparts(fileparts(fileparts(w))));
				spatial_frequency_tuning_obj = spatial_frequency_tuning_obj@ndi.calculator(session,'spatial_frequency_tuning_calc',...
					fullfile(parparparpar,'ndi_common','database_documents','calc','spatial_frequency_tuning_calc.json'));
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
					vlt.db.struct_name_value_search(parameters.depends_on,'stimulus_tuningcurve_id'),''));
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
				did.datastructures.assign(varargin{:});
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

				sf_props = vis.spatial_frequency_analysis(resp);

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
				abs_sf_props = vis.spatial_frequency_analysis(resp_abs);

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

	end; % methods()
end % spatial_frequency_tuning
