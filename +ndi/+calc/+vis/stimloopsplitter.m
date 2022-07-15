classdef stimloopsplitter < ndi.calculator
	methods

		function stimloopsplitter_obj = stimloopsplitter(session)
			% stimloopsplitter - an app that creates stimulus_presentation documents by splitting stimuli into substimuli 
			%
			% stimloopsplitter_OBJ = stimloopsplitter(SESSION)
			%
			% Creates a stimloopsplitter ndi.calculator object
			%
				ndi.globals;
				stimloopsplitter_obj = stimloopsplitter_obj@ndi.calculator(session,'stimloopsplitter_calc',...
					fullfile(ndi_globals.path.documentpath,'apps','calculators','stimloopsplitter_calc.json'));
		end; % stimloopsplitter()

		function doc = calculate(ndi_calculator_obj, parameters)
			% CALCULATE - perform the calculator for ndi.calc.example.stimloopsplitter
			%
			% DOC = CALCULATE(NDI_CALCULATOR_OBJ, PARAMETERS)
			%
			% Creates a stimloopsplitter_calc document given input parameters.
			%
			% The document that is created stimloopsplitter
			% by the input parameters.
				% check inputs
				if ~isfield(parameters,'input_parameters'), error(['parameters structure lacks ''input_parameters''.']); end;
				if ~isfield(parameters,'depends_on'), error(['parameters structure lacks ''depends_on''.']); end;
				
				% Step 1: set up the output structure
				stimloopsplitter_calc = parameters;

                % CHANGE to get the stimulus presentation document
				stim_response_doc = ndi_calculator_obj.session.database_search(ndi.query('ndi_document.id','exact_number',...
					vlt.db.struct_name_value_search(parameters.depends_on,'stimulus_response_scalar_id'),''));
				if numel(stim_response_doc)~=1, 
					error(['Could not find stimulus response doc..']);
				end;
				stim_response_doc = stim_response_doc{1};
			
				% Step 2: perform the calculator, which here creates a tuning curve from instructions

                % loop through fields of stimulus presentation
                
                new_stim_pres_struct.presentation_order = [ repmat([stims],reps,1) ];
                new_stim_pres_struct.presentation_time = vlt.data.emptystruct('clocktype','stimopen','onset','offset','stimclose','stimevents');
                new_stim_pres_struct.stimuli = vlt.data.emptystruct('parameters');

                % here use information in the old stim_presentation_doc to
                % build the new one

				doc = ndi.document(ndi_calculator_obj.doc_document_types{1},'stimloopsplitter_calc',stimloopsplitter_calc,...
                   'stimulus_presentation',new_stim_pres_struct);

				if numel(doc)==1,
					doc = doc{1};
				end;

		end; % calculate

		function parameters = default_search_for_input_parameters(ndi_calculator_obj)
			% DEFAULT_SEARCH_FOR_INPUT_PARAMETERS - default parameters for searching for inputs
			%
			% PARAMETERS = DEFAULT_SEARCH_FOR_INPUT_PARAMETERS(NDI_CALCULATOR_OBJ)
			%
			% Returns a list of the default search parameters for finding appropriate inputs
			% to the calculator. For stimloopsplitter_calc, there is no appropriate default parameters
			% so this search will yield empty.
			%
				parameters.input_parameters = struct('independent_label','','independent_parameter','','best_algorithm','empirical_maximum');
				parameters.input_parameters.selection = vlt.data.emptystruct('property','operation','value');
				parameters.depends_on = vlt.data.emptystruct('name','value');
				parameters.query = ndi_calculator_obj.default_parameters_query(parameters);
				parameters.query(end+1) = struct('name','will_fail','query',...
					ndi.query('ndi_document.id','exact_string','123',''));
					
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
			% For the ndi.calc.stimulus.stimloopsplitter_calc class, this looks for 
			% documents of type 'stimulus_response_scalar.json' with 'response_type' fields
			% the contain 'mean' or 'F1'.
			%
			%
				q1 = ndi.query('','isa','stimulus_presentation.json','');
				q_total = q1;

				query = struct('name','stimulus_presentation_id','query',q_total);
		end; % default_parameters_query()

		function b = is_valid_dependency_input(ndi_calculator_obj, name, value)
			% IS_VALID_DEPENDENCY_INPUT - is a potential dependency input actually valid for this calculator?
			%
			% B = IS_VALID_DEPENDENCY_INPUT(NDI_CALCULATOR_OBJ, NAME, VALUE)
			%
			% Tests whether a potential input to a calculator is valid.
			% The potential dependency name is provided in NAME and its ndi_document id is
			% provided in VALUE.
			%
			% The base class behavior of this function is simply to return true, but it
			% can be overriden if additional criteria beyond an ndi.query are needed to
			% assess if a document is an appropriate input for the calculator.
			%
				b = 1;
				return;
				% the below is wrong..this function does not take stimloopsplitters or stimloopsplitter_calc objects as inputs
				q1 = ndi.query('ndi_document.id','exact_string',value,'');
				q2 = ndi.query('','isa','stimloopsplitter_calc.json','');
				% can't also be a stimloopsplitter_calc document or we could have infinite recursion
				b = isempty(ndi_calculator_obj.session.database_search(q1&q2));
		end; % is_valid_dependency_input()

		function doc_about(ndi_calculator_obj)
			% ----------------------------------------------------------------------------------------------
			% NDI_CALCULATOR: stimloopsplitter_CALC
			% ----------------------------------------------------------------------------------------------
			%
			%   ------------------------
			%   | stimloopsplitter_CALC -- ABOUT |
			%   ------------------------
			%
			%   stimloopsplitter_CALC is a demonstration document. It simply produces the 'answer' that
			%   is provided in the input parameters. Each stimloopsplitter_CALC document 'depends_on' an
			%   NDI daq system.
			%
			%   Definition: apps/stimloopsplitter_calc.json
			%
				eval(['help ndi.calc.example.stimloopsplitter.doc_about']);
		end; %doc_about()

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

				tc = doc.document_properties.tuning_curve; % shorten our typing

				% if more than 2-d, complain
				
				if numel(tc.independent_variable_label)>2,
					a = axis;
					h.objects(end+1) = text(mean(a(1:2)),mean(a(3:4)),['Do not know how to plot with more than 2 independent axes.']);
					return;
				end;

				if numel(tc.independent_variable_label)==1,
					hold on;
					h_baseline = plot([min(tc.independent_variable_value) max(tc.independent_variable_value)],...
						[0 0],'k--','linewidth',1.0001);
					h_baseline.Annotation.LegendInformation.IconDisplayStyle = 'off';
					h.objects(end+1) = h_baseline;
					net_responses = tc.response_mean - tc.control_response_mean;
					[v,sortorder] = sort(tc.independent_variable_value);
					h_errorbar = errorbar(tc.independent_variable_value(sortorder(:)),...
						tc.response_mean(sortorder(:)),tc.response_stderr(sortorder(:)),tc.response_stderr(sortorder(:)));
					set(h_errorbar,'color',[0 0 0],'linewidth',1);
					h.objects = cat(2,h.objects,h_errorbar);
					if ~h.params.suppress_x_label,
						h.xlabel = xlabel(tc.independent_variable_label);
					end;
					if ~h.params.suppress_y_label,
						h.ylabel = ylabel(['Response (' tc.response_units ')']);
					end;
					box off;
				end;

				if numel(tc.independent_variable_label)==2,
					net_responses = tc.response_mean - tc.control_response_mean;
					first_dim = unique(tc.independent_variable_value(:,1));
					colormap = spring(numel(first_dim));
					h_baseline = plot([min(tc.independent_variable_value(:,2)) max(tc.independent_variable_value(:,2))],...
						[0 0],'k--','linewidth',1.0001);
					h_baseline.Annotation.LegendInformation.IconDisplayStyle = 'off';
					h.objects(end+1) = h_baseline;
					hold on;
					for i=1:numel(first_dim),
						indexes = find(tc.independent_variable_value(:,1)==first_dim(i));
						[v,sortorder] = sort(tc.independent_variable_value(indexes,2));
						h_errorbar = errorbar(tc.independent_variable_value(indexes(sortorder),2),...
							tc.response_mean(indexes(sortorder)),...
							tc.response_stderr(indexes(sortorder)), tc.response_stderr(indexes(sortorder)));
						set(h_errorbar,'color',colormap(i,:),'linewidth',1,...
							'DisplayName',...
							[tc.independent_variable_label{1} '=' num2str(tc.independent_variable_value(indexes(1),1))]);
						h.objects = cat(2,h.objects,h_errorbar);
					end;
					if ~h.params.suppress_x_label,
						h.xlabel = xlabel(tc.independent_variable_label{2});
					end;
					if ~h.params.suppress_y_label,
						h.ylabel = ylabel(['Response (' tc.response_units ')']);
					end;
					legend;
					box off;
				end;
		end; % plot()



		function [pva] = property_value_array(ndi_calculator_obj, stim_response_doc, property)
			% PROPERTY_VALUE_ARRAY - find all values of a stimulus property
			%
			% [PVA] = ndi.calc.stimulus.stimloopsplitter.property_value_array(NDI_CALC_STIMULUS_stimloopsplitter_OBJ, STIM_RESPONSE_DOC, PROPERTY)
			%
			% Given an ndi.document of type STIMULUS_RESPONSE_SCALAR, return all values of the parameter PROPERTY that were
			% used in the stimulus.
			%
			% Values will be returned in a cell array.
			%
			% If this function cannot find a stimulus presentation document for the STIM_RESPONSE_DOC, it produces
			% an error.
			%
				stim_pres_doc = ndi_calculator_obj.session.database_search(ndi.query('ndi_document.id', 'exact_string', ...
                                        stim_response_doc.dependency_value('stimulus_presentation_id'),''));

				if numel(stim_pres_doc)~=1, 
					error(['Could not find stimulus presentation doc for document ' stim_response_doc.id() '.']);
				end;
				stim_pres_doc = stim_pres_doc{1};

				pva = {};

				for i=1:numel(stim_pres_doc.document_properties.stimulus_presentation.stimuli),
					if isfield(stim_pres_doc.document_properties.stimulus_presentation.stimuli(i).parameters,property),
						% why not just use UNIQUE? because it cares about whether the data are numbers or strings, doesn't work on numbers
						value_here = getfield(stim_pres_doc.document_properties.stimulus_presentation.stimuli(i).parameters,property);
						match_already = 0;
						for k=1:numel(pva),
							if vlt.data.eqlen(pva{k},value_here),
								match_already = 1;
								break;
							end;
						end;
						if ~match_already,
							pva{end+1} = value_here;
						end;
					end;
				end;

		end; % property_value_array
	end; % methods()
end % stimloopsplitter
