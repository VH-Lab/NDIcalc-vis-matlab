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

                % CHANGE to get the stimulus presentation document: 
                %   1) get the stimulus presentation ID specified by
                %   parameters
                %   2) search for the specified stimulus_presentation
                %   document
				stim_presentation_doc = ndi_calculator_obj.session.database_search(ndi.query('ndi_document.id','exact_number',...
					vlt.db.struct_name_value_search(parameters.depends_on,'stimulus_presentation_id'),''));
				if numel(stim_presentation_doc)~=1, 
					error(['Could not find stimulus presentation doc..']);
				end;
				stim_presentation_doc = stim_presentation_doc{1};
			
				% Step 2: perform the calculator
                %   1) create new stimulus presentation struct with empty fields
                %   2) fill out fields based on old stimulus presentation
                %   3) turn new stimulus presentation struct into a
                %   document

                % loop through fields of stimulus presentation
                
                new_stim_pres_struct.presentation_order = [];
                new_stim_pres_struct.presentation_time = vlt.data.emptystruct('clocktype','stimopen','onset','offset','stimclose');%will be 0x0 struct array initially
                new_stim_pres_struct.stimuli = vlt.data.emptystruct('parameters');%will also be a 0x0 struct array initially 

                % here use information in the old stim_presentation_doc to
                % build the new one
                
                % go through each stimulus to check if it has loops
                % if it does, need to add stimuli, presentation_time, and
                % presentation_order (not the same amount necessarily)
                
                old_stim_pres_struct = stim_presentation_doc.document_properties.stimulus_presentation;%use stim_presentation_doc
                old_stimuli = old_stim_pres_struct.stimuli;
                old_pres_order = old_stim_pres_struct.presentation_order;
                old_pres_time = old_stim_pres_struct.presentation_time;
                new_stim_pres_struct.stimuli = old_stimuli;%initially will have the same number of stimuli
                
                stimulus_pairs = zeros(length(old_stimuli),1);%populated with either zero or the stimulus ID of the new stimulus created. Index corresponds to old_stim_pres_struct index
                for pres_orderInd = 1:numel(old_pres_order)%loop through presentation order
                    stimInd = old_pres_order(pres_orderInd);%presentation_order holds stimulus index
                    %get number of division parameters for current stimulus
                    stimulus = old_stimuli(stimInd);
                    %make a new helper method starting here? Can replace
                    %with another helper method if division parameter
                    %changes - the rest of this is dependent on loop being
                    %the division parameter
                    
                    %%% first need to check if the stimulus is not a
                    %%% control, because then it might not have the loops
                    %%% field (or whatever division parameter being used)
                    division_paramExists = isfield(stimulus.parameters,char(stimloopsplitter_calc.input_parameters.division_parameter)); %division parameter used in case we want to divide by something other than loops
                    %get number of loops using the division parameter
                    if (division_paramExists)
                        numLoops = eval(['stimulus.parameters.',char(stimloopsplitter_calc.input_parameters.division_parameter)]);%division parameter used in case we want to divide by something other than loops
                    end
                    %check if number of loops is greater than 0
                    if (~division_paramExists || numLoops==0)
                        %no difference here between old and new pres order and time
                        %fields
                        new_stim_pres_struct.presentation_order(end+1,1) = old_pres_order(pres_orderInd);
                        new_stim_pres_struct.presentation_time(end+1,1) = old_pres_time(pres_orderInd);
                    elseif (numLoops>0)
                        %check for pairing
                        if (stimulus_pairs(stimInd)==0)%0 means no pairing since no stimulus has ID 0
                            %create new stimulus by adjusting the current
                            %stimulus
                            eval(['stimulus.parameters.',char(stimloopsplitter_calc.input_parameters.parameter_to_split),'=',...
                            'stimulus.parameters.',char(stimloopsplitter_calc.input_parameters.parameter_to_split),'+',...
                            char(num2str(stimloopsplitter_calc.input_parameters.parameter_adjustment))]);
                            stimulus_adjusted = stimulus;
                            new_stim_pres_struct.stimuli(end+1)= stimulus_adjusted;%newly formed stimulus added to list of stimuli
                            %add it to pairing
                            new_stimInd = numel(new_stim_pres_struct.stimuli);%assumes the new stimulus is added at the end of the stimuli array
                            stimulus_pairs(stimInd)=new_stimInd;
                            %set loops to 0 for both stimuli in
                            %new stim pres struct
                            eval(['new_stim_pres_struct.stimuli(stimInd).parameters.',char(stimloopsplitter_calc.input_parameters.division_parameter),'=0']);
                            eval(['new_stim_pres_struct.stimuli(new_stimInd).parameters.',char(stimloopsplitter_calc.input_parameters.division_parameter),'=0']);
                        else %there is a pairing that's already been added - not sure there's anything to do 
                        end
                        %now set new presentation order and time fields:
                        %add to pres_order in correct order (old then new, alternating based on number of loops)
                        
                        new_stim_pres_struct.presentation_order(end+1,1) = stimInd; %stimulus index of stimulus that already existed
                        for loop = 1:numLoops %only add one more stimulus if loop = 1 (assumes no loops means loops = 0)
                            if mod(loop,2)==1 %if it's the 1st,3rd,etc. loop 
                                new_stimInd = stimulus_pairs(stimInd);
                                new_stim_pres_struct.presentation_order(end+1,1) = new_stimInd; %stimulus index of newly created stimulus
                            else
                                new_stim_pres_struct.presentation_order(end+1,1) = stimInd; %stimulus index of stimulus that already existed
                            end
                        end
                        %add to pres_time
                        %stimopen starts as previous stimulus' stimclose
                        %(except first stimopen, then old stimclose)
                        %onset is:
                        %onset_old + (loopInd-1)*(offset_old-onset_old)/(numLoops+1)
                        %offset is:
                        %onset_old + loopInd*(offset_old-onset_old)/(numLoops+1)
                        %stimclose is the same as offset, except the last: 
                        %the last one is the old stimclose
                        %clocktype is clocktype_old
                        for loopInd = 1:numLoops+1
                            %stimopen
                            numPresentations = numel(new_stim_pres_struct.presentation_time);
                            if (loopInd>1)
                                new_stim_pres_struct.presentation_time(end+1,1).stimopen = new_stim_pres_struct.presentation_time(numPresentations,1).stimclose;
                            else
                                new_stim_pres_struct.presentation_time(end+1,1).stimopen = old_pres_time(pres_orderInd).stimopen;
                            end
                            %clocktype
                            new_stim_pres_struct.presentation_time(end,1).clocktype = old_pres_time(pres_orderInd).clocktype;
                            %onset
                            onset_old = old_pres_time(pres_orderInd).onset;
                            offset_old = old_pres_time(pres_orderInd).offset;
                            new_stim_pres_struct.presentation_time(end,1).onset = onset_old + (loopInd-1)*(offset_old-onset_old)/(numLoops+1);
                            %offset
                            new_stim_pres_struct.presentation_time(end,1).offset = onset_old + (loopInd)*(offset_old-onset_old)/(numLoops+1);
                            %stimclose
                            new_stim_pres_struct.presentation_time(end,1).stimclose = new_stim_pres_struct.presentation_time(end,1).offset;
                        end
                        new_stim_pres_struct.presentation_time(end,1).stimclose = old_pres_time(pres_orderInd).stimclose;%last substimulus has same stimclose as old stimulus
                    else
                        error(['number of loops not valid']);
                    end
                end
%                 for stimInd = 1:numel(old_stimuli)
%                     stimulus = old_stimuli(stimInd);
%                     loopNum = eval(['stimulus.parameters.',stimloopsplitter_calc.input_parameters.division_parameter]);%division parameter used in case we want to divide by something other than loops
%                     if (loopNum~=0)%are there loops?
%                         %split by loops: add 1 stimulus, then add
%                         %presentation_times and presentation_orders by
%                         %calling a helper method
%                         stimulus_adjusted = eval(['stimulus.parameters.',stimloopsplitter_calc.input_parameters.parameter_to_split,'=',...
%                             'stimulus.parameters.',stimloopsplitter_calc.input_parameters.parameter_to_split,'+',...
%                             stimloopsplitter_calc.input_parameters.parameter_adjustment]);
%                         new_stim_pres_struct.stimuli(end+1)= stimulus_adjusted;%newly formed stimulus added to list of stimuli
%                         
%                     end
%                 end
				doc = ndi.document(ndi_calculator_obj.doc_document_types{1},'stimloopsplitter_calc',stimloopsplitter_calc,...
                   'stimulus_presentation',new_stim_pres_struct);

				if numel(doc)~=1,
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
                        % Produce a plot of the angles of stimuli in the order in which they're presented.
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
                sp = doc.document_properties.stimulus_presentation; % shorten our typing
                option = 'combined';
                switch option
                    case 'noarrow'
                        %option 1: using time as x-axis no arrows
                        for presInd = 1:numel(sp.presentation_order)
                            stimInd = sp.presentation_order(presInd);
                            if (~isfield(sp.stimuli(stimInd).parameters,'isblank'))
                                onset = sp.presentation_time(presInd).onset;
                                offset = sp.presentation_time(presInd).offset;
                                angle = sp.stimuli(stimInd).parameters.angle;
                                X = onset:.001:offset;
                                Y = angle*ones(1,numel(X));
                                hold on;plot(X,Y,'k');
                            end
                        end
                    case 'order'
                        %option 2: using order as x-axis with arrows
                        y_fixed = 0;
                        x_length = .1;
                        axis([-1*x_length/2 numel(sp.presentation_order)+x_length/2 -1*numel(sp.presentation_order)/2 numel(sp.presentation_order)/2])
                        for presInd = 1:numel(sp.presentation_order)
                            stimInd = sp.presentation_order(presInd);
                            if (~isfield(sp.stimuli(stimInd).parameters,'isblank'))
                                angle = sp.stimuli(stimInd).parameters.angle;
                                arrowplot(presInd,y_fixed,angle,x_length, 'linecolor',[rand rand rand],'headlength',1);
                                hold on
                            end
                        end
                    case 'time'
                        %option 3: using time as x-axis with arrows
                    case 'combined'
                        %option 4: using arrows and angle degrees
                        x_length = .1;
                        cmap = colormap(parula(360));
                        for presInd = 1:numel(sp.presentation_order)
                            stimInd = sp.presentation_order(presInd);
                            if (~isfield(sp.stimuli(stimInd).parameters,'isblank'))
                                onset = sp.presentation_time(presInd).onset;
                                offset = sp.presentation_time(presInd).offset;
                                angle = sp.stimuli(stimInd).parameters.angle;
                                X = mean([onset offset]);
                                Y = angle;
                                color = cmap(floor(angle)+1,:);
                                hold on;arrowplot(X,Y,angle,x_length,'linecolor',color);
                            end
                        end
                end 
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
