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
                w = which('ndi.calc.vis.stimloopsplitter');
				parparparpar = fileparts(fileparts(fileparts(fileparts(w))));                
				stimloopsplitter_obj = stimloopsplitter_obj@ndi.calculator(session,'stimloopsplitter_calc',...
					fullfile(parparparpar,'ndi_common','database_documents','calc','stimloopsplitter_calc.json'));
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
                tic
				stim_presentation_doc = ndi_calculator_obj.session.database_search(ndi.query('ndi_document.id','exact_number',...
					vlt.db.struct_name_value_search(parameters.depends_on,'stimulus_presentation_id'),''));
                toc
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
                new_stim_pres_struct.presentation_time = vlt.data.emptystruct('clocktype','stimopen','onset','offset','stimclose','stimevents');%will be 0x0 struct array initially
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
                            %stimulus using the parameter_adjustment field
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
                        
                        %option 1 (stimevents is a field):
                            %divide stimevents first
                            %use stimevents to set onset, offset,
                            % stimclose, stimopen
                        if isfield(old_pres_time,'stimevents')
                            frameTriggerMarker = 1;
                            new_stim_pres_struct = add_to_pres_time_with_stimevents(ndi_calculator_obj,old_pres_time,pres_orderInd,numLoops,new_stim_pres_struct,frameTriggerMarker);
                        else
                            new_stim_pres_struct = add_to_pres_time_without_stimevents(ndi_calculator_obj,old_pres_time,pres_orderInd,numLoops,new_stim_pres_struct);
                        end
                        
%                         stimevents_old = old_pres_time(pres_orderInd).stimevents;%old stimevents field    
%                         frame_events = find(stimevents_old(:,2)==frameTriggerMarker);%get indices of whatever number marks the frame being triggered in the stimevents second column - the 2's are the frame events
%                         total_frame_count = numel(frame_events);%get the number of frame events
%                          
%                         perLoop_frame_count = round(total_frame_count/(numLoops+1));%get the number of frame events per loop, round it in case not an integer
%                         for loopInd = 1:numLoops+1
%                             %go from the first frame in the loop to the
%                             %last:
%                             loop_start = (loopInd-1)*(perLoop_frame_count)+1;%the index of frame_events that is chosen as the first step of the individual loop
%                             startInd = frame_events(loop_start);
%                             loop_stop = loopInd*perLoop_frame_count;%alternatively: loop_stop = loop_start + perLoop_frame_count - 1
%                             if loop_stop>total_frame_count
%                                 stopInd = total_frame_count;%in case there's a rounding error, avoid having an index out of the array bounds
%                             else
%                                 stopInd = frame_events(loop_stop);
%                             end
%                             %how to use startInd and stopInd?
%                             %to get each new stimevents, stimopen,
%                             %stimclose, onset, offset fields
%                             
%                             %clocktype (doesn't use startInd or stopInd)
%                             new_stim_pres_struct.presentation_time(end+1,1).clocktype = old_pres_time(pres_orderInd).clocktype;
%                             %stimevents
%                             new_stim_pres_struct.presentation_time(end,1).stimevents = stimevents_old(startInd:stopInd,:);%new stimevents is subset of old stimevents
%                             %stimopen
%                             if (loopInd>1)
%                                 new_stim_pres_struct.presentation_time(end,1).stimopen = stimevents_old(startInd,1);%stimopen is the same time as onset if it's not the first substimulus
%                             else
%                                 new_stim_pres_struct.presentation_time(end,1).stimopen = old_pres_time(pres_orderInd).stimopen;%if it's the first substimulus, keep the same stimopen
%                             end
%                             %onset
%                             new_stim_pres_struct.presentation_time(end,1).onset = stimevents_old(startInd,1);%the time that the first frame is presented in each substimulus
%                             %offset
%                             new_stim_pres_struct.presentation_time(end,1).offset = stimevents_old(stopInd,1);%the time that the last frame is presented in each substimulus
%                             %stimclose
%                             if (loopInd==(numLoops+1))
%                                 new_stim_pres_struct.presentation_time(end,1).stimclose = old_pres_time(pres_orderInd).stimclose;%if it's the last substimulus, keep the same stimclose
%                             else
%                                 new_stim_pres_struct.presentation_time(end,1).stimclose = stimevents_old(stopInd,1);%if not last substimulus, stimclose is same as offset
%                             end
%                         end
                        
                        
                        %option 2 (stimevents is not a field):
                            %stimopen starts as previous stimulus' stimclose
                            %(except first stimopen, then old stimclose)
                            %onset is:
                            %onset_old + (loopInd-1)*(offset_old-onset_old)/(numLoops+1)
                            %offset is:
                            %onset_old + loopInd*(offset_old-onset_old)/(numLoops+1)
                            %stimclose is the same as offset, except the last: 
                            %the last one is the old stimclose
                            %clocktype is clocktype_old
                            %stimevents is split by the timing
                            
                            %lowerBoundInd = 1;%initialize splitInd for use inside the loop
%                         for loopInd = 1:numLoops+1
%                             %stimopen
%                             numPresentations = numel(new_stim_pres_struct.presentation_time);
%                             if (loopInd>1)
%                                 new_stim_pres_struct.presentation_time(end+1,1).stimopen = new_stim_pres_struct.presentation_time(numPresentations,1).stimclose;
%                             else
%                                 new_stim_pres_struct.presentation_time(end+1,1).stimopen = old_pres_time(pres_orderInd).stimopen;
%                             end
%                             %clocktype
%                             new_stim_pres_struct.presentation_time(end,1).clocktype = old_pres_time(pres_orderInd).clocktype;
%                             %onset
%                             onset_old = old_pres_time(pres_orderInd).onset;
%                             offset_old = old_pres_time(pres_orderInd).offset;
%                             new_stim_pres_struct.presentation_time(end,1).onset = onset_old + (loopInd-1)*(offset_old-onset_old)/(numLoops+1);
%                             %offset
%                             new_stim_pres_struct.presentation_time(end,1).offset = onset_old + (loopInd)*(offset_old-onset_old)/(numLoops+1);
%                             %stimclose
%                             new_stim_pres_struct.presentation_time(end,1).stimclose = new_stim_pres_struct.presentation_time(end,1).offset;
%                             %stimevents
% %                             if isfield(old_pres_time,'stimevents')%checks if this field is even needed
% %                                 stimevents_old = old_pres_time(pres_orderInd).stimevents;
% %                                 stimevents_toSplit = stimevents_old(:,1);%only use the time column
% %                                 timeOfSplit = new_stim_pres_struct.presentation_time(end,1).stimclose;
% %                                 startInd = lowerBoundInd;
% %                                 endInd = numel(stimevents_toSplit);%the last index of stimevents_toSplit
% %                                 upperBoundInd = split_stimevents(ndi_calculator_obj,stimevents_toSplit,timeOfSplit,startInd,endInd);%gives the index at which to split the stimevents array
% %                                 new_stim_pres_struct.presentation_time(end,1).stimevents = stimevents_old(lowerBoundInd:upperBoundInd,:);
% %                                 lowerBoundInd = upperBoundInd+1;%for the next round of splitting
% %                             end
%                         end
%                         new_stim_pres_struct.presentation_time(end,1).stimclose = old_pres_time(pres_orderInd).stimclose;%last substimulus has same stimclose as old stimulus
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
                
                % set dependency value to the dependency value of parameters
                for i=1:numel(parameters.depends_on),
                    % this doesn't work because doc is read-only: % doc.document_properties.depends_on(i) = parameters.depends_on(i)
                    doc = doc.set_dependency_value(parameters.depends_on(i).name,parameters.depends_on(i).value);
                end;
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
            % Needs to be plugged into the superclass search for input
            % parameters
			%
				parameters.input_parameters = struct('parameter_adjustment',180,'division_parameter','loops','parameter_to_split','angle');
				parameters.input_parameters.depends_on = struct('name','stimulus_presentation_id','value','');
				parameters.depends_on = vlt.data.emptystruct('name','value');
				parameters.query = ndi_calculator_obj.default_parameters_query(parameters);
% 				parameters.query(end+1) = struct('name','will_fail','query',...
% 					ndi.query('ndi_document.id','exact_string','123',''));
					
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
			% documents of type 'stimulus_presentation.json'. 
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
                %check whether the dependency is a stimulus presentation:
                % do a database search to make sure the dependency value
                % is an id for a stimulus presentation doc
                q1 = ndi.query('ndi_document.id','exact_string',value,'');
				q2 = ndi.query('','isa','stimulus_presentation.json','');
                %if there is a document that satisfies both queries and the
                %name is 'stimulus_presentation_id', the dependency is
                %valid
                b = ~isempty(ndi_calculator_obj.session.database_search(q1&q2))...
                    &strcmp(name,'stimulus_presentation_id');
				return;
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
			%   stimloopsplitter_CALC is a calculation document. It
			%   converts a set of stimuli into a new set of
			%   stimuli where a certain stimulus is split into multiple
			%   substimuli. This enriches the data by allowing more
			%   properties' effects to be considered.
			%
			%   Definition: stimloopsplitter_calc.json
			%

            
				eval(['help ndi.calc.vis.stimloopsplitter.doc_about']);
		end; %doc_about()

		function h=plot(ndi_calculator_obj, doc_or_parameters, varargin)
            % PLOT - provide a diagnostic plot to show the results of the calculator
            %
            % H=PLOT(NDI_CALCULATOR_OBJ, DOC_OR_PARAMETERS, ...)
            %
            % Produce a plot of the angles of stimuli in the order in which they're presented.
            % To edit: consider using drawshape (see help) to fill in times where stimuli are being presented 
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
                stim_pres = doc.document_properties.stimulus_presentation; % shorten our typing
                split_param_string = doc.document_properties.stimloopsplitter_calc.input_parameters.parameter_to_split;
                option = 'time';
                switch option
                    case 'noarrow'
                        %option 1: using time as x-axis no arrows
                        xlabel('time (s)')
                        ylabel([stim_param_string ' (\circ)'])
                        axis([-inf-10 inf+10 -10 370])
                        for presInd = 1:numel(stim_pres.presentation_order)
                            stimInd = stim_pres.presentation_order(presInd);
                            if (~isfield(stim_pres.stimuli(stimInd).parameters,'isblank'))
                                onset = stim_pres.presentation_time(presInd).onset;
                                offset = stim_pres.presentation_time(presInd).offset;
                                split_param = getfield(stim_pres.stimuli(stimInd).parameters,stim_param_string);%should use split parameter instead
                                X = onset:.001:offset;
                                Y = split_param*ones(1,numel(X));
                                hold on;plot(X,Y,'k');
                            end
                        end
                    case 'order'
                        %option 2: using order as x-axis with arrows
                        y_fixed = 0;
                        x_length = 1;
                        axis([-1*x_length/2 numel(stim_pres.presentation_order)+x_length/2 -1*numel(stim_pres.presentation_order)/2 numel(stim_pres.presentation_order)/2])
                        for presInd = 1:numel(stim_pres.presentation_order)
                            stimInd = stim_pres.presentation_order(presInd);
                            if (~isfield(stim_pres.stimuli(stimInd).parameters,'isblank'))
                                split_param = getfield(stim_pres.stimuli(stimInd).parameters,stim_param_string);%should use parameter_to_split field
                                arrowplot(presInd,y_fixed,split_param,x_length, 'linecolor',[rand rand rand],'headlength',.3);
                                hold on
                            end
                        end
                    case 'time'
                        %option 3: using time as x-axis with arrows
                        cmap(1:60,:) = [ones(60,1) zeros(60,1) (linspace(0,1,60))'];
                        cmap(61:120,:) = [(linspace(1,0,60))' zeros(60,1) ones(60,1)];
                        cmap(121:180,:) = [zeros(60,1) (linspace(0,1,60))' ones(60,1)];
                        cmap(181:240,:) = [zeros(60,1) ones(60,1) (linspace(1,0,60))'];
                        cmap(241:300,:) = [(linspace(0,1,60))' ones(60,1) zeros(60,1)];
                        cmap(301:360,:) = [ones(60,1) (linspace(1,0,60))' zeros(60,1)];
                        colormap(cmap);
                        cbar = colorbar;cbar.Label.String = [split_param_string ' (degrees)'];
                        cbar.TickLabels = num2cell(linspace(0,360,11));
                        y_fixed = 0;
                        x_length = 1;
                        linethickness = 1;
                        for presInd = 1:numel(stim_pres.presentation_order)
                            stimInd = stim_pres.presentation_order(presInd);
                            if (~isfield(stim_pres.stimuli(stimInd).parameters,'isblank'))
                                onset = stim_pres.presentation_time(presInd).onset;
                                offset = stim_pres.presentation_time(presInd).offset;
                                split_param_actual = getfield(stim_pres.stimuli(stimInd).parameters,split_param_string);
                                split_param_compass = vlt.math.cartesian2compass(split_param_actual,0);
                                X = mean([onset offset]);
                                
                                hold on;
                                if (split_param_actual == 0)
                                    split_param_actual = 360;%cmap doesn't take 0 as an array index
                                end
                                arrowplot(X,y_fixed,split_param_compass,x_length,'linethickness',linethickness,'linecolor',cmap(split_param_actual,:));
                                %consider using drawshape (see help) to fill in times where stimuli are being presented 
                            end
                        end
                        ylim(xlim - diff(xlim)/2);
                        h.xlabel = xlabel('time (s)');
                        yticks([]);%y axis doesn't have any meaning, so ticks should be hidden
                        
                    case 'combined'
                        %option 4: using arrows and angle degrees with time
                        xlabel('time (s)')
                        ylabel([split_param_string ' (degrees)'])
                        axis([-inf-10 inf+10 -10 370])
                        x_length = 3;
                        linethickness = 1;
                        %set cmap to a 360 degree wheel
                        %old code: %cmap = colormap(parula(360));
                        cmap(1:60,:) = [ones(60,1) zeros(60,1) (linspace(0,1,60))'];
                        cmap(61:120,:) = [(linspace(1,0,60))' zeros(60,1) ones(60,1)];
                        cmap(121:180,:) = [zeros(60,1) (linspace(0,1,60))' ones(60,1)];
                        cmap(181:240,:) = [zeros(60,1) ones(60,1) (linspace(1,0,60))'];
                        cmap(241:300,:) = [(linspace(0,1,60))' ones(60,1) zeros(60,1)];
                        cmap(301:360,:) = [ones(60,1) (linspace(1,0,60))' zeros(60,1)];
                        colormap(cmap);colorbar;
                        for presInd = 1:numel(stim_pres.presentation_order)
                            stimInd = stim_pres.presentation_order(presInd);
                            if (~isfield(stim_pres.stimuli(stimInd).parameters,'isblank'))
                                onset = stim_pres.presentation_time(presInd).onset;
                                offset = stim_pres.presentation_time(presInd).offset;
                                split_param = getfield(stim_pres.stimuli(stimInd).parameters,split_param_string);
                                X = mean([onset offset]);
                                Y = split_param;
                                color = cmap(floor(split_param)+1,:);
                                hold on;
                                arrowplot(X,Y,split_param,x_length,'linecolor',color,...
                                    'linethickness',linethickness);
                            end
                        end
                    case 'stimuli'
                        %option 5: show each stimulus' split parameter in
                        %order of index
                        
                        split_param_vals = property_value_array(ndi_calculator_obj,doc,split_param_string);
                        y_fixed = 0;
                        x_length = 1;
                        headangle = 15;
                        for valInd = 1:numel(split_param_vals)
                            split_param = split_param_vals{valInd};
                            hold on;
                            arrowplot(valInd,y_fixed,split_param,x_length,'headangle',headangle);
                        end
                        ylim(xlim - diff(xlim)/2);
                end 
		end; % plot()

         % NEW functions in tuningcurve_calc that are not overriding any superclass functions
        function [new_stim_pres_struct] = add_to_pres_time_with_stimevents(ndi_calculator_obj,old_pres_time,pres_orderInd,numLoops,new_stim_pres_struct,frameTriggerMarker)
            % ADD_TO_PRES_TIME_WITH_STIMEVENTS - add a new row or set of
            % rows to the presentation time field in the new stimulus presentation struct, calculated from
            % a single row in the presentation time field of the old
            % stimulus presentation, assuming that field contains the
            % subfield "stimevents"
			%
			% [NEW_STIM_PRES_STRUCT] = 
			% add_to_pres_time_with_stimevents(old_pres_time,pres_orderInd,numLoops,new_stim_pres_struct,frameTriggerMarker)
			%
			% identify the starting index and ending index for each loop in the current presentation by using the
			% specified frameTriggerMarker and dividing stimevents into
			% substimuli by the number of loops.
			%
			% 
            stimevents_old = old_pres_time(pres_orderInd).stimevents;%old stimevents field    
            frame_events = find(stimevents_old(:,2)==frameTriggerMarker);%get indices of whatever number marks the frame being triggered in the stimevents second column - the 2's are the frame events
            total_frame_count = numel(frame_events);%get the number of frame events

            perLoop_frame_count = round(total_frame_count/(numLoops+1));%get the number of frame events per loop, round it in case not an integer
            for loopInd = 1:numLoops+1
                %go from the first frame in the loop to the
                %last:
                loop_start = (loopInd-1)*(perLoop_frame_count)+1;%the index of frame_events that is chosen as the first step of the individual loop
                startInd = frame_events(loop_start);
                loop_stop = loopInd*perLoop_frame_count;%alternatively: loop_stop = loop_start + perLoop_frame_count - 1
                if loop_stop>total_frame_count
                    stopInd = total_frame_count;%in case there's a rounding error, avoid having an index out of the array bounds
                else
                    stopInd = frame_events(loop_stop);
                end
                %how to use startInd and stopInd?
                %to get each new stimevents, stimopen,
                %stimclose, onset, offset fields

                %clocktype (doesn't use startInd or stopInd)
                new_stim_pres_struct.presentation_time(end+1,1).clocktype = old_pres_time(pres_orderInd).clocktype;
                %stimevents
                new_stim_pres_struct.presentation_time(end,1).stimevents = stimevents_old(startInd:stopInd,:);%new stimevents is subset of old stimevents
                %stimopen
                if (loopInd>1)
                    new_stim_pres_struct.presentation_time(end,1).stimopen = stimevents_old(startInd,1);%stimopen is the same time as onset if it's not the first substimulus
                else
                    new_stim_pres_struct.presentation_time(end,1).stimopen = old_pres_time(pres_orderInd).stimopen;%if it's the first substimulus, keep the same stimopen
                end
                %onset
                new_stim_pres_struct.presentation_time(end,1).onset = stimevents_old(startInd,1);%the time that the first frame is presented in each substimulus
                %offset
                new_stim_pres_struct.presentation_time(end,1).offset = stimevents_old(stopInd,1);%the time that the last frame is presented in each substimulus
                %stimclose
                if (loopInd==(numLoops+1))
                    new_stim_pres_struct.presentation_time(end,1).stimclose = old_pres_time(pres_orderInd).stimclose;%if it's the last substimulus, keep the same stimclose
                else
                    new_stim_pres_struct.presentation_time(end,1).stimclose = stimevents_old(stopInd,1);%if not last substimulus, stimclose is same as offset
                end
            end
        end % add_to_pres_time_with_stimevents
        function new_stim_pres_struct = add_to_pres_time_without_stimevents(ndi_calculator_obj,old_pres_time,pres_orderInd,numLoops,new_stim_pres_struct)
            % ADD_TO_PRES_TIME_WITHOUT_STIMEVENTS - add a new row or set of
            % rows to the presentation time field in the new stimulus presentation struct, calculated from
            % a single row in the presentation time field of the old
            % stimulus presentation, assuming that field does NOT contain the
            % subfield "stimevents"
			%
			% [NEW_STIM_PRES_STRUCT] = 
			% add_to_pres_time_without_stimevents(old_pres_time,pres_orderInd,numLoops,new_stim_pres_struct,frameTriggerMarker)
			%
			% use the stimopen, stimclose, onset, and offset fields in old_pres_time, as well as numLoops, to
			% calculate new stimopen, stimclose, onset, and offset fields
			% that are set in a new row or set of rows in
			% new_stim_pres_struct
			% 
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
        end % add_to_pres_time_without_stimevents
        function [upperBoundInd] = split_stimevents(ndi_calculator_obj,stimEventTimes, timeMarker, lowerInd, upperInd)
            % SPLIT_STIMEVENTS - get the upper bound index with which to
            % split the stimEventTimes array, a field found in a stimulus presentation
            % presentation_times field
			%
			% [UPPERBOUNDIND] =
			% splitStimEvents(stimEventTimes,timeMarker,lowerInd,upperInd)
			%
			% take the stimEventTimes array and find the index within lowerInd and upperInd 
			% with the closest time less than the time in timeMarker. This will be the upper
            % bound index of the newly split stimulus stimEvents array
			%
			%
                if upperInd==lowerInd %base case
                    if timeMarker>= stimEventTimes(lowerInd) 
                        upperBoundInd = lowerInd;
                        return;
                    else
                        upperBoundInd = lowerInd-1;%the array upper bound should be lower than the given time of split
                        return;
                    end
                end
                meanInd = floor(mean([lowerInd,upperInd+1]));
                if timeMarker>=stimEventTimes(meanInd)
                    upperBoundInd = split_stimevents(ndi_calculator_obj,stimEventTimes,timeMarker,meanInd,upperInd);%keep the same upperInd but switch the lowerInd to the mean index calculated earlier
                else
                    upperBoundInd = split_stimevents(ndi_calculator_obj,stimEventTimes,timeMarker,lowerInd,meanInd);%keep the same lowerInd but switch upperInd to the mean index
                end
        end %splitStimEvents

		function [pva] = property_value_array(ndi_calculator_obj, stimloopsplitter_calc_doc, property)
			% PROPERTY_VALUE_ARRAY - find all values of a stimulus property
			%
			% [PVA] = ndi.calc.stimulus.stimloopsplitter.property_value_array(NDI_CALCULATOR_OBJ, STIMLOOPSPLITTER_CALC_DOC, PROPERTY)
			%
			% Given an ndi.document of type STIMLOOPSPLITTER_CALC, return all values of the parameter PROPERTY that were
			% used in the stimulus.
			%
			% Values will be returned in a cell array.
			%
			% If this function cannot find a stimulus presentation document for the STIM_RESPONSE_DOC, it produces
			% an error.
			%
				
                stim_pres_doc = ndi_calculator_obj.session.database_search(ndi.query('ndi_document.id', 'exact_string', ...
                                        stimloopsplitter_calc_doc.dependency_value('stimulus_presentation_id'),''));
				if numel(stim_pres_doc)~=1, 
					error(['Could not find stimulus presentation doc for document ' stimloopsplitter_calc_doc.id() '.']);
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
