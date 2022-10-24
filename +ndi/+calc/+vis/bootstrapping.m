classdef bootstrapping < ndi.calculator

	methods
        function bootstrapping_obj = bootstrapping(session)
			% bootstrapping - ndi.calculator object that
			% creates bootstrapped orientation and direction tuning curves from spike
			% elements
			%
			% bootstrapping_OBJ = BOOTSTRAPPING(SESSION)
			%
			% Creates a bootstrapping ndi.calculator object
			%
				ndi.globals;
				w = which('ndi.calc.vis.bootstrapping');
				parparparpar = fileparts(fileparts(fileparts(fileparts(w))));                
				bootstrapping_obj = bootstrapping_obj@ndi.calculator(session,'bootstrapping',...
					fullfile(parparparpar,'ndi_common','database_documents','calc','bootstrapping_calc.json'));
		end % bootstrapping() creator

		function doc = calculate(ndi_calculator_obj, parameters)
			% CALCULATE - perform the calculator for
            % ndi.calc.example.bootstrapping
			%
			% DOC = CALCULATE(NDI_CALCULATOR_OBJ, PARAMETERS)
			%
			% Creates a bootstrapping_calc document given input parameters.
			%
			% The document that is created bootstrap
			% by the input parameters.
				% check inputs
				if ~isfield(parameters,'input_parameters'), error('parameters structure lacks ''input_parameters''.'); end
				if ~isfield(parameters,'depends_on'), error('parameters structure lacks ''depends_on''.'); end
			
				% Step 1: set up the output structure
				bootstrapping_calc = parameters;
                
                %looking for oridir_tuning docs
                index = strcmp('orientation_direction_tuning_id',{parameters.depends_on.name});

                if isempty(index)
                    error('No orientation_direction_tuning_id specified.');
                end

                orientation_direction_tuning_id = parameters.depends_on(index).value;

				stim_response_doc = ndi_calculator_obj.session.database_search(ndi.query('ndi_document.id','exact_string',orientation_direction_tuning_id,''));
keyboard
 				if numel(stim_response_doc)~=1 
 					error('Could not find stimulus response doc..');
 				end
				stim_response_doc = stim_response_doc{1};
			
				% Step 2: perform the calculator, which here creates a tuning curve from instructions

				% extract individual trial stimulus responses
				responses = stim_response_doc.document_properties.orientation_direction_tuning.tuning_curve.individual;
                stimuli = stim_response_doc.document_properties.orientation_direction_tuning.tuning_curve.direction;
                
                %the data we will fit
                data2fit = responses;
                
                %gather indices from stim response doc
                ntrials = size(responses,2);
                nstim  = size(responses,1);
                niter = parameters.input_parameters.iterations;

                %create random inds
                randind = ceil(ntrials*rand(ntrials,nstim,niter));

                %assign response values
                for i_iter = 1:niter
                      myinds = (randind(:,:,i_iter)-1)*nstim+repmat(1:nstim,ntrials,1);
                      data2fit(:,:,i_iter) = responses(myinds)';
                end

                dp = zeros(niter,5); % [Rsp Rp Op sigm Rn] x niter
                fits = zeros(niter,360); %save fits for plotting
                %perform fits for each bootstrap
                for i = 1:niter
                    %prepare the inputs for the fitting function
                    bootstrap_res = data2fit(:,:,i)'; %isolate a single bootstrap iteration
                    response_mean = mean(bootstrap_res); %get the mean response at each stimuli
                    response_stddev = std(bootstrap_res); %std 
                    response_stderr = response_stddev/sqrt(ntrials); %SEM
                    %build the input struct
 				    response.curve = ...
 					    [stimuli'; ...
 						    response_mean; ...
 						    response_stddev; ...
 						    response_stderr];
 				    response.ind = bootstrap_res;
                    %run the fit
       				fi = vlt.neuro.vision.oridir.index.oridir_fitindexes(response);
                    %extract [Rsp Rp Op sigm Rn] to add to NDI doc
                    %they are the 5 values under fi.fit_parameters
                    dp(i,:) = [fi.fit_parameters];
                    fits(i,:)= [fi.fit(2,:)];
                end

                %set oui_itertputs
                parameters_here = struct();
                fitparams = struct('bootstrapping_fit_params',dp);
                fitparams.values = '[Rsp Rp Op sigm Rn]';
                
                %add outputs to a struct that will be loaded into the doc
                parameters_here.fit_parameters = fitparams;
                parameters_here.orientation_direction_tuning = stim_response_doc.document_properties.orientation_direction_tuning;
                parameters_here.input_parameters = parameters;
                parameters_here.bootstrapping_fits = fits;
                %create bootstrapping ndi document
                doc = ndi.document(ndi_calculator_obj.doc_document_types{1},'bootstrapping_calc', parameters_here);

                % set dependency value to the dependency value of parameters
                for i=1:numel(parameters.depends_on)
                     doc = doc.set_dependency_value(parameters.depends_on(i).name,parameters.depends_on(i).value);
                end
				if numel(doc)~=1
					doc = doc{1};
                end

        end % calculate

		function parameters = default_search_for_input_parameters(ndi_calculator_obj)
            % DEFAULT_SEARCH_FOR_INPUT_PARAMETERS - default parameters for searching for inputs
            %
            % PARAMETERS = DEFAULT_SEARCH_FOR_INPUT_PARAMETERS(NDI_CALCULATION_OBJ)
            %
            % Returns a list of the default search parameters for finding appropriate inputs
            % to the calculator.
            %
            % search for stimulus_tuningcurve_id
            parameters.input_parameters.iterations = 100;
            parameters.input_parameters.depends_on = vlt.data.emptystruct('name','value');
            
            % parameters.input_parameters = struct('independent_label','','independent_parameter','','best_algorithm','empirical_maximum');
            % parameters.input_parameters.selection = vlt.data.emptystruct('property','operation','value');
            
            parameters.depends_on = vlt.data.emptystruct('name','value');
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
			% For the ndi.calc.vis.bootstrapping_calc class, this looks for 
			% documents of type 'oridir_tuning_calc.json' with 'response_type' fields
			% the contain 'mean'.
			%
			%

				q1 = ndi.query('','isa','orientation_direction_tuning','');
				q2 = ndi.query('orientation_direction_tuning.properties.response_type','contains_string','mean','');
                q_total = q1 & q2;

				query = struct('name','orientation_direction_tuning_id','query',q_total);
		end % default_parameters_query()

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
	            %look for oridirtuning_calc documents
				q1 = ndi.query('ndi_document.id','exact_string',value,'');
				q2 = ndi.query('','isa','oridirtuning_calc.json','');
				b = isempty(ndi_calculator_obj.session.database_search(q1&q2));
		end % is_valid_dependency_input()

		function doc_about(ndi_calculator_obj)
			% ----------------------------------------------------------------------------------------------
			% NDI_CALCULATOR: BOOTSTRAPPING_CALC
			% ----------------------------------------------------------------------------------------------
			%
			%   -------------------------------
			%   | BOOTSTRAPPING_CALC -- ABOUT |
			%   -------------------------------
			%
			%  bootstrapping_calc is an NDI calculator that creates n
            %  shuffled trials from oridir_tuningcurves
            %  using the bootstrp technique. After bootstrapping, double 
            %  gaussian fits are performed on each trial. The data is
            %  output into a boostrapping document.
			%
			%   Definition: apps/bootstrapping_calc.json
			%
				eval(['help ndi.calc.example.bootstrapping.doc_about']);
		end %doc_about()

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

                        %Rsp+Rp*EXP(-(X-Op)^2)/(2*sig^2))+Rn*EXP(-(X-Op+180)^2/(2*sig^2))
                        %plot gaussian for each iteration
keyboard

				% call superclass plot method to set up axes
				h=plot@ndi.calculator(ndi_calculator_obj, doc_or_parameters, varargin{:});

				if isa(doc_or_parameters,'ndi.document'),
					doc = doc_or_parameters;
				else,
					error(['Do not know how to proceed without an ndi document for doc_or_parameters.']);
				end

				bsc = doc_or_parameters.document_properties.bootstrapping_calc; % shorten our typing
                

                h.xlabel = 'Direction';






		end % plot()

		
	end % methods()
end % bootstrapping
