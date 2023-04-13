classdef oridir_bootstrap < ndi.calculator

	methods
		function oridir_bootstrap_obj = oridir_bootstrap(session)
			% oridir_bootstrap - ndi.calculator object that
			% creates bootstrapped orientation and direction tuning curves from spike
			% elements
			%
			% oridir_bootstrap_OBJ = BOOTSTRAPPING(SESSION)
			%
			% Creates a oridir_bootstrap ndi.calculator object
			%
				ndi.globals;
				w = which('ndi.calc.vis.oridir_bootstrap');
				parparparpar = fileparts(fileparts(fileparts(fileparts(w))));                
				oridir_bootstrap_obj = oridir_bootstrap_obj@ndi.calculator(session,'oridir_bootstrap',...
					fullfile(parparparpar,'ndi_common','database_documents','calc','oridir_bootstrap_calc.json'));
		end % oridir_bootstrap() creator

		function doc = calculate(ndi_calculator_obj, parameters)
			% CALCULATE - perform the calculator for ndi.calc.example.oridir_bootstrap
			%
			% DOC = CALCULATE(NDI_CALCULATOR_OBJ, PARAMETERS)
			%
			% Creates a oridir_bootstrap_calc document given input parameters.
			%
			% The document that is created bootstrap
			% by the input parameters.
				% check inputs
				if ~isfield(parameters,'input_parameters'), error('parameters structure lacks ''input_parameters''.'); end
				if ~isfield(parameters,'depends_on'), error('parameters structure lacks ''depends_on''.'); end
			
				% Step 1: set up the output structure
				oridir_bootstrap_calc = parameters;
                
				%looking for oridir_tuning docs
				index = strcmp('orientation_direction_tuning_id',{parameters.depends_on.name});

				if isempty(index)
					error('No orientation_direction_tuning_id specified.');
				end

				orientation_direction_tuning_id = parameters.depends_on(index).value;

				stim_response_doc = ndi_calculator_obj.session.database_search(ndi.query('ndi_document.id','exact_string',orientation_direction_tuning_id,''));

				if numel(stim_response_doc)~=1 
					error('Could not find stimulus response doc..');
				end
				stim_response_doc = stim_response_doc{1};

				ndi.globals
			
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
				% make a local store of our global variables
				logger = ndi_globals.log;

				parfor i = 1:niter
					logger.msg('debug',5,['oridir_bootstrap: Analyzing simulation ' int2str(i) ' of ' int2str(niter) '.']);
					%prepare the inputs for the fitting function
					bootstrap_res = data2fit(:,:,i)'; %isolate a single bootstrap iteration
					response_mean = mean(bootstrap_res); %get the mean response at each stimuli
					response_stddev = std(bootstrap_res); %std 
					response_stderr = response_stddev/sqrt(ntrials); %SEM
					%build the input struct
					response = [];
					response.curve = [stimuli'; ...
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
				end % looping over bootstraps

				%set outputs
				parameters_here = struct();
				fitparams = struct('oridir_bootstrap_fit_params',dp);
				fitparams.values = '[Rsp Rp Op sigm Rn]';
                
				%add outputs to a struct that will be loaded into the doc
				parameters_here.fit_parameters = fitparams;
				parameters_here.input_parameters = parameters;
				parameters_here.oridir_bootstrap_fits = fits;
				%create oridir_bootstrap ndi document
				doc = ndi.document(ndi_calculator_obj.doc_document_types{1},'oridir_bootstrap_calc', parameters_here);

				% set dependency value to the dependency value of parameters
				for i=1:numel(parameters.depends_on)
					doc = doc.set_dependency_value(parameters.depends_on(i).name,parameters.depends_on(i).value);
				end
				if numel(doc)~=1
					doc = doc{1};
				end
		end % calculate method

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
			% For the ndi.calc.vis.oridir_bootstrap_calc class, this looks for 
			% documents of type 'oridir_tuning_calc.json' with 'response_type' fields
			% the contain 'mean'.
			%
			%
				q1 = ndi.query('','isa','orientation_direction_tuning','');
				q2_a = ndi.query('orientation_direction_tuning.properties.response_type','contains_string','mean','');
				q2_b = ndi.query('orientation_direction_tuning.properties.response_type','contains_string','F1','');
				q3 = ndi.query('orientation_direction_tuning.tuning_curve.direction','hasmember',270,'');
				q4 = ndi.query('orientation_direction_tuning.significance.visual_response_anova_p','lessthan',0.05,'');
				q_total = q1 & (q2_a | q2_b) & q3 & q4;

				query = struct('name','orientation_direction_tuning_id','query',q_total);
		end % default_parameters_query()

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
				%plot gaussian for each iteratio

				% call superclass plot method to set up axes
				h=plot@ndi.calculator(ndi_calculator_obj, doc_or_parameters, varargin{:});

				if isa(doc_or_parameters,'ndi.document'),
					doc = doc_or_parameters;
				else,
					error(['Do not know how to proceed without an ndi document for doc_or_parameters.']);
				end

				bsc = doc_or_parameters.document_properties.oridir_bootstrap_calc; % shorten our typing
                
				ha = plot(bsc.oridir_bootstrap_fits','Color',[0,0,0]+0.7);
                
				for i = 1:(size(ha,1))
					h.objects(end+1) = ha(i,1);
				end

				if ~h.params.suppress_x_label,
					h.xlabel = xlabel('Direction (\circ)');
				end;
				if ~h.params.suppress_y_label,
					h.ylabel = ylabel('Response (Hz)');
				end;

				if ~h.params.suppress_title,
					h.title = title(doc_or_parameters.document_properties.ndi_document.id);
				end;
                
				h.axes.XLim =[0 360];
		end % plot()

	end % methods()
end % oridir_bootstrap

