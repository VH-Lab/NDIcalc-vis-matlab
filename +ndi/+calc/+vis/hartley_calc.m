classdef hartley_calc < ndi.calculator

    methods
        function hartley_calc_obj = hartley_calc(session)
            % HARTLEY_CALC - a hartley_calc demonstration of an ndi.calculator object
            %
            % HARTLEY_CALC_OBJ = HARTLEY_CALC(SESSION)
            %
            % Creates a HARTLEY_CALC ndi.calculator object
            %
            hartley_calc_obj = hartley_calc_obj@ndi.calculator(session,'hartley_calc',...
                'hartley_calc');
        end % hartley_calc()

        function doc = calculate(ndi_calculator_obj, parameters)
            % CALCULATE - perform the calculator for ndi.calc.example.hartley_calc
            %
            % DOC = CALCULATE(NDI_CALCULATOR_OBJ, PARAMETERS)
            %
            % Creates a hartley_calc_calc document given input parameters.
            %
            % The document that is created hartley_calc
            % by the input parameters.
            arguments
                ndi_calculator_obj
                parameters (1,1) struct {ndi.validators.mustHaveFields(parameters,{'input_parameters','depends_on'})}
            end

            doc = {};

            % Step 1: set up the output structure, and load the element_id and stimulus_presentation_doc
            hartley_calc = parameters;

            element_doc = ndi_calculator_obj.session.database_search(ndi.query('base.id','exact_string',...
                did.db.struct_name_value_search(parameters.depends_on,'element_id'),''));
            if numel(element_doc)~=1
                error('Could not find element doc..');
            end
            element_doc = element_doc{1};
            element = ndi.database.fun.ndi_document2ndi_object(element_doc, ndi_calculator_obj.session);

            q1 = ndi.query('','isa','stimulus_presentation','');
            stimulus_presentation_docs = ndi_calculator_obj.session.database_search(q1);
            ndi_decoder = ndi.app.stimulus.decoder(ndi_calculator_obj.session);

            for i=1:numel(stimulus_presentation_docs)
                presentation_time = ndi_decoder.load_presentation_time(stimulus_presentation_docs{i});
                [b,stimids] = ndi.calc.vis.hartley_calc.ishartleystim(stimulus_presentation_docs{i});

                if ~b, continue; end

                stimulus_element = ndi.database.fun.ndi_document2ndi_object(dependency_value(stimulus_presentation_docs{i},'stimulus_element_id'),...
                    ndi_calculator_obj.session);

                % Step 2: do we have a stimulus presentation that has Hartley stims in it? Was it running at the same time as our element?

                et = stimulus_element.epochtable();

                % ASSUMPTION: each stimulus epoch will overlap a single element epoch
                stim_timeref = ndi.time.timereference(stimulus_element,...
                    ndi.time.clocktype(presentation_time(1).clocktype),...
                    stimulus_presentation_docs{i}.document_properties.epochid.epochid,...
                    presentation_time(1).onset);
                [ts_epoch_t0_out, ts_epoch_timeref, msg] = ndi_calculator_obj.session.syncgraph.time_convert(stim_timeref,...
                    0, element, ndi.time.clocktype('dev_local_time'));
                % time is 0 because stim_timeref is relative to 1st stim

                if ~isempty(ts_epoch_t0_out)% we have a match

                    % Step 3: now to calculate

                    % Step 3a: set up variables for returning

                            % EDIT THIS SO IT TAKES ALL STIMIDS, MORE THAN JUST THE FIRST HARTLEY STIM
                    hartley_reverse_correlation.stimulus_properties = ndi.calc.vis.hartley_calc.hartleystimdocstruct(...
                        stimulus_presentation_docs{i}.document_properties.stimulus_presentation.stimuli(stimids(1)).parameters);
                    hartley_reverse_correlation.reconstruction_properties = struct(...
                        'T_coords', parameters.input_parameters.T,...
                        'X_coords', 1:parameters.input_parameters.X_sample:hartley_reverse_correlation.stimulus_properties.M,...
                        'Y_coords', 1:parameters.input_parameters.Y_sample:hartley_reverse_correlation.stimulus_properties.M);
                    reverse_correlation.method = 'Hartley';
                    reverse_correlation.dimension_labels = '';

                    % Step 3b: load the spike times and spike parameters

                    frameTimes_indexes = find(presentation_time.stimevents(:,2)==1);
                    frameTimes = presentation_time.stimevents(frameTimes_indexes,1);
                    [spike_values,spike_times] = element.readtimeseries(stim_timeref, ...
                        presentation_time(1).onset, ...
                        presentation_time(end).offset);

                    % load the hartley states
                    P = stimulus_presentation_docs{i}.document_properties.stimulus_presentation.stimuli(stimids(1)).parameters;
                    P.rect = P.rect(:)';
                    P.dispprefs = {};
                    P.chromhigh = P.chromhigh(:)';
                    P.chromlow = P.chromlow(:)';
                    hartleys = hartleystim(P);
                    [S,KXV,KYV,ORDER] = hartleynumbers(hartleys);

                    hartley_reverse_correlation.hartley_numbers.S = S;
                    hartley_reverse_correlation.hartley_numbers.KXV = KXV;
                    hartley_reverse_correlation.hartley_numbers.KYV = KYV;
                    hartley_reverse_correlation.hartley_numbers.ORDER = ORDER;

                    hartley_reverse_correlation.spiketimes = spike_times;
                    hartley_reverse_correlation.frameTimes = frameTimes;

                    % write JSON file here for now

                    mystring = vlt.data.jsonencodenan(hartley_reverse_correlation);
                    mypath = fullfile(ndi_calculator_obj.session.path,'hartley');
                    if ~isfolder(mypath)
                        mkdir(mypath);
                    end
                    myfile = fullfile(mypath,...
                        [stimulus_presentation_docs{i}.document_properties.epochid.epochid '_' ...
                            element.elementstring() '_hartley.json']);
                    mystring = char(vlt.data.prettyjson(mystring));
                    vlt.file.str2text(myfile,mystring);

                    [sta,p_val, rescale, cmap] = vis.revcorr.sta_pipeline(hartley_reverse_correlation.hartley_numbers.S,...
                        hartley_reverse_correlation.hartley_numbers.KXV, ...
                        hartley_reverse_correlation.hartley_numbers.KYV, ...
                        hartley_reverse_correlation.frameTimes, ...
                        hartley_reverse_correlation.spiketimes, ...
                        hartley_reverse_correlation.reconstruction_properties.T_coords, ...
                        hartley_reverse_correlation.reconstruction_properties.X_coords, ...
                        hartley_reverse_correlation.reconstruction_properties.Y_coords);

                    sta = sta(:,:,end:-1:1); % make it match T
                    p_val = p_val(:,:,end:-1:1); % make it match T

                    % Step 3c: actually make the document
                    if 1%
                        ngridp.data_size = 8;
                        ngridp.data_type = 'double';
                        ngridp.data_dim = [size(sta) 2];
                        ngridp.coordinates = [hartley_reverse_correlation.reconstruction_properties.T_coords(:); ...
                            hartley_reverse_correlation.reconstruction_properties.X_coords(:); ...
                            hartley_reverse_correlation.reconstruction_properties.Y_coords(:)];

                        doc{end+1} = ndi.document(ndi_calculator_obj.doc_document_types{1},...
                            'hartley_calc',parameters,...
                            'hartley_reverse_correlation',hartley_reverse_correlation,...
                            'reverse_correlation',reverse_correlation,'ngrid',ngridp) + ndi_calculator_obj.newdocument();
                        doc{end} = doc{end}.set_dependency_value('element_id',element_doc.id());
                        doc{end} = doc{end}.set_dependency_value('stimulus_presentation_id', stimulus_presentation_docs{i}.id());

                        % open the ngrid file
                         % TODO: update when new database available
                         % this is one problem that needs to be fixed, we need to be able to specify binary
                         % data before adding to the database; we will use this temporary workaround of writing
                         % to the hartley directory
                        myfile = fullfile(mypath,[doc{end}.id() '.ngrid']);
                        fid = fopen(myfile,'w','ieee-le');
                        if fid<0
                            error(['Could not open file ' myfile '.']);
                        end
                        % write the ngrid file
                        fwrite(fid,cat(4,sta,p_val),'double');
                        fclose(fid);

                        doc{end} = doc{end}.add_file('hartley_results.ngrid',myfile);
                    end
                end
            end
        end % calculate

        function parameters = default_search_for_input_parameters(ndi_calculator_obj)
            % DEFAULT_SEARCH_FOR_INPUT_PARAMETERS - default parameters for searching for inputs
            %
            % PARAMETERS = DEFAULT_SEARCH_FOR_INPUT_PARAMETERS(NDI_CALCULATOR_OBJ)
            %
            % Returns a list of the default search parameters for finding appropriate inputs
            % to the calculator. For hartley_calc_calc, there is no appropriate default parameters
            % so this search will yield empty.
            %
            parameters.input_parameters = struct(...
                'T', -0.100:0.010:0.250, ...
                'X_sample', 1, ...
                'Y_sample', 1);
            parameters.depends_on = did.datastructures.emptystruct('name','value');
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
            % For the ndi.calc.stimulus.hartley_calc_calc class, this looks for
            % documents of type 'stimulus_response_scalar.json' with 'response_type' fields
            % the contain 'mean' or 'F1'.
            %
            %
            q_total = ndi.query('element.type','exact_string','spikes','');

            query = struct('name','element_id','query',q_total);
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

            % call superclass plot method to set up axes
            h=plot@ndi.calculator(ndi_calculator_obj, doc_or_parameters, varargin{:});

            if isa(doc_or_parameters,'ndi.document')
                doc = doc_or_parameters;
            else
                error('Do not know how to proceed without an ndi document for doc_or_parameters.');
            end

            plotrows = 6;
            plotcols = 6;

            [sta,pval] = read_sta(ndi_calculator_obj, doc);
            clim = [-1 1] * max(abs([min(sta(:)) max(sta(:))]));
            significance_plot = vis.revcorr.rescale_p_image(pval);
            cmap = vis.revcorr.get_cmap();

            size_param = min(36,size(sta,3));
            H = vlt.image.stack2tile(sta(:,:,1:size_param),plotrows,plotcols);
            Hp = vlt.image.stack2tile(significance_plot(:,:,1:size_param),plotrows,plotcols);

            ax(1) = subplot(2,2,1);
            imshow(H{1},clim);

            ax(2)=subplot(2,2,3);
            image(Hp{1});
            colormap(ax(2),cmap);
            axis equal off;

            linkaxes([ax(1) ax(2)]);

            rp = doc.document_properties.hartley_reverse_correlation.reconstruction_properties;
            [theta,sta_r,pval_r] = vis.revcorr.rotate_sta(rp.T_coords, sta, significance_plot);
            [t_profile,t_profile_pval] = vis.revcorr.peak_time_profile(rp.T_coords,sta,significance_plot);

            ax(3) = subplot(2,2,2);
            imshow(t_profile,clim,'XData',rp.T_coords,'YData',rp.Y_coords);
            axis normal on;
            grid on;

            ax(4) = subplot(2,2,4);
            image(t_profile_pval,'XData',rp.T_coords,'YData',rp.Y_coords);
            colormap(ax(4),cmap);
            axis normal on;
            grid on;

            linkaxes([ax(3) ax(4)]);

            h.objects = cat(1,h.objects,ax);

        end % plot()

        function [sta,pval] = read_sta(ndi_calculator_obj, doc_or_id)
            % READ_STA - read the spike-triggered-average file from disk
            %
            % [STA, PVAL] = READ_STA(NDI_CALCULATOR_OBJ, DOC_OR_ID)
            %
            % Reads the spike-triggered average and p-values from disk.
            % NDI_CALCULATOR_OBJ should be an ndi.calc.vis.hartley_calc object and
            % DOC_OR_ID should either be an ndi.document or the id of the hartley
            % document object to be read.
            %

            if ischar(doc_or_id)
                doc = ndi_calculator_obj.database_search(ndi.query('base.id','exact_string',doc_or_id,''));
                if numel(doc)~=1
                    % there cannot be two documents with same id
                    error(['No document with id ' doc_or_id ' found.']);
                end
                doc = doc{1};
            else
                doc = doc_or_id;
            end

            if 0
                mypath = fullfile(ndi_calculator_obj.session.path,'hartley');

                myfile = fullfile(mypath,[doc.id() '.ngrid']);
                fid = fopen(myfile,'r','ieee-le');
                if fid<0
                    error(['Could not open file ' myfile '.']);
                end
                read the ngrid file

            end

            myfile = ndi_calculator_obj.session.database_openbinarydoc(doc,'hartley_results.ngrid');
            fulldata = fread(myfile,prod(doc.document_properties.ngrid.data_dim),doc.document_properties.ngrid.data_type);
            ndi_calculator_obj.session.database_closebinarydoc(myfile);

            fulldata = reshape(fulldata,vlt.data.rowvec(doc.document_properties.ngrid.data_dim));
            sta = fulldata(:,:,:,1);
            pval = fulldata(:,:,:,2);
        end % read_sta()

    end % methods()

    methods (Static)
        function [b,stimids] = ishartleystim(stim_presentation_doc)
            % ISHARTLEYSTIM - does a stimulus presentation doc contain a Hartley stimulus?
            %
            % [B,STIMIDS] = ndi.calc.hartley_calc.ishartleystim(STIM_PRESENTATION_DOC)
            %
            % Returns 1 iff STIM_PRESENTATION_DOC contains Hartley stimuli. Returns
            % the STIMIDS of any Hartley stimuli.
            %
            stimids = [];
            S = stim_presentation_doc.document_properties.stimulus_presentation.stimuli;
            for i=1:numel(S)
                thestruct = S(i).parameters;
                if isfield(thestruct,'M')&isfield(thestruct,'chromhigh')&isfield(thestruct,'K_absmax')% hartley
                    stimids(end+1) = i;
                end
            end
            b = ~isempty(stimids);
        end % ishartleystim

        function hartleydocinfo = hartleystimdocstruct(stimstruct)
            % HARTLEYSTIMDOCSTRUCT - return the fields of the Hartley stimulus necessary for the hartley_reverse_correlation document
            %
            % HARTLEYDOCINFO = HARTLEYSTIMDOCSTRUCT(STIMSTRUCT)
            %
            % Returns the fields of the Hartley stim that are needed for the
            % NDI hartley_reverse_correlation document:
            %
            % Fields: M, L_max, K_max, sf_max, fps, color_high, color_low, rect
            %
            %
            fields_out = {'M','L_max','K_max','sf_max','fps','color_high','color_low','rect'};
            fields_names = {'M', 'L_absmax','K_absmax','sfmax','fps','chromhigh','chromlow','rect'};
            hartleydocinfo = did.datastructures.emptystruct(fields_out{:});
            hartleydocinfo(1).M = stimstruct.M;
            for i=1:numel(fields_out)
                if ~isfield(stimstruct,fields_names{i})
                    error(['STIMSTRUCT has no field ' fields_names{i} '.']);
                end
                hartleydocinfo = setfield(hartleydocinfo, fields_out{i}, getfield(stimstruct,fields_names{i}));
            end
        end % hartleystimdocinfo

    end % static methods

end % hartley_calc
