classdef hartley < ndi.calculator

    methods
        function hartley_obj = hartley(session)
            % HARTLEY - a hartley demonstration of an ndi.calculator object
            %
            % HARTLEY_OBJ = HARTLEY(SESSION)
            %
            % Creates a HARTLEY ndi.calculator object
            %
            hartley_obj = hartley_obj@ndi.calculator(session,'hartley',...
                'hartley_calc');
            hartley_obj.numberOfSelfTests = 1;
        end % hartley()

        function [docs, doc_output, doc_expected_output] = generate_mock_docs(obj, scope, number_of_tests, kwargs)
            % GENERATE_MOCK_DOCS - generate mock documents for testing
            %
            % [DOCS, DOC_OUTPUT, DOC_EXPECTED_OUTPUT] = GENERATE_MOCK_DOCS(OBJ, SCOPE, NUMBER_OF_TESTS, ...)
            %
            % Generates mock documents and runs the calculator for testing purposes.
            %
            % Inputs:
            %   OBJ - the calculator object
            %   SCOPE - 'highSNR' or 'lowSNR'
            %   NUMBER_OF_TESTS - number of tests to run (only 1 supported for now)
            %   KWARGS - keyword arguments:
            %     'generate_expected_docs' (boolean, default false) - if true, save expected output
            %     'specific_test_inds' (vector, default []) - run only specific test indices
            %

            arguments
                obj
                scope (1,:) char {mustBeMember(scope,{'highSNR','lowSNR'})}
                number_of_tests
                kwargs.generate_expected_docs (1,1) logical = false
                kwargs.specific_test_inds (1,:) double = []
            end

            docs = cell(1,number_of_tests);
            doc_output = cell(1,number_of_tests);
            doc_expected_output = cell(1,number_of_tests);

            p = fileparts(mfilename('fullpath'));
            mock_dir = fullfile(p, 'mock', 'hartley');
            if ~isfolder(mock_dir)
                mkdir(mock_dir);
            end

            for i = 1:number_of_tests
                if ~isempty(kwargs.specific_test_inds) && ~ismember(i, kwargs.specific_test_inds)
                    continue;
                end

                % Use Object Session
                S = obj.session;

                % 1. Create mock elements (Subject, Stimulator, Spikes)
                mock_data = ndi.mock.fun.subject_stimulator_neuron(S);
                nde = mock_data.mock_spikes;
                nde_stim = mock_data.mock_stimulator;

                current_docs = {mock_data.mock_subject, nde_stim.load_element_doc(), nde.load_element_doc()};

                % 2. Calculate Hartley Response (imitating vis.revcorr.test)

                % Construct P fully inline to ensure reproducibility and remove external dependencies
                P = struct();
                P.M = 200;
                P.L_absmax = 20;
                P.K_absmax = 20;
                P.sfmax = 3;
                P.fps = 10;
                P.chromhigh = [255 255 255];
                P.chromlow = [0 0 0];
                P.rect = [0 -100 800 700];
                P.reps = 5;
                P.windowShape = 0;
                P.distance = 30;
                P.contrast = 1;
                P.background = 0.5;
                P.backdrop = 0.5;
                P.dispprefs = {};

                % Use a deterministic random state for reproducibility
                rng(1);
                P.randState = rand(35,1);

                % Generate the consistent stimulus sequence
                hartleys = hartleystim(P);
                [s, kx_v, ky_v, ORDER] = hartleynumbers(hartleys);

                % Generate frame times
                % s has one entry per frame.
                num_frames = numel(s);
                frameTimes = (0:num_frames-1)' / P.fps;

                % RF parameters
                rfTimeRange = 0.2;
                max_TimeBlock_StartTime = 500;
                rfDeltaT = 0.005;
                rfNumTimeSteps = 40;
                responseDeltaT = 0.01;
                threshold = 1;
                Verbose = true;

                [rf,rfTimeLags] = vis.revcorr.setRF(P.M, rfNumTimeSteps, rfDeltaT);

                % Calculate Spike Times
                [response, spiketimes] = vis.revcorr.calculateHartleyResponse(s, kx_v, ky_v, frameTimes, rf, ...
                    'rfDeltaT', rfDeltaT, 'rfNumTimeSteps', rfNumTimeSteps, 'responseDeltaT', responseDeltaT, ...
                    'max_TimeBlock_StartTime', max_TimeBlock_StartTime, 'threshold', threshold, 'rfTimeRange', rfTimeRange, 'Verbose', Verbose);

                % 3. Update the spikes element and stimulator with the epoch

                % Define epoch
                epoch_name = 'mockepoch';
                clock_type = 'dev_local_time,utc';

                % Time range
                t_start = frameTimes(1) - 1;
                t_end = frameTimes(end) + 5;
                t0_t1 = [t_start t_end; t_start t_end];

                nde.addepoch(epoch_name, clock_type, t0_t1, spiketimes(:), ones(size(spiketimes(:))));
                nde_stim.addepoch(epoch_name, clock_type, t0_t1, [], []);

                % 4. Generate Stimulus Presentation Document

                % Create stimulus presentation structure
                stim_pres_struct = struct();
                stim_pres_struct.presentation_order = 1; % One big stimulus
                % The calculator expects one stimulus entry that describes the whole Hartley sequence
                stim_pres_struct.stimuli(1).parameters = P;

                % Create presentation_time structure
                pt_struct = vlt.data.emptystruct('clocktype','stimopen','onset','offset','stimclose','stimevents');
                pt_struct(1).clocktype = 'dev_local_time';
                pt_struct(1).onset = frameTimes(1);
                pt_struct(1).offset = frameTimes(end);
                pt_struct(1).stimopen = frameTimes(1);
                pt_struct(1).stimclose = frameTimes(end);
                pt_struct(1).stimevents = [frameTimes(:), ones(size(frameTimes(:)))];

                % Write presentation time to binary file
                pt_filename = ndi.file.temp_name();
                ndi.database.fun.write_presentation_time_structure(pt_filename, pt_struct);

                epochid_struct.epochid = epoch_name;

                stim_pres_doc = ndi.document('stimulus_presentation', 'stimulus_presentation', stim_pres_struct, 'epochid', epochid_struct) + ...
                    obj.session.newdocument();
                stim_pres_doc = stim_pres_doc.set_dependency_value('stimulus_element_id', nde_stim.id());
                stim_pres_doc = stim_pres_doc.add_file('presentation_time.bin', pt_filename);

                S.database_add(stim_pres_doc);
                current_docs{end+1} = stim_pres_doc;

                docs{i} = current_docs;

                % 6. Run the calculator
                hc = ndi.calc.vis.hartley(S);
                parameters = hc.default_search_for_input_parameters();
                parameters.query.query = ndi.query('element.type','exact_string','spikes','');

                hc.run('Replace', parameters);

                % Get Output
                q = ndi.query('','isa','hartley_calc','');
                out_docs = S.database_search(q);

                if ~isempty(out_docs)
                    doc_output{i} = out_docs{1};
                else
                    error('No output document generated.');
                end

                % Handle Expected Docs
                if kwargs.generate_expected_docs
                    doc_expected_output{i} = doc_output{i};
                    vlt.file.str2text(fullfile(mock_dir, ['mock.' int2str(i) '.json']), jsonencode(doc_expected_output{i}.document_properties));

                    % Also copy the binary results file
                    % Try to locate it in `obj.session.path/hartley`.
                    generated_file = fullfile(obj.session.path, 'hartley', [doc_output{i}.id() '.ngrid']);
                    if isfile(generated_file)
                        expected_binary_name = ['mock_' int2str(i) '_hartley_results.ngrid'];
                        expected_binary_path = fullfile(mock_dir, expected_binary_name);
                        copyfile(generated_file, expected_binary_path);
                    end
                else
                    expected_file = fullfile(mock_dir, ['mock.' int2str(i) '.json']);
                    if isfile(expected_file)
                        doc_expected_output{i} = ndi.document(jsondecode(fileread(expected_file)));
                    else
                        doc_expected_output{i} = [];
                    end
                end
            end
        end

        function [b, report] = compare_mock_docs(ndi_calculator_obj, expected_doc, actual_doc, scope, docCompare)
            % COMPARE_MOCK_DOCS - compare expected and actual documents
            %
            % [B, REPORT] = COMPARE_MOCK_DOCS(NDI_CALCULATOR_OBJ, EXPECTED_DOC, ACTUAL_DOC, SCOPE, DOCCOMPARE)
            %
            % Compares the EXPECTED_DOC and ACTUAL_DOC.

            % Call parent method first
            [b, report] = compare_mock_docs@ndi.calculator(ndi_calculator_obj, expected_doc, actual_doc, scope, docCompare);

            if ~b
                return;
            end

            % Compare binary data
            p = fileparts(mfilename('fullpath'));
            mock_dir = fullfile(p, 'mock', 'hartley');

            % Assuming test index 1 for now.
            % Ideally we should infer index from expected_doc filename if possible,
            % but standard ndi.mock framework doesn't explicitly pass it.
            % Given user instructions for this task, checking 'mock_1_hartley_results.ngrid' is the target.
            expected_binary_path = fullfile(mock_dir, 'mock_1_hartley_results.ngrid');

            if ~isfile(expected_binary_path)
                % If we can't find specific file, warn but don't fail basic doc check
                warning(['Expected binary file ' expected_binary_path ' not found. Skipping binary comparison.']);
                return;
            end

            % Read expected data
            fid_exp = fopen(expected_binary_path, 'r', 'ieee-le');
            if fid_exp < 0
                error(['Could not open expected binary file ' expected_binary_path]);
            end
            % We need dimensions. Expected doc should have them.
            dims_exp = expected_doc.document_properties.ngrid.data_dim;
            type_exp = expected_doc.document_properties.ngrid.data_type;

            [sta_exp, ~] = ndi.calc.vis.hartley.readStaFromFile(fid_exp, dims_exp, type_exp);
            fclose(fid_exp);

            % Read actual data
            [sta_act, ~] = ndi_calculator_obj.read_sta(actual_doc);

            % Compare 10th frame (3rd dimension)
            if size(sta_exp, 3) < 10 || size(sta_act, 3) < 10
                b = 0;
                report = 'STA data has fewer than 10 frames.';
                return;
            end

            frame_exp = sta_exp(:, :, 10);
            frame_act = sta_act(:, :, 10);

            diff = abs(frame_exp - frame_act);
            max_diff = max(diff(:));

            if max_diff > 1e-2
                b = 0;
                report = sprintf('STA data mismatch at frame 10. Max diff: %g', max_diff);
            else
                b = 1;
                report = ''; % Success
            end
        end

        function doc = calculate(ndi_calculator_obj, parameters)
            % CALCULATE - perform the calculator for ndi.calc.vis.hartley
            %
            % DOC = CALCULATE(NDI_CALCULATOR_OBJ, PARAMETERS)
            %
            % Creates a hartley_calc document given input parameters.
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
                [b,stimids] = ndi.calc.vis.hartley.ishartleystim(stimulus_presentation_docs{i});

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
                    hartley_reverse_correlation.stimulus_properties = ndi.calc.vis.hartley.hartleystimdocstruct(...
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
            % to the calculator. For hartley_calc, there is no appropriate default parameters
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
            % For the ndi.calc.stimulus.hartley_calc class, this looks for
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
            % NDI_CALCULATOR_OBJ should be an ndi.calc.vis.hartley object and
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
            [sta,pval] = ndi.calc.vis.hartley.readStaFromFile(myfile, doc.document_properties.ngrid.data_dim, doc.document_properties.ngrid.data_type);
            ndi_calculator_obj.session.database_closebinarydoc(myfile);
        end % read_sta()

    end % methods()

    methods (Static)
        function [sta,pval] = readStaFromFile(fileid_or_obj, data_dim, data_type)
            % READSTAFROMFILE - read STA and PVAL from a file identifier or object
            %
            % [STA, PVAL] = READSTAFROMFILE(FILEID_OR_OBJ, DATA_DIM, DATA_TYPE)
            %
            % Reads the spike-triggered average (STA) and p-value (PVAL) data from
            % a file identifier or object.
            %
            % Inputs:
            %   FILEID_OR_OBJ - A file identifier (integer) or a file object (e.g.,
            %                   did.file.fileobj or vlt.file.fileobj) that supports fread.
            %   DATA_DIM - The dimensions of the data to be read.
            %   DATA_TYPE - The precision of the data (e.g., 'double', 'single').
            %
            % Outputs:
            %   STA - The spike-triggered average data.
            %   PVAL - The p-values associated with the STA.
            %

            fulldata = fread(fileid_or_obj, prod(data_dim), data_type);

            fulldata = reshape(fulldata, vlt.data.rowvec(data_dim));
            sta = fulldata(:,:,:,1);
            pval = fulldata(:,:,:,2);
        end % readStaFromFile

        function [b,stimids] = ishartleystim(stim_presentation_doc)
            % ISHARTLEYSTIM - does a stimulus presentation doc contain a Hartley stimulus?
            %
            % [B,STIMIDS] = ndi.calc.vis.hartley.ishartleystim(STIM_PRESENTATION_DOC)
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

end % hartley
