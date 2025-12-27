classdef spike_shape < ndi.calculator

    methods

        function spike_shape_obj = spike_shape(session)
            % SPIKE_SHAPE_CALC - calculator that produces spike waveform shapes in an epoch
            %
            % SPIKE_SHAPE_CALC_OBJ = SPIKE_SHAPE_CALC(SESSION)
            %
            % Creates a SPIKE_SHAPE_CALC ndi.calculator object
            %
            spike_shape_obj = spike_shape_obj@ndi.calculator(session,'spike_shape_calc',...
                                    'spike_shape_calc');
        end % spike_shape()

        function doc = calculate(ndi_calculator_obj, parameters)
            % CALCULATE - perform the calculator for ndi.calc.example.spike_shape
            %
            % DOC = CALCULATE(NDI_CALCULATOR_OBJ, PARAMETERS)
            %
            % Creates a spike_shape_calc document given input parameters.
            %
            % The document that is created spike_shape has an 'answer' that is given
            % by the input parameters.
            arguments
                ndi_calculator_obj
                parameters (1,1) struct {ndi.validators.mustHaveFields(parameters,{'input_parameters','depends_on'})}
            end

            % Step 1: set up
            element_id = did.db.struct_name_value_search(parameters.depends_on,'element_id');
            spiking_element = ndi.database.fun.ndi_document2ndi_object(element_id, ndi_calculator_obj.session);

            % Step 2: perform the calculator for each recording epoch of the element

            dirpath = [ndi_calculator_obj.session.path filesep 'ndiobjects'];
            if ~isfolder(dirpath)
                mkdir(dirpath);
            end

            doc = {};

            et = spiking_element.epochtable(); % find all the epochs for this element
            for i=1:numel(et)
                disp(['Now processing epoch ' int2str(i) ' of ' int2str(numel(et)) '.']);
                spike_shape = parameters;
                q1 = ndi.query('','depends_on','element_id',element_id);
                q2 = ndi.query('epochid.epochid','exact_string',et(i).epoch_id,'');
                epoch_id_doc = ndi_calculator_obj.session.database_search(q1&q2);
                if numel(epoch_id_doc)~=1
                    error(['Could not find exactly 1 epoch id doc for ' et(i).epoch_id '.']);
                end
                epoch_id_doc = epoch_id_doc{1};

                [mean_waves,std_waves,output_parameters] = ndi.fun.spiketrains.mean_spike_waveforms(spiking_element, et(i).epoch_id,...
                    'spike_window_before_time',parameters.input_parameters.spike_window_before_time,...
                    'spike_window_beforeafter_time',parameters.input_parameters.spike_window_after_time,...
                    'averaging_window', parameters.input_parameters.averaging_window,...
                    'averaging_window_step', parameters.input_parameters.averaging_window_step,...
                    'filter_padding', parameters.input_parameters.filter_padding',...
                    'cheby_order', parameters.input_parameters.cheby_order,...
                    'cheby_R', parameters.input_parameters.cheby_R,...
                    'cheby_cutoff', parameters.input_parameters.cheby_cutoff');

                spike_shape.interval_center_times = output_parameters.interval_center_times;
                spike_shape.number_of_spikes_per_interval = output_parameters.number_of_spikes_per_interval;
                spike_shape.sample_times = output_parameters.sample_times;

                % Step 3: place the results of the calculator into an NDI document
                doc{end+1} = ndi.document(ndi_calculator_obj.doc_document_types{1},'spike_shape_calc',spike_shape);
                for k=1:numel(parameters.depends_on)
                    doc{end} = doc{end}.set_dependency_value(parameters.depends_on(k).name,parameters.depends_on(k).value);
                end
                doc{end} = doc{end}.set_dependency_value('element_epoch_id', epoch_id_doc.document_properties.base.id);

                % files
                fileparameters.numchannels = size(mean_waves,2);
                fileparameters.S0 = output_parameters.s0;
                fileparameters.S1 = output_parameters.s1;
                fileparameters.name = '';
                fileparameters.ref = 0;
                fileparameters.comment = '';
                fileparameters.samplingrate = output_parameters.sample_rate;

                fname = [dirpath filesep doc{end}.document_properties.base.id '.vsw'];
                file_obj = vlt.file.fileobj('fullpathfilename',fname,'permission','w','machineformat','ieee-le');
                file_obj = file_obj.fopen();

                fid = vlt.file.custom_file_formats.newvhlspikewaveformfile(file_obj,fileparameters);
                for j=1:size(mean_waves,3)
                    vlt.file.custom_file_formats.addvhlspikewaveformfile(file_obj, mean_waves(:,:,j));
                    vlt.file.custom_file_formats.addvhlspikewaveformfile(file_obj, std_waves(:,:,j));
                end
                file_obj.fclose();

                doc{end} = doc{end}.add_file('spikewaves.vsw',fname);
            end

        end % calculate

        function parameters = default_search_for_input_parameters(ndi_calculator_obj)
            % DEFAULT_SEARCH_FOR_INPUT_PARAMETERS - default parameters for searching for inputs
            %
            % PARAMETERS = DEFAULT_SEARCH_FOR_INPUT_PARAMETERS(NDI_CALCULATOR_OBJ)
            %
            % Returns a list of the default search parameters for finding appropriate inputs
            % to the calculator.
            %
            parameters.input_parameters = struct('spike_window_before_time',-0.001,...
                'spike_window_after_time', 0.002, ...
                'averaging_window', 60, ...
                'averaging_window_step', 300, ...
                'filter_padding', 0.100, ...
                'cheby_order', 4,...
                'cheby_R', 0.5,...
                'cheby_cutoff', 300);
            parameters.depends_on = did.datastructures.emptystruct('name','value');
            parameters.query = struct('name','element_id','query',ndi.query('element.type','exact_string','spikes',''));
        end % default_search_for_input_parameters

        function [mean_waves,std_waves,sample_times] = load(ndi_calculator_obj, doc)
            % LOAD - load binary data from spike_shape_calc document
            %
            % [MEAN_WAVES, STD_WAVES, SAMPLE_TIMES] = LOAD(NDI_CALCULATOR_OBJ, DOC_ID)
            %
            % Loads the mean waveforms MEAN_WAVES and standard deviation waveforms STD_WAVES
            % from document with id DOC_ID.
            %
            % MEAN_WAVES has the form MxCxT, where M is the number of samples per mean spike waveform,
            %   C is the number of channels, and T is the time measurement.
            % STD_WAVES has the form MxCxT, where M is the number of samples per standard deviation spike waveform,
            %   C is the number of channels, and T is the time measurement.
            % SAMPLES_TIMES is an Mx1 vector with the sample times of each spike waveform.
            %
            %dirpath = [ndi_calculator_obj.session.path filesep 'ndiobjects'];
            %fname = [dirpath filesep doc.document_properties.base.id '.vsw'];

            if ~isa(doc,'ndi.document')
                doc = ndi_calculator_obj.session.database_search(ndi.query('base.id','exact_string',doc));
            end

            myfile = ndi_calculator_obj.session.database_openbinarydoc(doc, 'spikewaves.vsw');

            [waveforms, header] = vlt.file.custom_file_formats.readvhlspikewaveformfile(myfile);
            mean_waves = NaN(size(waveforms,1),size(waveforms,2),size(waveforms,3)/2);
            std_waves = NaN(size(waveforms,1),size(waveforms,2),size(waveforms,3)/2);
            for i=1:size(waveforms,3)/2
                mean_waves(:,:,i) = waveforms(:,:,1+2*(i-1));
                std_waves(:,:,i) = waveforms(:,:,2+2*(i-1));
            end
            sample_times = (header.S0:header.S1)/header.samplingrate;

            ndi_calculator_obj.session.database_closebinarydoc(myfile);
        end  % load()

        function h = plot(ndi_calculator_obj, doc_or_parameters, varargin)
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

            [mean_waves,std_waves,sample_times] = ndi_calculator_obj.load(doc);

            delta = 1.3*(sample_times(end)-sample_times(1));

            for i=1:size(mean_waves,3)
                hnew = ndi.fun.plot.multichan(mean_waves(:,:,i),(i-1)*delta+sample_times,30);
                h.objects = cat(1,h.objects,hnew);
            end
        end % plot()

    end % methods()

end % spike_shape
