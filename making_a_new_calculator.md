# Today I made a new calculator. This is what I did.

My calculator is going to be named 'spike_shape'. It computes the mean spike shape for a neuron_extracellular object in all of its recorded epochs, and at regular intervals within epochs.

# 1. Thinking before coding

The first thing I did was to think carefully about its input and output.

I decided I needed the calculator to generate an output document for each epoch. Why not a single document with all the epoch data in a single place, so we don't have a massive proliferation of documents? Well, what if we record more data later? If that happened, we would have to delete the document and recompute it, because there is NO UPDATING of values in NDI databases (only deletion and adding a new document).

OK, the output will be a document that records the epoch examined, and records the mean spike shape over all channels and standard deviation at a) the beginning of the record, b) the end of the record, and c) every n seconds during the record.

I am doing this while we still haven't updated to the database version with multiple binary files, so I am going to put them all into a single file.

For the binary file format, I am going to choose `vhlspikewaveformfile` from our vhlab-toolbox-matlab library (see [vlt.file.custom_file_formats.readvhlspikewaveformfile](https://github.com/VH-Lab/vhlab-toolbox-matlab/blob/master/%2Bvlt/%2Bfile/%2Bcustom_file_formats/readvhlspikewaveformfile.m), for example). It is made to store spike waveforms and has nice reader/writer functions that we wrote previously.

# 2. Setting things up

I know I need a `spike_shape_calc` document type and a `spike_shape.m` file with the code.

I am going to design in NDIcalc-vis-matlab. Really, we ought to make a new NDIcalc-extracellular-matlab package because this has applications well beyond vision. I will pledge to do it but I have a grant deadline right now.

I make copies of the `simple` calculator so I have templates to start from.  

1. I copy `[...]/NDI-matlab/+ndi/+calc/+example/simple.m` to `[...]/NDIcalc-vis-matlab/+ndi/+calc/+extracellular/spike_shape.m`
2. I copy `[...]/NDI-matlab/+ndi/+calc/+example/docs/simple.docs*.txt .m` to a new directory that I made `[...]/NDIcalc-vis-matlab/+ndi/+calc/+extracellular/docs`. Then I rename them all so that instead of simple.docs._something_ it says spike_shape.docs._something_.
3. I copy `[...]/NDI-matlab/ndi_common/database_documents/apps/calculators/simple_calc.json` to `[...]/NDIcalc-vis-matlab/ndi_common/database_documents/calc/spike_shape_calc.json` 

# 3. Designing the output document type

I want to start with the output document. I want to have fields for selecting how much time around each spike should be selected (I'll call `spike_window_time_before`, `spike_window_time_after`), and the time interval over which spikes should be averaged (`averaging_interval`). The spike waves will be stored in our binary file of type .vsw, but we want to store some metadata to tell the user the center time of each time bin (`interval_center_times` and the number of spikes that were found in each time bin `number_of_spikes_per_interval`.

```json
{
        "document_class": {
                "definition":                                           "$NDICALCDOCUMENTPATH\/calc\/spike_shape_calc.json",
                "validation":                                           "$NDICALCSCHEMAPATH\/calc\/spike_shape_calc_shema.json",
                "class_name":                                           "spike_shape_calc",
                "property_list_name":                                   "spike_shape_calc",
                "class_version":                                        1,
                "superclasses": [
                        { "definition":                                 "$NDIDOCUMENTPATH\/ndi_document.json" },
                        { "definition":                                 "$NDIDOCUMENTPATH\/ndi_document_app.json" }
                ]
        },
        "depends_on": [
                {       "name": "element_epoch_id",
                        "value": 0
                },
                {       "name": "element_id",
                        "value": 0,
                }
        ],
        "spike_shape_calc": {
                "input_parameters": {
                        "spike_window_before_time":                                     -0.001,
                        "spike_window_after_time":                                      0.002,
                        "averaging_interval":                                           60
                }
                "interval_center_times":                                                [  ]
                "number_of_spikes_per_interval":                                        [ ]
        }
}
```

# 4. Writing the code

1. First, I start with my template, and I just do a find and replace on "simple" to change it to "spike_shape".
2. Next, I update the creator with some documentation and add code that adds the database document type to the object by calling the superclass creator:
```matlab
                function spike_shape_obj = spike_shape(session)
                        % SPIKE_SHAPE_CALC - calculator that produces spike waveform shapes in an epoch 
                        %
                        % SPIKE_SHAPE_CALC_OBJ = SPIKE_SHAPE_CALC(SESSION)
                        %
                        % Creates a SPIKE_SHAPE_CALC ndi.calculator object
                        %
                                w = which('ndi.calc.vis.spike_shape');
                                parparparpar = fileparts(fileparts(fileparts(fileparts(w))));
                                spike_shape_obj = spike_shape_obj@ndi.calculator(session,'spike_shape_calc',...
                                        fullfile(parparparpar,'ndi_common','database_documents','calc','spike_shape_calc.json'));
                end; % spike_shape()
```
3. The next function I edit is the function that searches for input parameters. Here I want it to find all spiking neurons as input:
```matlab
                function parameters = default_search_for_input_parameters(ndi_calculator_obj)
                        % DEFAULT_SEARCH_FOR_INPUT_PARAMETERS - default parameters for searching for inputs
                        %
                        % PARAMETERS = DEFAULT_SEARCH_FOR_INPUT_PARAMETERS(NDI_CALCULATOR_OBJ)
                        %
                        % Returns a list of the default search parameters for finding appropriate inputs
                        % to the calculator.
                        %
                                parameters.input_parameters = struct('spike_window_before_time',-0.001,...
                                        'spike_window_after_time', 0.002, 'averaging_interval', 60);
                                parameters.depends_on = vlt.data.emptystruct('name','value');
                                parameters.query = struct('name','element_epoch_id','query',ndi.query('element.type','exact_string','spikes',''));
                end; % default_search_for_input_parameters
```





