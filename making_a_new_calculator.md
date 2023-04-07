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

I want to start with the output document. I want to have fields for selecting how much time around each spike should be selected (I'll call spike_window_time_before, spike_window_time_after), and the time interval over which spikes should be averaged (averaging_interval).

