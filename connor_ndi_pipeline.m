function doc = connor_ndi_pipeline(S, varargin)
% SOPHIE_NDI_PIPELINE - pipeline for Sophie's contrast sensitivity analysis
%
% SOPHIE_NDI_PIPELINE(S, ...)
% Performs Sophie's analysis pipeline on ndi.session S
%
% Accepts name/value pairs that turn on/off parts of the pipeline
% 
% -----------------------------------------------------------
% Parameter (default)      | Description
% -----------------------------------------------------------
% VERBOSE (1)              | 0/1 tell what is going on
% STIM_DECODE (1)          | 0/1 run stimulus decoding?
% STIM_DECODE_FORCE (0)    | 0/1 re-run stimulus decoding?
% STIM_RESPONSE (1)        | 0/1 run stimulus response calculations?
% STIM_RESPONSE_FORCE(0)   | 0/1 re-run stimulus response calculations?
% CONTRAST_TUNING_CURVES   | 0/1 Segment big stimulus runs into smaller
%  (1)                     |       tuning curves?
% CONT_TC_FORCE (0)        | 0/1 Force re-run of contrast tuning curves?
% CONTRAST_FITS (1)        | 0/1 Perform contrast tuning curve fits?
% CONTRAST_FITS_FORCE(0)   | 0/1 Redo contrast tuning curve fits?
% CSENSITIVITY_GROUPING(1) | 0/1 Group analyses of contrast sensitivity? 
% ""_FORCE (0)             | 0/1 Force re-run of group analyses of contrast
%                          |       sensitivity?


VERBOSE = 1;
STIM_DECODE = 1;
STIM_DECODE_FORCE = 1;
STIM_RESPONSE = 1;
STIM_RESPONSE_FORCE = 1;
% CONTRAST_TUNING_CURVES = 0;
% CONT_TC_FORCE = 0;
% CONTRAST_FITS = 0;
% CONTRAST_FITS_FORCE = 0;
% CSENSITIVITY_GROUPING = 0;
% CSENSITIVITY_GROUPING_FORCE = 0;
% ORI_TUNING_CURVES = 1;
% OT_TC_FORCE = 0;
% ORIENTATION_FITS = 1;
% ORIENTATION_FITS_FORCE = 0;
DIR_TC_FORCE = 1;
DIR_TUNING_CURVES = 0;
DIR_FITS = 0;
DIR_FITS_FORCE = 1;
BOOTSTRAPPING = 1;
BOOTSTRAP_FORCE = 1;
assign(varargin{:}); % set the user's preferences

workspace2struct,

%Similar to the NDI ttorial, we start out our pipeline by extracting
%probes and stimuli

if STIM_DECODE

    %create 
    rapp = ndi.app.stimulus.tuning_response(S);
    if VERBOSE
        disp('STEP STIM_DECODE: About to decode stimuli.');
    end
        
    decoder = ndi.app.stimulus.decoder(S);
    p = S.getprobes('type','stimulator');
    if isempty(p)
        error('Could not find stimulator.') 
    end
    decoder.parse_stimuli(p{1},STIM_DECODE_FORCE);
    if VERBOSE
        disp(['STEP STIM_DECODE: Done decoding stimuli.']);
    end
    if VERBOSE,
        disp(['STEP STIM_DECODE: About to label control stimuli.']);
    end;
    cs_doc = rapp.label_control_stimuli(p{1},STIM_DECODE_FORCE);
    if VERBOSE,
        disp(['STEP STIM_DECODE: Done labeling control stimuli.']);
    end;
end;



if STIM_RESPONSE,
    rapp = ndi.app.stimulus.tuning_response(S);
    if VERBOSE,
        disp(['STEP STIM_RESPONSES: About to calculate stimulus responses.']);
    end;
    p = S.getprobes('type','stimulator');
    if isempty(p), error(['Could not find stimulator.']); end;
    stimprobe = p{1};

    if STIM_RESPONSE_FORCE,
        %delete all existing documents
        d = S.database_search(ndi.query('','isa','stimulus_response',''));
        S.database_rm(d);
    end;

    if VERBOSE,
        disp(['STEP STIM_RESPONSES: About to search for spiking neurons...']);
    end;
    e = S.getelements('element.type','spikes');
    if VERBOSE,
        disp(['STEP STIM_RESPONSES: found ' int2str(numel(e)) ' spiking elements, now calculating responses']);
    end;
    for j=1:numel(e),
        if VERBOSE,
            disp(['STEP STIM_RESPONSES: element (' int2str(j) '/' int2str(numel(e)) '), ' e{j}.elementstring() '...']);
        end;
        rapp.stimulus_responses(stimprobe, e{j}, STIM_RESPONSE_FORCE);
    end;    
end;

% if CONTRAST_TUNING_CURVES,
%     c = ndi.calc.stimulus.tuningcurve(S);
%     % break out contrast tuning curves
%     clear parameters;
%     parameters.input_parameters.independent_label='Contrast';
%     parameters.input_parameters.independent_parameter = 'contrast';
%     parameters.input_parameters.best_algorithm = 'empirical_maximum';
%     parameters.input_parameters.depends_on = struct('name','stimulus_response_scalar_id','value','');
%     parameters.input_parameters.selection = struct('property','angle','operation','exact_number','value','best');
%     parameters.input_parameters.selection(2) = struct('property','sFrequency','operation','exact_number','value','deal');
%     parameters.input_parameters.selection(3) = struct('property','sFrequency','operation','hasfield','value','varies');
%     parameters.input_parameters.selection(4) = struct('property','contrast','operation','hasfield','value','varies');
%     parameters.input_parameters.selection(5) = struct('property','angle','operation','hasfield','value','varies');
%     parameters.depends_on = vlt.data.emptystruct('name','value');
%     docs = c.run(vlt.data.conditional(CONT_TC_FORCE,'NoAction','Replace'),parameters);
% end;
% 
% if CONTRAST_FITS,
%     cc = ndi.calc.vis.contrast_tuning(S);
%     cc.run(vlt.data.conditional(CONTRAST_FITS_FORCE,'NoAction','Replace'));
% end;
% 
% if CSENSITIVITY_GROUPING,
%     ccc = ndi.calc.vis.contrast_sensitivity(S);
%     clear parameters;
%     parameters = ccc.default_search_for_input_parameters();
%     parameters.query.query = ndi.query('element.type','exact_string','spikes','');
%     ccc.run(vlt.data.conditional(CSENSITIVITY_GROUPING_FORCE,'NoAction','Replace'),parameters);
% end;
% 
% if ORI_TUNING_CURVES,
%     c = ndi.calc.stimulus.tuningcurve(S);
%     % break out contrast tuning curves
%     clear parameters;
%     parameters.input_parameters.independent_label='Orientation';
%     parameters.input_parameters.independent_parameter = 'angle';
%     parameters.input_parameters.best_algorithm = 'empirical_maximum';
%     parameters.input_parameters.depends_on = struct('name','stimulus_response_scalar_id','value','');
%     parameters.input_parameters.selection = struct('property','sFrequency','operation','exact_number','value','best');
%     parameters.input_parameters.selection(2) = struct('property','sFrequency','operation','hasfield','value','varies');
%     parameters.input_parameters.selection(3) = struct('property','contrast','operation','exact_number','value',1);
%     parameters.input_parameters.selection(4) = struct('property','angle','operation','hasfield','value','varies');
%     parameters.depends_on = vlt.data.emptystruct('name','value');
%     docs = c.run(vlt.data.conditional(OT_TC_FORCE,'NoAction','Replace'),parameters);
% end;
% 
% if ORIENTATION_FITS,
%     co = ndi.calc.vis.oridir_tuning(S);
%     parameters = co.default_search_for_input_parameters();
%     parameters.query.query = ndi.query('tuning_curve.independent_variable_label','exact_string_anycase','Orientation','');
%     co.run(vlt.data.conditional(ORIENTATION_FITS_FORCE,'NoAction','Replace'),parameters);
% end;

if DIR_TUNING_CURVES,
    c_d = ndi.calc.stimulus.tuningcurve(S);
    % break out contrast tuning curves
    clear parameters;
    parameters.input_parameters.independent_label='Direction';
    parameters.input_parameters.independent_parameter = 'angle';
    parameters.input_parameters.best_algorithm = 'empirical_maximum';
    parameters.input_parameters.depends_on = struct('name','stimulus_response_scalar_id','value','');
    parameters.input_parameters.selection = struct('property','angle','operation','hasfield','value','varies');
    parameters.depends_on = vlt.data.emptystruct('name','value');
    c_d.run(vlt.data.conditional(DIR_TC_FORCE,'NoAction','Replace'),parameters);
end;

if DIR_FITS,
    co = ndi.calc.vis.oridir_tuning(S);
    parameters = co.default_search_for_input_parameters();
    parameters.query.query = ndi.query('tuning_curve.independent_variable_label','exact_string_anycase','Direction','');
    co.run(vlt.data.conditional(0,'NoAction','Replace'),parameters);
end;


if BOOTSTRAPPING
    bs = ndi.calc.vis.bootstrapping(S);
    parameters = bs.default_search_for_input_parameters();
%     doc = bs.run(vlt.data.conditional(1,'NoAction','Replace'),parameters);
end;