function out = extract_tuning_curve(session, tuning_doc)
% EXTRACT_TUNING_CURVE - pull the speed-tuning curve from a tuning curve document
%
%  OUT = vis.speed.extract_tuning_curve(SESSION, TUNING_DOC)
%
%  Extracts the pieces that both ndi.calc.vis.speed_tuning and
%  ndi.calc.vis.speed_tuning_bootstrap read from a stimulus_tuningcurve
%  document, so the two calculators cannot drift in how they interpret their
%  shared input.
%
%  OUT is a struct with fields:
%    .properties         - struct with fields response_units and response_type
%    .resp               - the vhlab response structure from
%                          ndi.app.stimulus.tuning_response.tuningcurvedoc2vhlabrespstruct.
%                          resp.ind{k} holds the individual single-trial
%                          responses for condition k -- the resampling unit for
%                          the bootstrap -- aligned with resp.curve columns and
%                          with the spatial/temporal frequency rows below.
%    .spatial_frequency  - column vector, one entry per stimulus condition
%    .temporal_frequency - column vector, one entry per stimulus condition
%    .tuning_curve       - the tuning_curve struct stored by speed_tuning
%    .significance       - struct with visual_response_anova_p and
%                          across_stimuli_anova_p
%
%  See also: ndi.calc.vis.speed_tuning, ndi.calc.vis.speed_tuning_bootstrap,
%            vis.speed.fit_priebe_triplet

    arguments
        session
        tuning_doc
    end

    properties = struct();
    properties.response_units = tuning_doc.document_properties.stimulus_tuningcurve.response_units;

    stim_response_doc = session.database_search(ndi.query('base.id', ...
        'exact_string', tuning_doc.dependency_value('stimulus_response_scalar_id'), ''));
    if numel(stim_response_doc) ~= 1
        error('Could not find stimulus response scalar document.');
    end
    if iscell(stim_response_doc)
        stim_response_doc = stim_response_doc{1};
    end
    properties.response_type = stim_response_doc.document_properties.stimulus_response_scalar.response_type;

    resp = ndi.app.stimulus.tuning_response.tuningcurvedoc2vhlabrespstruct(tuning_doc);

    [anova_across_stims, anova_across_stims_blank] = neural_response_significance(resp);

    spatial_frequency  = vlt.data.colvec(tuning_doc.document_properties.stimulus_tuningcurve.independent_variable_value(:, 1));
    temporal_frequency = vlt.data.colvec(tuning_doc.document_properties.stimulus_tuningcurve.independent_variable_value(:, 2));

    tuning_curve = struct( ...
        'spatial_frequency', spatial_frequency, ...
        'temporal_frequency', temporal_frequency, ...
        'mean', vlt.data.colvec(resp.curve(2, :)), ...
        'stddev', vlt.data.colvec(resp.curve(3, :)), ...
        'stderr', vlt.data.colvec(resp.curve(4, :)), ...
        'individual', vlt.data.cellarray2mat(resp.ind), ...
        'control_stddev', resp.blankresp(2), ...
        'control_stderr', resp.blankresp(3));

    significance = struct('visual_response_anova_p', anova_across_stims_blank, ...
        'across_stimuli_anova_p', anova_across_stims);

    out = struct();
    out.properties         = properties;
    out.resp               = resp;
    out.spatial_frequency  = spatial_frequency;
    out.temporal_frequency = temporal_frequency;
    out.tuning_curve       = tuning_curve;
    out.significance       = significance;

end % extract_tuning_curve()
