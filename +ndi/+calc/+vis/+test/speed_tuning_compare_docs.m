function [b, errormsg] = speed_tuning_compare_docs(document_expected, document_actual, scope)
% Comparing Speed_tuning
%
% [B, ERRORMSG] = ndi.calc.vis.test.speed_tuning_compare_docs(DOC_EXPECTED, DOC_ACTUAL, SCOPE)
%
%


b = 1;
errormsg = '';

% start comparison 

doc_e = document_expected.document_properties.speed_tuning;
doc_a = document_actual.document_properties.speed_tuning;


% Comparing tuning_curve variables
% Comparing spatial_frequency
tol = %put in appropriate tolorance;

spatial_frequency_match = vlt.data.sizeeq(doc_e.tuning_curve.spatial_frequency(:),doc_a.tuning_curve.spatial_frequency(:));
if spatial_frequency_match
   spatial_frequency_match = max(abs(doc_e.tuning_curve.spatial_frequency(:) - doc_a.tuning_curve.spatial_frequency(:))) < tol;
end

if ~spatial_frequency_match
   b = 0;
   errormsg = ['Spatial Frequency differed by greater than ' num2str(tol) '.'];
   return;
end

% Comparing temporal_frequency
tol = %put in appropriate tolorance;

temporal_frequency_match = vlt.data.sizeeq(doc_e.tuning_curve.temporal_frequency(:),doc_a.tuning_curve.temporal_frequency(:));
if temporal_frequency_match
   temporal_frequency_match = max(abs(doc_e.tuning_curve.temporal_frequency(:) - doc_a.tuning_curve.temporal_frequency(:))) < tol;
end

if ~temporal_frequency_match
   b = 0;
   errormsg = ['Temporal Frequency differed by greater than ' num2str(tol) '.'];
   return;
end

% Comparing mean
tol = %put in appropriate tolorance;

mean_match = vlt.data.sizeeq(doc_e.tuning_curve.mean(:),doc_a.tuning_curve.mean(:));
if mean_match
   mean_match = max(abs(doc_e.tuning_curve.mean(:) - doc_a.tuning_curve.mean(:))) < tol;
end

if ~mean_match
   b = 0;
   errormsg = ['Mean differed by greater than ' num2str(tol) '.'];
   return;
end

% Comparing stddev
tol = %put in appropriate tolorance;

stddev_match = vlt.data.sizeeq(doc_e.tuning_curve.stddev(:),doc_a.tuning_curve.stddev(:));
if stddev_match
   stddev_match = max(abs(doc_e.tuning_curve.stddev(:) - doc_a.tuning_curve.stddev(:))) < tol;
end

if ~stddev_match
   b = 0;
   errormsg = ['Stddev differed by greater than ' num2str(tol) '.'];
   return;
end

% Comparing stderr
tol = %put in appropriate tolorance;

stderr_match = vlt.data.sizeeq(doc_e.tuning_curve.stderr(:),doc_a.tuning_curve.stderr(:));
if stderr_match
   stderr_match = max(abs(doc_e.tuning_curve.stderr(:) - doc_a.tuning_curve.stderr(:))) < tol;
end

if ~stderr_match
   b = 0;
   errormsg = ['Stderr differed by greater than ' num2str(tol) '.'];
   return;
end


% Comparing significance variables
% Comparing visual_response_anova_p
tol = %put in appropriate tolorance;

visual_response_anova_p_match = vlt.data.sizeeq(doc_e.significance.visual_response_anova_p(:),doc_a.significance.visual_response_anova_p(:));
if visual_response_anova_p_match
   visual_response_anova_p_match = max(abs(doc_e.significance.visual_response_anova_p(:) - doc_a.significance.visual_response_anova_p(:))) < tol;
end

if ~visual_response_anova_p_match
   b = 0;
   errormsg = ['Visual Response Anova P by greater than ' num2str(tol) '.'];
   return;
end

% Comparing across_stimuli_anova_p
tol = %put in appropriate tolorance;

across_stimuli_anova_p_match = vlt.data.sizeeq(doc_e.significance.across_stimuli_anova_p(:),doc_a.significance.across_stimuli_anova_p(:));
if across_stimuli_anova_p_match
   across_stimuli_anova_p_match = max(abs(doc_e.significance.across_stimuli_anova_p(:) - doc_a.significance.across_stimuli_anova_p(:))) < tol;
end

if ~across_stimuli_anova_p_match
   b = 0;
   errormsg = ['Across Stimuli Anova P by greater than ' num2str(tol) '.'];
   return;
end


% Comparing fit variables
% Comparing Priebe_fit_parameters
tol = %put in appropriate tolorance;

Priebe_fit_parameters_match = vlt.data.sizeeq(doc_e.fit.Priebe_fit_parameters(:),doc_a.fit.Priebe_fit_parameters(:));
if Priebe_fit_parameters_match
   Priebe_fit_parameters_match = max(abs(doc_e.fit.Priebe_fit_parameters(:) - doc_a.fit.Priebe_fit_parameters(:))) < tol;
end

if ~Priebe_fit_parameters_match
   b = 0;
   errormsg = ['Priebe fit parameters by greater than ' num2str(tol) '.'];
   return;
end

% Comparing Priebe_fit_spatial_frequencies
tol = %put in appropriate tolorance;

Priebe_fit_spatial_frequencies_match = vlt.data.sizeeq(doc_e.fit.Priebe_fit_spatial_frequencies(:),doc_a.fit.Priebe_fit_spatial_frequencies(:));
if Priebe_fit_spatial_frequencies_match
   Priebe_fit_spatial_frequencies_match = max(abs(doc_e.fit.Priebe_fit_spatial_frequencies(:) - doc_a.fit.Priebe_fit_spatial_frequencies(:))) < tol;
end

if ~Priebe_fit_spatial_frequencies_match
   b = 0;
   errormsg = ['Priebe fit spatial frequency by greater than ' num2str(tol) '.'];
   return;
end

% Comparing Priebe_fit_temporal_frequencies
tol = %put in appropriate tolorance;

Priebe_fit_temporal_frequencies_match = vlt.data.sizeeq(doc_e.fit.Priebe_fit_temporal_frequencies(:),doc_a.fit.Priebe_fit_temporal_frequencies(:));
if Priebe_fit_temporal_frequencies_match
   Priebe_fit_temporal_frequencies_match = max(abs(doc_e.fit.Priebe_fit_temporal_frequencies(:) - doc_a.fit.Priebe_fit_temporal_frequencies(:))) < tol;
end

if ~Priebe_fit_temporal_frequencies_match
   b = 0;
   errormsg = ['Priebe fit temporal frequencies by greater than ' num2str(tol) '.'];
   return;
end

% Comparing Priebe_fit_values
tol = %put in appropriate tolorance;

Priebe_fit_values_match = vlt.data.sizeeq(doc_e.fit.Priebe_fit_values(:),doc_a.fit.Priebe_fit_values(:));
if Priebe_fit_values_match
   Priebe_fit_values_match = max(abs(doc_e.fit.Priebe_fit_values(:) - doc_a.fit.Priebe_fit_values(:))) < tol;
end

if ~Priebe_fit_values_match
   b = 0;
   errormsg = ['Priebe fit values by greater than ' num2str(tol) '.'];
   return;
end

% Comparing Priebe_fit_speed_tuning_index
tol = %put in appropriate tolorance;

Priebe_fit_speed_tuning_index_match = vlt.data.sizeeq(doc_e.fit.Priebe_fit_speed_tuning_index(:),doc_a.fit.Priebe_fit_speed_tuning_index(:));
if Priebe_fit_speed_tuning_index_match
   Priebe_fit_speed_tuning_index_match = max(abs(doc_e.fit.Priebe_fit_speed_tuning_index(:) - doc_a.fit.Priebe_fit_speed_tuning_index(:))) < tol;
end

if ~Priebe_fit_speed_tuning_index_match
   b = 0;
   errormsg = ['Priebe fit speed tuning index by greater than ' num2str(tol) '.'];
   return;
end

% Comparing Priebe_fit_spatial_frequency_preference
tol = %put in appropriate tolorance;

Priebe_fit_spatial_frequency_preference_match = vlt.data.sizeeq(doc_e.fit.Priebe_fit_spatial_frequency_preference(:),doc_a.fit.Priebe_fit_spatial_frequency_preference(:));
if Priebe_fit_spatial_frequency_preference_match
   Priebe_fit_spatial_frequency_preference_match = max(abs(doc_e.fit.Priebe_fit_spatial_frequency_preference(:) - doc_a.fit.Priebe_fit_spatial_frequency_preference(:))) < tol;
end

if ~Priebe_fit_spatial_frequency_preference_match
   b = 0;
   errormsg = ['Priebe fit spatial frequency preference by greater than ' num2str(tol) '.'];
   return;
end

% Comparing Priebe_fit_temporal_frequency_preference
tol = %put in appropriate tolorance;

Priebe_fit_temporal_frequency_preference_match = vlt.data.sizeeq(doc_e.fit.Priebe_fit_temporal_frequency_preference(:),doc_a.fit.Priebe_fit_temporal_frequency_preference(:));
if Priebe_fit_temporal_frequency_preference_match
   Priebe_fit_temporal_frequency_preference_match = max(abs(doc_e.fit.Priebe_fit_temporal_frequency_preference(:) - doc_a.fit.Priebe_fit_temporal_frequency_preference(:))) < tol;
end

if ~Priebe_fit_temporal_frequency_preference_match
   b = 0;
   errormsg = ['Priebe fit temporal frequency preference by greater than ' num2str(tol) '.'];
   return;
end


