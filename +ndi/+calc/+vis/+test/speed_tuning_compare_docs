function [b, errormsg] = speed_tuning_compare_docs(document_expected, document_actual, scope)
% Speed_Tuning_Compare_Docs
%
% [B, ERRORMSG] = ndi.calc.vis.test.speed_tuning_compare_docs(DOC_EXPECTED, DOC_ACTUAL, SCOPE)
%


b = 1;
errormsg = '';

% start comparison

doc_e = document_expected.document_properties.speed_tuning;
doc_a = document_actual.document_properties.speed_tuning;

% Comparing Properties
% Comparing Properties.response_units

properties_match = strcmpi(doc_e.properties.response_units, doc_a.properties.response_units);
if ~properties_match
   b = 0;
   errormsg = ['Expected response units in ' doc_e.properties.response_units ' but observed ' doc_a.properties.response_units];
   return;
end

% Comparing Properties.response_type

properties_match = strcmpi(doc_e.properties.response_type, doc_a.properties.response_type);
if ~properties_match
   b = 0;
   errormsg = ['Expected response type of ' doc_e.properties.response_units ' but observed ' doc_a.properties.response_units];
   return;
end

% Comparing tuning_curve
%   spatial_frequency				
%   temporal_frequency				
%   mean					
%   stddev						
%   stderr				
%   individual				
%   raw_individual					
%   control_individual			

[b,errormsg] = ndi.test.values_within_tolerance(doc_e.tuning_curve.spatial_frequency, doc_a.tuning_curve.spatial_frequency, tolerance, 'spatial frequency');
[b,errormsg] = ndi.test.values_within_tolerance(doc_e.tuning_curve.temporal_frequency, doc_a.tuning_curve.temporal_frequency, tolerance, 'temporal frequency');
[b,errormsg] = ndi.test.values_within_tolerance(doc_e.tuning_curve.mean, doc_a.tuning_curve.mean, tolerance, 'mean');
[b,errormsg] = ndi.test.values_within_tolerance(doc_e.tuning_curve.stddev, doc_a.tuning_curve.stddev, tolerance, 'stddev');
[b,errormsg] = ndi.test.values_within_tolerance(doc_e.tuning_curve.stderr, doc_a.tuning_curve.stderr, tolerance, 'stderr');
[b,errormsg] = ndi.test.values_within_tolerance(doc_e.tuning_curve.individual, doc_a.tuning_curve.individual, tolerance, 'individual');
[b,errormsg] = ndi.test.values_within_tolerance(doc_e.tuning_curve.raw_individual, doc_a.tuning_curve.raw_individual, tolerance, 'raw individual');
[b,errormsg] = ndi.test.values_within_tolerance(doc_e.tuning_curve.control_individual, doc_a.tuning_curve.control_individual, tolerance, 'control individual');

% Comparing significance
%   visual_response_anova_p
%   across_stimuli_anova_p

[b,errormsg] = ndi.test.values_within_tolerance(doc_e.significance.visual_response_anova_p, doc_a.significance.visual_response_anova_p, tolerance, 'visual response anova p');
[b,errormsg] = ndi.test.values_within_tolerance(doc_e.significance.across_stimuli_anova_p, doc_a.significance.across_stimuli_anova_p, tolerance, 'across stimuli anova p');

% Comparing fit:
%   Priebe_fit_parameters
%   Priebe_fit_spatial_frequencies
%   Priebe_fit_temporal_frequencies
%   Priebe_fit_values
%   Priebe_fit_speed_tuning_index
%   Priebe_fit_spatial_frequency_preference
%   Priebe_fit_temporal_frequency_preference

[b,errormsg] = ndi.test.values_within_tolerance(doc_e.fit.Priebe_fit_parameters, doc_a.fit.Priebe_fit_parameters, tolerance, 'Priebe fit parameters');
[b,errormsg] = ndi.test.values_within_tolerance(doc_e.fit.Priebe_fit_spatial_frequencies, doc_a.fit.Priebe_fit_spatial_frequencies, tolerance, 'Priebe fit spatial frequencies');
[b,errormsg] = ndi.test.values_within_tolerance(doc_e.fit.Priebe_fit_temporal_frequencies, doc_a.fit.Priebe_fit_temporal_frequencies, tolerance, 'Priebe fit temporal frequencies');
[b,errormsg] = ndi.test.values_within_tolerance(doc_e.fit.Priebe_fit_values, doc_a.fit.Priebe_fit_values, tolerance, 'Priebe fit values');
[b,errormsg] = ndi.test.values_within_tolerance(doc_e.fit.Priebe_fit_speed_tuning_index, doc_a.fit.Priebe_fit_speed_tuning_index, tolerance, 'Priebe fit speed tuning index');
[b,errormsg] = ndi.test.values_within_tolerance(doc_e.fit.Priebe_fit_spatial_frequency_preference, doc_a.fit.Priebe_fit_spatial_frequency_preference, tolerance, 'Priebe fit spatial frequency preference');
[b,errormsg] = ndi.test.values_within_tolerance(doc_e.fit.Priebe_fit_temporal_frequency_preference, doc_a.fit.Priebe_fit_temporal_frequency_preference, tolerance, 'Priebe fit temporal frequency preference');

end

