function [b_, errormsg_] = speed_tuning_compare_docs(document_expected, document_actual, scope)
% Speed_Tuning_Compare_Docs
%
% [B, ERRORMSG] = ndi.calc.vis.test.speed_tuning_compare_docs(DOC_EXPECTED, DOC_ACTUAL, SCOPE)
%

% Initialize b_ as a row vector of ones for tracking comparison results.
% Initialize errormsg_ as an empty cell array to hold error messages.

b_ = ones(1,19);
errormsg_ = cell(1,19);

% start comparison

doc_e = document_expected.document_properties.speed_tuning;
doc_a = document_actual.document_properties.speed_tuning;

% Comparing Properties
% Comparing Properties.response_units

properties_match = strcmpi(char(doc_e.properties.response_units), char(doc_a.properties.response_units));
if ~properties_match
   b_(18) = 0;
   errormsg_{18} = ['Expected response units in ' doc_e.properties.response_units ' but observed ' doc_a.properties.response_units];
   return;
end

% Comparing Properties.response_type

properties_match = strcmpi(char(doc_e.properties.response_type), char(doc_a.properties.response_type));
if ~properties_match
   b_(19) = 0;
   errormsg_{19} = ['Expected response type of ' doc_e.properties.response_type ' but observed ' doc_a.properties.response_type];
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

tolerance = 1e-6; %placeholder

[b_(1),errormsg_{1}] = ndi.test.values_within_tolerance(doc_e.tuning_curve.spatial_frequency, doc_a.tuning_curve.spatial_frequency, 1e-16, 'spatial frequency');
[b_(2),errormsg_{2}] = ndi.test.values_within_tolerance(doc_e.tuning_curve.temporal_frequency, doc_a.tuning_curve.temporal_frequency, 1e-16, 'temporal frequency');
[b_(3),errormsg_{3}] = ndi.test.values_within_tolerance(doc_e.tuning_curve.mean, doc_a.tuning_curve.mean, tolerance, 'mean');
[b_(4),errormsg_{4}] = ndi.test.values_within_tolerance(doc_e.tuning_curve.stddev, doc_a.tuning_curve.stddev, tolerance, 'stddev');
[b_(5),errormsg_{5}] = ndi.test.values_within_tolerance(doc_e.tuning_curve.stderr, doc_a.tuning_curve.stderr, tolerance, 'stderr');
[b_(6),errormsg_{6}] = ndi.test.values_within_tolerance(doc_e.tuning_curve.individual, doc_a.tuning_curve.individual, tolerance, 'individual');
%[b_(7),errormsg_{7}] = ndi.test.values_within_tolerance(doc_e.tuning_curve.raw_individual, doc_a.tuning_curve.raw_individual, tolerance, 'raw individual');
%[b_(8),errormsg_{8}] = ndi.test.values_within_tolerance(doc_e.tuning_curve.control_individual, doc_a.tuning_curve.control_individual, tolerance, 'control individual');
[b_(7),errormsg_{7}] = ndi.test.values_within_tolerance(doc_e.tuning_curve.control_stddev, doc_a.tuning_curve.control_stddev, tolerance, 'individual');
[b_(8),errormsg_{8}] = ndi.test.values_within_tolerance(doc_e.tuning_curve.control_stderr, doc_a.tuning_curve.control_stderr, tolerance, 'individual');


% Comparing significance
%   visual_response_anova_p
%   across_stimuli_anova_p

[b_(9),errormsg_{9}] = ndi.test.values_within_tolerance(doc_e.significance.visual_response_anova_p, doc_a.significance.visual_response_anova_p, tolerance, 'visual response anova p');
[b_(10),errormsg_{10}] = ndi.test.values_within_tolerance(doc_e.significance.across_stimuli_anova_p, doc_a.significance.across_stimuli_anova_p, tolerance, 'across stimuli anova p');

% Comparing fit:
%   Priebe_fit_parameters
%   Priebe_fit_spatial_frequencies
%   Priebe_fit_temporal_frequencies
%   Priebe_fit_values
%   Priebe_fit_speed_tuning_index
%   Priebe_fit_spatial_frequency_preference
%   Priebe_fit_temporal_frequency_preference

[b_(11),errormsg_{11}] = ndi.test.values_within_tolerance(doc_e.fit.Priebe_fit_parameters, doc_a.fit.Priebe_fit_parameters, tolerance, 'Priebe fit parameters');
[b_(12),errormsg_{12}] = ndi.test.values_within_tolerance(doc_e.fit.Priebe_fit_spatial_frequencies, doc_a.fit.Priebe_fit_spatial_frequencies, tolerance, 'Priebe fit spatial frequencies');
[b_(13),errormsg_{13}] = ndi.test.values_within_tolerance(doc_e.fit.Priebe_fit_temporal_frequencies, doc_a.fit.Priebe_fit_temporal_frequencies, tolerance, 'Priebe fit temporal frequencies');
[b_(14),errormsg_{14}] = ndi.test.values_within_tolerance(doc_e.fit.Priebe_fit_values, doc_a.fit.Priebe_fit_values, tolerance, 'Priebe fit values');
[b_(15),errormsg_{15}] = ndi.test.values_within_tolerance(doc_e.fit.Priebe_fit_speed_tuning_index, doc_a.fit.Priebe_fit_speed_tuning_index, tolerance, 'Priebe fit speed tuning index');
[b_(16),errormsg_{16}] = ndi.test.values_within_tolerance(doc_e.fit.Priebe_fit_spatial_frequency_preference, doc_a.fit.Priebe_fit_spatial_frequency_preference, tolerance, 'Priebe fit spatial frequency preference');
[b_(17),errormsg_{17}] = ndi.test.values_within_tolerance(doc_e.fit.Priebe_fit_temporal_frequency_preference, doc_a.fit.Priebe_fit_temporal_frequency_preference, tolerance, 'Priebe fit temporal frequency preference');


% The following code does three things:
% 1. Identify the b_ values with unmatched results
% 2. Update b_ to only include those
% 3. Update the corresponding errormsg_ messages

if any(b_==0),
    error_indices = find(b_==0);
    b_ = b_(error_indices);
    errormsg_ = errormsg_(error_indices);
end
if all(b_) && all(cellfun('isempty', errormsg_))
    disp('No error');
end

end

