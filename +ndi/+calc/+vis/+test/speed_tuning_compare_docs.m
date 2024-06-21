function [b_, errormsg_] = speed_tuning_compare_docs(document_expected, document_actual, scope)
% Speed_Tuning_Compare_Docs
%
% [B, ERRORMSG] = ndi.calc.vis.test.speed_tuning_compare_docs(DOC_EXPECTED, DOC_ACTUAL, SCOPE)
%

% Initialize b_ as a row vector of ones for tracking comparison results.
% Initialize errormsg_ as an empty cell array to hold error messages.

b_ = ones(1,19);
errormsg_ = cell(1,19);

 % establish scope-dependent tolerances
switch(scope),
    case 'standard',
    
       tolerance.tuning_curve.spatial_frequency = 1e-6;				
       tolerance.tuning_curve.temporal_frequency = 1e-6;				
       tolerance.tuning_curve.mean = 1e-6;					
       tolerance.tuning_curve.stddev = 1e-6;			
       tolerance.tuning_curve.stderr = 1e-6;			
       tolerance.tuning_curve.individual = 1e-6;					
       tolerance.tuning_curve.control_stddv = 1e-6;
       tolerance.tuning_curve.control_stderr = 1e-6;
       tolerance.significance.visual_response_anova_p = 1e-6;	
       tolerance.significance.across_stimuli_anova_p = 1e-6;	
       tolerance.fit.Priebe_fit_parameters = [1e-3 1e-3 1 1e-1 1e-2 1e-3 1e-2];
       tolerance.fit.Priebe_fit_spatial_frequencies = 1e-6;
       tolerance.fit.Priebe_fit_temporal_frequencies = 1e-6;
       tolerance.fit.Priebe_fit_values = 1e-6;
       tolerance.fit.Priebe_fit_speed_tuning_index = 1e-6;
       tolerance.fit.Priebe_fit_spatial_frequency_preference = 1e-6;
       tolerance.fit.Priebe_fit_temporal_frequency_preference = 1e-6;

    case 'low_noise',

       tolerance.tuning_curve.spatial_frequency = 1e-6;				
       tolerance.tuning_curve.temporal_frequency = 1e-6;				
       tolerance.tuning_curve.mean = 1e-6;					
       tolerance.tuning_curve.stddev = 1e-6;			
       tolerance.tuning_curve.stderr = 1e-6;			
       tolerance.tuning_curve.individual = 1e-6;					
       tolerance.tuning_curve.control_stddev = 1e-6;
       tolerance.tuning_curve.control_stderr = 1e-6;
       tolerance.significance.visual_response_anova_p = 1e-6;	
       tolerance.significance.across_stimuli_anova_p = 1e-6;	
       tolerance.fit.Priebe_fit_parameters = [1e-3 1e-3 1 1e-1 1e-2 1e-3 1e-2];
       tolerance.fit.Priebe_fit_spatial_frequencies = 1e-6;
       tolerance.fit.Priebe_fit_temporal_frequencies = 1e-6;
       tolerance.fit.Priebe_fit_values = 1e-6;
       tolerance.fit.Priebe_fit_speed_tuning_index = 1e-6;
       tolerance.fit.Priebe_fit_spatial_frequency_preference = 1e-6;
       tolerance.fit.Priebe_fit_temporal_frequency_preference = 1e-6;

    case 'high_noise',

       tolerance.tuning_curve.spatial_frequency = 1e-6;				
       tolerance.tuning_curve.temporal_frequency = 1e-6;				
       tolerance.tuning_curve.mean = 1e-6;					
       tolerance.tuning_curve.stddev = 1e-6;			
       tolerance.tuning_curve.stderr = 1e-6;			
       tolerance.tuning_curve.individual = 1e-6;					
       tolerance.tuning_curve.control_stddev = 1e-6;
       tolerance.tuning_curve.control_stderr = 1e-6;
       tolerance.significance.visual_response_anova_p = 1e-6;	
       tolerance.significance.across_stimuli_anova_p = 1e-6;	
       tolerance.fit.Priebe_fit_parameters = [1e-3 1e-3 1 1e-1 1e-2 1e-3 1e-2];
       tolerance.fit.Priebe_fit_spatial_frequencies = 1e-6;
       tolerance.fit.Priebe_fit_temporal_frequencies = 1e-6;
       tolerance.fit.Priebe_fit_values = 1e-6;
       tolerance.fit.Priebe_fit_speed_tuning_index = 1e-6;
       tolerance.fit.Priebe_fit_spatial_frequency_preference = 1e-6;
       tolerance.fit.Priebe_fit_temporal_frequency_preference = 1e-6;

    otherwise,
       error(['Unknown scope ' scope '.']);
end;

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

[b_(1),errormsg_{1}] = ndi.test.values_within_tolerance(doc_e.tuning_curve.spatial_frequency, doc_a.tuning_curve.spatial_frequency, tolerance.tuning_curve.spatial_frequency, 'spatial frequency');
[b_(2),errormsg_{2}] = ndi.test.values_within_tolerance(doc_e.tuning_curve.temporal_frequency, doc_a.tuning_curve.temporal_frequency, tolerance.tuning_curve.temporal_frequency, 'temporal frequency');
[b_(3),errormsg_{3}] = ndi.test.values_within_tolerance(doc_e.tuning_curve.mean, doc_a.tuning_curve.mean, tolerance.tuning_curve.mean, 'mean');
[b_(4),errormsg_{4}] = ndi.test.values_within_tolerance(doc_e.tuning_curve.stddev, doc_a.tuning_curve.stddev, tolerance.tuning_curve.stddev, 'stddev');
[b_(5),errormsg_{5}] = ndi.test.values_within_tolerance(doc_e.tuning_curve.stderr, doc_a.tuning_curve.stderr, tolerance.tuning_curve.stderr, 'stderr');
[b_(6),errormsg_{6}] = ndi.test.values_within_tolerance(doc_e.tuning_curve.individual, doc_a.tuning_curve.individual, tolerance.tuning_curve.individual, 'individual');
%[b_(7),errormsg_{7}] = ndi.test.values_within_tolerance(doc_e.tuning_curve.raw_individual, doc_a.tuning_curve.raw_individual, tolerance.tuning_curve.raw_individual, 'raw individual');
%[b_(8),errormsg_{8}] = ndi.test.values_within_tolerance(doc_e.tuning_curve.control_individual, doc_a.tuning_curve.control_individual, tolerance.tuning_curve.control_individual, 'control individual');
[b_(7),errormsg_{7}] = ndi.test.values_within_tolerance(doc_e.tuning_curve.control_stddev, doc_a.tuning_curve.control_stddev, tolerance.tuning_curve.control_stddev, 'control stddev');
[b_(8),errormsg_{8}] = ndi.test.values_within_tolerance(doc_e.tuning_curve.control_stderr, doc_a.tuning_curve.control_stderr, tolerance.tuning_curve.control_stderr, 'control stderr');


% Comparing significance
%   visual_response_anova_p
%   across_stimuli_anova_p

[b_(9),errormsg_{9}] = ndi.test.values_within_tolerance(doc_e.significance.visual_response_anova_p, doc_a.significance.visual_response_anova_p, tolerance.significance.visual_response_anova_p, 'visual response anova p');
[b_(10),errormsg_{10}] = ndi.test.values_within_tolerance(doc_e.significance.across_stimuli_anova_p, doc_a.significance.across_stimuli_anova_p, tolerance.significance.across_stimuli_anova_p, 'across stimuli anova p');

% Comparing fit:
%   Priebe_fit_parameters
%   Priebe_fit_spatial_frequencies
%   Priebe_fit_temporal_frequencies
%   Priebe_fit_values
%   Priebe_fit_speed_tuning_index
%   Priebe_fit_spatial_frequency_preference
%   Priebe_fit_temporal_frequency_preference

[b_(11),errormsg_{11}] = ndi.test.values_within_tolerance(doc_e.fit.Priebe_fit_parameters(1), doc_a.fit.Priebe_fit_parameters(1), tolerance.fit.Priebe_fit_parameters(1), 'Priebe fit parameters');
[b_(12),errormsg_{12}] = ndi.test.values_within_tolerance(doc_e.fit.Priebe_fit_parameters(2), doc_a.fit.Priebe_fit_parameters(2), tolerance.fit.Priebe_fit_parameters(2), 'Priebe fit parameters');
[b_(13),errormsg_{13}] = ndi.test.values_within_tolerance(doc_e.fit.Priebe_fit_parameters(3), doc_a.fit.Priebe_fit_parameters(3), tolerance.fit.Priebe_fit_parameters(3), 'Priebe fit parameters');
[b_(14),errormsg_{14}] = ndi.test.values_within_tolerance(doc_e.fit.Priebe_fit_parameters(4), doc_a.fit.Priebe_fit_parameters(4), tolerance.fit.Priebe_fit_parameters(4), 'Priebe fit parameters');
[b_(15),errormsg_{15}] = ndi.test.values_within_tolerance(doc_e.fit.Priebe_fit_parameters(5), doc_a.fit.Priebe_fit_parameters(5), tolerance.fit.Priebe_fit_parameters(5), 'Priebe fit parameters');
[b_(16),errormsg_{16}] = ndi.test.values_within_tolerance(doc_e.fit.Priebe_fit_parameters(6), doc_a.fit.Priebe_fit_parameters(6), tolerance.fit.Priebe_fit_parameters(6), 'Priebe fit parameters');
[b_(17),errormsg_{17}] = ndi.test.values_within_tolerance(doc_e.fit.Priebe_fit_spatial_frequencies, doc_a.fit.Priebe_fit_spatial_frequencies, tolerance.fit.Priebe_fit_spatial_frequencies, 'Priebe fit spatial frequencies');
[b_(18),errormsg_{18}] = ndi.test.values_within_tolerance(doc_e.fit.Priebe_fit_temporal_frequencies, doc_a.fit.Priebe_fit_temporal_frequencies, tolerance.fit.Priebe_fit_temporal_frequencies, 'Priebe fit temporal frequencies');
[b_(19),errormsg_{19}] = ndi.test.values_within_tolerance(doc_e.fit.Priebe_fit_values, doc_a.fit.Priebe_fit_values, tolerance.fit.Priebe_fit_values, 'Priebe fit values');
[b_(20),errormsg_{20}] = ndi.test.values_within_tolerance(doc_e.fit.Priebe_fit_speed_tuning_index, doc_a.fit.Priebe_fit_speed_tuning_index, tolerance.fit.Priebe_fit_speed_tuning_index, 'Priebe fit speed tuning index');
[b_(21),errormsg_{21}] = ndi.test.values_within_tolerance(doc_e.fit.Priebe_fit_spatial_frequency_preference, doc_a.fit.Priebe_fit_spatial_frequency_preference, tolerance.fit.Priebe_fit_spatial_frequency_preference, 'Priebe fit spatial frequency preference');
[b_(22),errormsg_{22}] = ndi.test.values_within_tolerance(doc_e.fit.Priebe_fit_temporal_frequency_preference, doc_a.fit.Priebe_fit_temporal_frequency_preference, tolerance.fit.Priebe_fit_temporal_frequency_preference, 'Priebe fit temporal frequency preference');

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

