function [b_, errormsg_] = temporal_frequency_tuning_compare_docs(document_expected, document_actual, scope)
% Comparing temporal_frequency_tuning
%
% [B, ERRORMSG] = ndi.calc.vis.test.temporal_frequency_tuning_compare_docs(DOC_EXPECTED, DOC_ACTUAL, SCOPE)
%
%

% Initialize b_ as a row vector of ones for tracking comparison results
% Initialize errormsg_ as an empty cell array to hold error messages

b_ = [];
errormsg_ = cell(0);

% establish scope-dependent tolerances
switch(scope),
    case 'standard',
    
       tol_tuning_curve = 1e-6;
       tol_significance = 1e-6;
       tol_fitless = 0.1;
       tol_fit_spline = 0.1;
       tol_fit_dog = 0.1;

    case 'low_noise',
    
       tol_tuning_curve = 1e-6;
       tol_significance = 1e-6;
       tol_fitless = 0.1;
       tol_fit_spline = 0.1;
       tol_fit_dog = 0.1;

    case 'high_noise',
    
       tol_tuning_curve = 1e-6;
       tol_significance = 1e-6;
       tol_fitless = 0.1;
       tol_fit_spline = 0.1;
       tol_fit_dog = 0.1;

    otherwise,
       error(['Unknown scope ' scope '.']);
end;


% start comparison

doc_e = document_expected.document_properties.temporal_frequency_tuning;
doc_a = document_actual.document_properties.temporal_frequency_tuning;

% Comparing properties
%   Response Units

properties_match = strcmpi(char(doc_e.properties.response_units), char(doc_a.properties.response_units));
if ~properties_match
   b_(end+1) = 0;
   errormsg_{end+1} = ['Expected response units in ' doc_e.properties.response_units ' but observed ' doc_a.properties.response_units];
   return;
end

%   Response Type

properties_match = strcmpi(char(doc_e.properties.response_type), char(doc_a.properties.response_type));
if ~properties_match
   b_(end+1) = 0;
   errormsg_{end+1} = ['Expected response type of ' doc_e.properties.response_type ' but observed ' doc_a.properties.response_type];
   return;
end

% Comparing Tuning_curve
%	temporal frequency                               
%	mean                                   
%	stddev                                 
%	stderr                                
%	individual                             
%	control_stddev                        
%	control_stderr
%   control_mean
%   control_mean_stddev
%   control_mean_stderr

[b_(end+1),errormsg_{end+1}] = ndi.test.values_within_tolerance(doc_e.tuning_curve.temporal_frequency, doc_a.tuning_curve.temporal_frequency, tol_tuning_curve, 'temporal_frequency');
[b_(end+1),errormsg_{end+1}] = ndi.test.values_within_tolerance(doc_e.tuning_curve.mean, doc_a.tuning_curve.mean, tol_tuning_curve, 'mean');
[b_(end+1),errormsg_{end+1}] = ndi.test.values_within_tolerance(doc_e.tuning_curve.stddev, doc_a.tuning_curve.stddev, tol_tuning_curve, 'stddev');
[b_(end+1),errormsg_{end+1}] = ndi.test.values_within_tolerance(doc_e.tuning_curve.stderr, doc_a.tuning_curve.stderr, tol_tuning_curve, 'stderr');
%[b_(end+1),errormsg_{end+1}] = ndi.test.values_within_tolerance(doc_e.tuning_curve.individual, doc_a.tuning_curve.individual, tol_tuning_curve, 'individual');
for i = 1:size(doc_a.tuning_curve.individual,1) %number of rows in individual corresponds to the number of specified reps
    [b_(end+1),errormsg_{end+1}] = ndi.test.values_within_tolerance(doc_e.tuning_curve.individual(1,:), doc_a.tuning_curve.individual(i,:), tol_tuning_curve, ['individual data row ',num2str(i)]);
end
[b_(end+1),errormsg_{end+1}] = ndi.test.values_within_tolerance(doc_e.tuning_curve.control_stddev, doc_a.tuning_curve.control_stddev, tol_tuning_curve, 'control_stddev');
[b_(end+1),errormsg_{end+1}] = ndi.test.values_within_tolerance(doc_e.tuning_curve.control_stderr, doc_a.tuning_curve.control_stderr, tol_tuning_curve, 'control_stderr');
[b_(end+1),errormsg_{end+1}] = ndi.test.values_within_tolerance(doc_e.tuning_curve.control_mean, doc_a.tuning_curve.control_mean, tol_tuning_curve, 'control_mean');
[b_(end+1),errormsg_{end+1}] = ndi.test.values_within_tolerance(doc_e.tuning_curve.control_mean_stddev, doc_a.tuning_curve.control_mean_stddev, tol_tuning_curve, 'control_mean_stddev');
[b_(end+1),errormsg_{end+1}] = ndi.test.values_within_tolerance(doc_e.tuning_curve.control_mean_stderr, doc_a.tuning_curve.control_mean_stderr, tol_tuning_curve, 'control_mean_stderr');



% Comparing Significance
%   visual_response_anova_p
%   across_stimuli_anova_p

[b_(end+1),errormsg_{end+1}] = ndi.test.values_within_tolerance(doc_e.significance.visual_response_anova_p, doc_a.significance.visual_response_anova_p, tol_significance, 'visual response anova p');
[b_(end+1),errormsg_{end+1}] = ndi.test.values_within_tolerance(doc_e.significance.across_stimuli_anova_p, doc_a.significance.across_stimuli_anova_p, tol_significance, 'across stimuli anova p');

% Comparing fitless
%   H50
%   Pref
%   L50
%   bandwidth
%   low_pass_index
%   high_pass_index

[b_(end+1),errormsg_{end+1}] = ndi.test.values_within_tolerance(doc_e.fitless.H50, doc_a.fitless.H50, tol_fitless, 'H50');
[b_(end+1),errormsg_{end+1}] = ndi.test.values_within_tolerance(doc_e.fitless.Pref, doc_a.fitless.Pref, tol_fitless, 'Pref');
[b_(end+1),errormsg_{end+1}] = ndi.test.values_within_tolerance(doc_e.fitless.L50, doc_a.fitless.L50, tol_fitless, 'L50');
[b_(end+1),errormsg_{end+1}] = ndi.test.values_within_tolerance(doc_e.fitless.bandwidth, doc_a.fitless.bandwidth, tol_fitless, 'bandwidth');
[b_(end+1),errormsg_{end+1}] = ndi.test.values_within_tolerance(doc_e.fitless.low_pass_index, doc_a.fitless.low_pass_index, tol_fitless, 'low_pass_index');
[b_(end+1),errormsg_{end+1}] = ndi.test.values_within_tolerance(doc_e.fitless.high_pass_index, doc_a.fitless.high_pass_index, tol_fitless, 'high_pass_index');

% Comparing fit_spline
%   fit
%   H50
%   Pref
%   L50
%   values
%   bandwidth

[b_(end+1),errormsg_{end+1}] = ndi.test.values_within_tolerance(doc_e.fit_spline.fit, doc_a.fit_spline.fit, tol_fit_spline, 'fit');
[b_(end+1),errormsg_{end+1}] = ndi.test.values_within_tolerance(doc_e.fit_spline.H50, doc_a.fit_spline.H50, tol_fit_spline, 'H50');
[b_(end+1),errormsg_{end+1}] = ndi.test.values_within_tolerance(doc_e.fit_spline.Pref, doc_a.fit_spline.Pref, tol_fit_spline, 'Pref');
[b_(end+1),errormsg_{end+1}] = ndi.test.values_within_tolerance(doc_e.fit_spline.L50, doc_a.fit_spline.L50, tol_fit_spline, 'L50');
[b_(end+1),errormsg_{end+1}] = ndi.test.values_within_tolerance(doc_e.fit_spline.values, doc_a.fit_spline.values, tol_fit_spline, 'values');
[b_(end+1),errormsg_{end+1}] = ndi.test.values_within_tolerance(doc_e.fitless.bandwidth, doc_a.fitless.bandwidth, tol_fitless, 'bandwidth');

% Comparing fit_dog
%   fit
%   H50
%   Pref
%   L50
%   values
%   parameters
%   bandwidth

[b_(end+1),errormsg_{end+1}] = ndi.test.values_within_tolerance(doc_e.fit_dog.fit, doc_a.fit_dog.fit, tol_fit_dog, 'fit');
[b_(end+1),errormsg_{end+1}] = ndi.test.values_within_tolerance(doc_e.fit_dog.H50, doc_a.fit_dog.H50, tol_fit_dog, 'H50');
[b_(end+1),errormsg_{end+1}] = ndi.test.values_within_tolerance(doc_e.fit_dog.Pref, doc_a.fit_dog.Pref, tol_fit_dog, 'Pref');
[b_(end+1),errormsg_{end+1}] = ndi.test.values_within_tolerance(doc_e.fit_dog.L50, doc_a.fit_dog.L50, tol_fit_dog, 'L50');
[b_(end+1),errormsg_{end+1}] = ndi.test.values_within_tolerance(doc_e.fit_dog.values, doc_a.fit_dog.values, tol_fit_dog, 'values');
[b_(end+1),errormsg_{end+1}] = ndi.test.values_within_tolerance(doc_e.fit_dog.parameters, doc_a.fit_dog.parameters, tol_fit_dog, 'parameters');
[b_(end+1),errormsg_{end+1}] = ndi.test.values_within_tolerance(doc_e.fit_dog.bandwidth, doc_a.fit_dog.bandwidth, tol_fit_dog, 'bandwidth');

% Comparing fit_gausslog
%   fit
%   H50
%   Pref
%   L50
%   values
%   parameters
%   bandwidth

% Comparing fit_dog
%   fit
%   H50
%   Pref
%   L50
%   values
%   parameters
%   bandwidth

% Comparing fit_movshon
%   fit
%   H50
%   Pref
%   L50
%   values
%   parameters
%   R2
%   bandwidth

% Comparing fit_movshon_c
%   fit
%   H50
%   Pref
%   L50
%   values
%   parameters
%   R2
%   bandwidth

% Comparing abs

% Identify the b_ values with unmatched results
% Update b_ to only include those
% Update the corresponding errormsg_ messages

if any(b_==0),
    error_indices = find(b_==0);
    b_ = b_(error_indices);
    errormsg_ = errormsg_(error_indices);
end

if all(b_) && all(cellfun('isempty', errormsg_))
    disp('No error');
end


end
