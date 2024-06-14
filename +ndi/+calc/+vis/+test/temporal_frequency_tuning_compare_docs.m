function [b_, errormsg_] = temporal_frequency_tuning_compare_docs(document_expected, document_actual, scope)
% Comparing temporal_frequency_tuning
%
% [B, ERRORMSG] = ndi.calc.vis.test.temporal_frequency_tuning_compare_docs(DOC_EXPECTED, DOC_ACTUAL, SCOPE)
%
%

% Initialize b_ as a row vector of ones for tracking comparison results
% Initialize errormsg_ as an empty cell array to hold error messages

b_ = ones(1,39);
errormsg_ = cell(1,39);

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
   b_(38) = 0;
   errormsg_{38} = ['Expected response units in ' doc_e.properties.response_units ' but observed ' doc_a.properties.response_units];
   return;
end

%   Response Type

properties_match = strcmpi(char(doc_e.properties.response_type), char(doc_a.properties.response_type));
if ~properties_match
   b_(39) = 0;
   errormsg_{39} = ['Expected response type of ' doc_e.properties.response_type ' but observed ' doc_a.properties.response_type];
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

[b_(1),errormsg_{1}] = ndi.test.values_within_tolerance(doc_e.tuning_curve.temporal_frequency, doc_a.tuning_curve.temporal_frequency, tol_tuning_curve, 'temporal_frequency');
[b_(2),errormsg_{2}] = ndi.test.values_within_tolerance(doc_e.tuning_curve.mean, doc_a.tuning_curve.mean, tol_tuning_curve, 'mean');
[b_(3),errormsg_{3}] = ndi.test.values_within_tolerance(doc_e.tuning_curve.stddev, doc_a.tuning_curve.stddev, tol_tuning_curve, 'stddev');
[b_(4),errormsg_{4}] = ndi.test.values_within_tolerance(doc_e.tuning_curve.stderr, doc_a.tuning_curve.stderr, tol_tuning_curve, 'stderr');
[b_(5),errormsg_{5}] = ndi.test.values_within_tolerance(doc_e.tuning_curve.individual, doc_a.tuning_curve.individual, tol_tuning_curve, 'individual');
[b_(6),errormsg_{6}] = ndi.test.values_within_tolerance(doc_e.tuning_curve.control_stddev, doc_a.tuning_curve.control_stddev, tol_tuning_curve, 'control_stddev');
[b_(7),errormsg_{7}] = ndi.test.values_within_tolerance(doc_e.tuning_curve.control_stderr, doc_a.tuning_curve.control_stderr, tol_tuning_curve, 'control_stderr');

% Comparing Significance
%   visual_response_anova_p
%   across_stimuli_anova_p

[b_(8),errormsg_{8}] = ndi.test.values_within_tolerance(doc_e.significance.visual_response_anova_p, doc_a.significance.visual_response_anova_p, tol_significance, 'visual response anova p');
[b_(9),errormsg_{9}] = ndi.test.values_within_tolerance(doc_e.significance.across_stimuli_anova_p, doc_a.significance.across_stimuli_anova_p, tol_significance, 'across stimuli anova p');

% Comparing Fitless
%   H50
%   Pref
%   L50

[b_(10),errormsg_{10}] = ndi.test.values_within_tolerance(doc_e.fitless.H50, doc_a.fitless.H50, tol_fitless, 'H50');
[b_(11),errormsg_{11}] = ndi.test.values_within_tolerance(doc_e.fitless.Pref, doc_a.fitless.Pref, tol_fitless, 'Pref');
[b_(12),errormsg_{12}] = ndi.test.values_within_tolerance(doc_e.fitless.L50, doc_a.fitless.L50, tol_fitless, 'L50');

% Comparing fit_spline
%   fit
%   H50
%   Pref
%   L50
%   values

[b_(13),errormsg_{13}] = ndi.test.values_within_tolerance(doc_e.fit_spline.H50, doc_a.fit_spline.H50, tol_fit_spline, 'fit');
[b_(14),errormsg_{14}] = ndi.test.values_within_tolerance(doc_e.fit_spline.H50, doc_a.fit_spline.H50, tol_fit_spline, 'H50');
[b_(15),errormsg_{15}] = ndi.test.values_within_tolerance(doc_e.fit_spline.Pref, doc_a.fit_spline.Pref, tol_fit_spline, 'Pref');
[b_(16),errormsg_{16}] = ndi.test.values_within_tolerance(doc_e.fit_spline.L50, doc_a.fit_spline.L50, tol_fit_spline, 'L50');
[b_(17),errormsg_{17}] = ndi.test.values_within_tolerance(doc_e.fit_spline.L50, doc_a.fit_spline.L50, tol_fit_spline, 'values');

% Comparing fit_dog
%   fit
%   H50
%   Pref
%   L50
%   values
%   parameters

[b_(18),errormsg_{18}] = ndi.test.values_within_tolerance(doc_e.fit_dog.fit, doc_a.fit_dog.fit, tol_fit_dog, 'fit');
[b_(19),errormsg_{19}] = ndi.test.values_within_tolerance(doc_e.fit_dog.fit, doc_a.fit_dog.fit, tol_fit_dog, 'H50');
[b_(20),errormsg_{20}] = ndi.test.values_within_tolerance(doc_e.fit_dog.fit, doc_a.fit_dog.fit, tol_fit_dog, 'Pref');
[b_(21),errormsg_{21}] = ndi.test.values_within_tolerance(doc_e.fit_dog.fit, doc_a.fit_dog.fit, tol_fit_dog, 'L50');
[b_(22),errormsg_{22}] = ndi.test.values_within_tolerance(doc_e.fit_dog.fit, doc_a.fit_dog.fit, tol_fit_dog, 'values');
[b_(23),errormsg_{23}] = ndi.test.values_within_tolerance(doc_e.fit_dog.fit, doc_a.fit_dog.fit, tol_fit_dog, 'parameters');
                                                                    
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
