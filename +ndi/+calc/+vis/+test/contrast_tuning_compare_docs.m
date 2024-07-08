function [b_, errormsg_] = contrast_tuning_compare_docs(document_expected, document_actual, scope)
% Comparing Contrast_tuning
%
% [B, ERRORMSG] = ndi.calc.vis.test.contrast_tuning_compare_docs(DOC_EXPECTED, DOC_ACTUAL, SCOPE)
%
%

% Initialize b_ as a row vector of ones for tracking comparison results
% Initialize errormsg_ as an empty cell array to hold error messages

b_ = ones(1,39);
errormsg_ = cell(1,39);


% establish scope-dependent tolerances
switch(scope),
    case 'standard',
    
       tolerance.tuning_curve.contrast = 1e-6;
       tolerance.tuning_curve.mean = 1e-6;
       tolerance.tuning_curve.stddev = 1e-6;
       tolerance.tuning_curve.stderr = 1e-6;
       tolerance.tuning_curve.individual = 0.5;
       tolerance.tuning_curve.control_stddev = 1e-6;
       tolerance.tuning_curve.control_stderr = 1e-6;
       tolerance.significance.visual_response_anova_p = 0.1;
       tolerance.significance.across_stimuli_anova_p = 0.1;
       tolerance.fitless.interpolated_c50 = 0.1;
       tolerance.fit.naka_rushton_RB_parameters(1) = 5; %max response (0 to ~50)
       tolerance.fit.naka_rushton_RB_parameters(2) = 0.1; %value of contrasts (0 to 1) such that R(C50) = 0.5 * max(R)
       tolerance.fit.naka_rushton_RB_contrast = 0.1; %0 to 1
       tolerance.fit.naka_rushton_RB_values = 5; %responses (0 to ~50)
       tolerance.fit.naka_rushton_RB_pref = 0.1;	%0 to 1
       tolerance.fit.naka_rushton_RB_empirical_c50 = 0.1; %0 to 1
       tolerance.fit.naka_rushton_RB_r2 = 0.1; %0 to 1
       tolerance.fit.naka_rushton_RB_relative_max_gain = 5; %max(dr/dc) where dc is <1 and dr can be up to ~50
       tolerance.fit.naka_rushton_RB_saturation_index = 0.1; %0 to 1 if R100 > R0, otherwise can theoretically go up to inf
       tolerance.fit.naka_rushton_RB_sensitivity = 10; %units of 1/contrast (1 to up to 1000)
       tolerance.fit.naka_rushton_RBN_parameters = 0.1;
       tolerance.fit.naka_rushton_RBN_contrast = 0.1;
       tolerance.fit.naka_rushton_RBN_values = 0.1;
       tolerance.fit.naka_rushton_RBN_pref = 0.1;
       tolerance.fit.naka_rushton_RBN_empirical_c50 = 0.1;
       tolerance.fit.naka_rushton_RBN_r2 = 0.1;
       tolerance.fit.naka_rushton_RBN_relative_max_gain = 0.1;
       tolerance.fit.naka_rushton_RBN_saturation_index = 0.1;
       tolerance.fit.naka_rushton_RBN_sensitivity = 0.1;
       tolerance.fit.naka_rushton_RBNS_parameters = 0.05;
       tolerance.fit.naka_rushton_RBNS_contrast = 0.1;
       tolerance.fit.naka_rushton_RBNS_values = 0.1;
       tolerance.fit.naka_rushton_RBNS_pref = 0.1; 
       tolerance.fit.naka_rushton_RBNS_empirical_c50 = 0.1;
       tolerance.fit.naka_rushton_RBNS_r2 = 0.1;
       tolerance.fit.naka_rushton_RBNS_relative_max_gain = 0.1;
       tolerance.fit.naka_rushton_RBNS_saturation_index = 0.1;
       tolerance.fit.naka_rushton_RBNS_sensitivity = 0.1;
     
       
    case 'low_noise',
    
       tolerance.tuning_curve.contrast = 1e-6;
       tolerance.tuning_curve.mean = 1e-6;
       tolerance.tuning_curve.stddev = 1e-6;
       tolerance.tuning_curve.stderr = 1e-6;
       tolerance.tuning_curve.individual = 0.5;
       tolerance.tuning_curve.control_stddev = 1e-6;
       tolerance.tuning_curve.control_stderr = 1e-6;
       tolerance.significance.visual_response_anova_p = 0.1;
       tolerance.significance.across_stimuli_anova_p = 0.1;
       tolerance.fitless.interpolated_c50 = 0.1;
       tolerance.fit.naka_rushton_RB_parameters(1) = 5; %max response (0 to ~50)
       tolerance.fit.naka_rushton_RB_parameters(2) = 0.1; %value of contrasts (0 to 1) such that R(C50) = 0.5 * max(R)
       tolerance.fit.naka_rushton_RB_contrast = 0.1; %0 to 1
       tolerance.fit.naka_rushton_RB_values = 5; %responses (0 to ~50)
       tolerance.fit.naka_rushton_RB_pref = 0.1;	%0 to 1
       tolerance.fit.naka_rushton_RB_empirical_c50 = 0.1; %0 to 1
       tolerance.fit.naka_rushton_RB_r2 = 0.1; %0 to 1
       tolerance.fit.naka_rushton_RB_relative_max_gain = 5; %max(dr/dc) where dc is <1 and dr can be up to ~50
       tolerance.fit.naka_rushton_RB_saturation_index = 0.1; %0 to 1 if R100 > R0, otherwise can theoretically go up to inf
       tolerance.fit.naka_rushton_RB_sensitivity = 10; %units of 1/contrast (1 to up to 1000)
       tolerance.fit.naka_rushton_RBN_parameters = 0.1;
       tolerance.fit.naka_rushton_RBN_contrast = 0.1;
       tolerance.fit.naka_rushton_RBN_values = 0.1;
       tolerance.fit.naka_rushton_RBN_pref = 0.1;
       tolerance.fit.naka_rushton_RBN_empirical_c50 = 0.1;
       tolerance.fit.naka_rushton_RBN_r2 = 0.1;
       tolerance.fit.naka_rushton_RBN_relative_max_gain = 0.1;
       tolerance.fit.naka_rushton_RBN_saturation_index = 0.1;
       tolerance.fit.naka_rushton_RBN_sensitivity = 0.1;
       tolerance.fit.naka_rushton_RBNS_parameters = 0.1;
       tolerance.fit.naka_rushton_RBNS_contrast = 0.1;
       tolerance.fit.naka_rushton_RBNS_values = 0.1;
       tolerance.fit.naka_rushton_RBNS_pref = 0.1; 
       tolerance.fit.naka_rushton_RBNS_empirical_c50 = 0.1;
       tolerance.fit.naka_rushton_RBNS_r2 = 0.1;
       tolerance.fit.naka_rushton_RBNS_relative_max_gain = 0.1;
       tolerance.fit.naka_rushton_RBNS_saturation_index = 0.1;
       tolerance.fit.naka_rushton_RBNS_sensitivity = 0.1;
       
    case 'high_noise',
    
       tolerance.tuning_curve.contrast = 1e-6;
       tolerance.tuning_curve.mean = 1e-6;
       tolerance.tuning_curve.stddev = 1e-6;
       tolerance.tuning_curve.stderr = 1e-6;
       tolerance.tuning_curve.individual = 0.5;
       tolerance.tuning_curve.control_stddev = 1e-6;
       tolerance.tuning_curve.control_stderr = 1e-6;
       tolerance.significance.visual_response_anova_p = 0.1;
       tolerance.significance.across_stimuli_anova_p = 0.1;
       tolerance.fitless.interpolated_c50 = 0.1;
       tolerance.fit.naka_rushton_RB_parameters(1) = 5; %max response (0 to ~50)
       tolerance.fit.naka_rushton_RB_parameters(2) = 0.1; %value of contrasts (0 to 1) such that R(C50) = 0.5 * max(R)
       tolerance.fit.naka_rushton_RB_contrast = 0.1; %0 to 1
       tolerance.fit.naka_rushton_RB_values = 5; %responses (0 to ~50)
       tolerance.fit.naka_rushton_RB_pref = 0.1;	%0 to 1
       tolerance.fit.naka_rushton_RB_empirical_c50 = 0.1; %0 to 1
       tolerance.fit.naka_rushton_RB_r2 = 0.1; %0 to 1
       tolerance.fit.naka_rushton_RB_relative_max_gain = 5; %max(dr/dc) where dc is <1 and dr can be up to ~50
       tolerance.fit.naka_rushton_RB_saturation_index = 0.1; %0 to 1 if R100 > R0, otherwise can theoretically go up to inf
       tolerance.fit.naka_rushton_RB_sensitivity = 10; %units of 1/contrast (1 to up to 1000)
       tolerance.fit.naka_rushton_RBN_parameters = 0.1;
       tolerance.fit.naka_rushton_RBN_contrast = 0.1;
       tolerance.fit.naka_rushton_RBN_values = 0.1;
       tolerance.fit.naka_rushton_RBN_pref = 0.1;
       tolerance.fit.naka_rushton_RBN_empirical_c50 = 0.1;
       tolerance.fit.naka_rushton_RBN_r2 = 0.1;
       tolerance.fit.naka_rushton_RBN_relative_max_gain = 0.1;
       tolerance.fit.naka_rushton_RBN_saturation_index = 0.1;
       tolerance.fit.naka_rushton_RBN_sensitivity = 0.1;
       tolerance.fit.naka_rushton_RBNS_parameters = 0.1;
       tolerance.fit.naka_rushton_RBNS_contrast = 0.1;
       tolerance.fit.naka_rushton_RBNS_values = 0.1;
       tolerance.fit.naka_rushton_RBNS_pref = 0.1; 
       tolerance.fit.naka_rushton_RBNS_empirical_c50 = 0.1;
       tolerance.fit.naka_rushton_RBNS_r2 = 0.1;
       tolerance.fit.naka_rushton_RBNS_relative_max_gain = 0.1;
       tolerance.fit.naka_rushton_RBNS_saturation_index = 0.1;
       tolerance.fit.naka_rushton_RBNS_sensitivity = 0.1;
        
    otherwise,
       error(['Unknown scope ' scope '.']);
end;


% start comparison

doc_e = document_expected.document_properties.contrast_tuning;
doc_a = document_actual.document_properties.contrast_tuning;

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
%	contrast                               
%	mean                                   
%	stddev                                 
%	stderr                                
%	individual                             
%	control_stddev                        
%	control_stderr

% tol: 1e -6

[b_(1),errormsg_{1}] = ndi.test.values_within_tolerance(doc_e.tuning_curve.contrast, doc_a.tuning_curve.contrast, tolerance.tuning_curve.contrast, 'contrast');
[b_(2),errormsg_{2}] = ndi.test.values_within_tolerance(doc_e.tuning_curve.mean, doc_a.tuning_curve.mean, tolerance.tuning_curve.mean, 'mean');
[b_(3),errormsg_{3}] = ndi.test.values_within_tolerance(doc_e.tuning_curve.stddev, doc_a.tuning_curve.stddev, tolerance.tuning_curve.stddev, 'stddev');
[b_(4),errormsg_{4}] = ndi.test.values_within_tolerance(doc_e.tuning_curve.stderr, doc_a.tuning_curve.stderr, tolerance.tuning_curve.stderr, 'stderr');
[b_(5),errormsg_{5}] = ndi.test.values_within_tolerance(doc_e.tuning_curve.individual, doc_a.tuning_curve.individual, tolerance.tuning_curve.individual, 'individual');
[b_(6),errormsg_{6}] = ndi.test.values_within_tolerance(doc_e.tuning_curve.control_stddev, doc_a.tuning_curve.control_stddev, tolerance.tuning_curve.control_stddev, 'control_stddev');
[b_(7),errormsg_{7}] = ndi.test.values_within_tolerance(doc_e.tuning_curve.control_stderr, doc_a.tuning_curve.control_stderr, tolerance.tuning_curve.control_stderr, 'control_stderr');

% Comparing Significance
%   visual_response_anova_p
%   across_stimuli_anova_p

[b_(8),errormsg_{8}] = ndi.test.values_within_tolerance(doc_e.significance.visual_response_anova_p, doc_a.significance.visual_response_anova_p, tolerance.significance.visual_response_anova_p, 'visual response anova p');
[b_(9),errormsg_{9}] = ndi.test.values_within_tolerance(doc_e.significance.across_stimuli_anova_p, doc_a.significance.across_stimuli_anova_p, tolerance.significance.across_stimuli_anova_p, 'across stimuli anova p');

% Comparing Fitless
%   interpolated_c50

% tol = 0.1

[b_(10),errormsg_{10}] = ndi.test.values_within_tolerance(doc_e.fitless.interpolated_c50, doc_a.fitless.interpolated_c50, tolerance.fitless.interpolated_c50, 'interpolated c50');


% tol = 0.1

% Comparing Fit
%   naka_rushton_RB_parameters
%   naka_rushton_RB_contrast
%	naka_rushton_RB_values
%	naka_rushton_RB_pref
%	naka_rushton_RB_empirical_c50							
%	naka_rushton_RB_r2
%	naka_rushton_RB_relative_max_gain
%	naka_rushton_RB_saturation_index
%	naka_rushton_RB_sensitivity
%	naka_rushton_RBN_parameters
%	naka_rushton_RBN_contrast
%	naka_rushton_RBN_values
%	naka_rushton_RBN_pref
%	naka_rushton_RBN_empirical_c50
%	naka_rushton_RBN_r2
%	naka_rushton_RBN_relative_max_gain
%	naka_rushton_RBN_saturation_index
%	naka_rushton_RBN_sensitivity
%	naka_rushton_RBNS_parameters
%	naka_rushton_RBNS_contrast
%	naka_rushton_RBNS_values
%	naka_rushton_RBNS_pref
%	naka_rushton_RBNS_empirical_c50     	
%	naka_rushton_RBNS_r2
%	naka_rushton_RBNS_relative_max_gain
%	naka_rushton_RBNS_saturation_index
%	naka_rushton_RBNS_sensitivity
                                                              
[b_(11),errormsg_{11}] = ndi.test.values_within_tolerance(doc_e.fit.naka_rushton_RB_parameters, doc_a.fit.naka_rushton_RB_parameters, tolerance.fit.naka_rushton_RB_parameters, 'naka rushton RB parameters');
[b_(12),errormsg_{12}] = ndi.test.values_within_tolerance(doc_e.fit.naka_rushton_RB_contrast, doc_a.fit.naka_rushton_RB_contrast, tolerance.fit.naka_rushton_RB_contrast, 'naka rushton RB contrast');
[b_(13),errormsg_{13}] = ndi.test.values_within_tolerance(doc_e.fit.naka_rushton_RB_values, doc_a.fit.naka_rushton_RB_values, tolerance.fit.naka_rushton_RB_values, 'naka rushton RB values');
[b_(14),errormsg_{14}] = ndi.test.values_within_tolerance(doc_e.fit.naka_rushton_RB_pref, doc_a.fit.naka_rushton_RB_pref, tolerance.fit.naka_rushton_RB_pref, 'naka rushton RB pref');
[b_(15),errormsg_{15}] = ndi.test.values_within_tolerance(doc_e.fit.naka_rushton_RB_empirical_c50, doc_a.fit.naka_rushton_RB_empirical_c50, tolerance.fit.naka_rushton_RB_empirical_c50, 'naka rushton RB empirical c50');
%[b_(16),errormsg_{16}] = ndi.test.values_within_tolerance(doc_e.fit.naka_rushton_RB_r2, doc_a.fit.naka_rushton_RB_r2, tolerance.fit.naka_rushton_RB_r2, 'naka rushton RB r2');
%[b_(17),errormsg_{17}] = ndi.test.values_within_tolerance(doc_e.fit.naka_rushton_RB_relative_max_gain, doc_a.fit.naka_rushton_RB_relative_max_gain, tolerance.fit.naka_rushton_RB_relative_max_gain, 'naka rushton RB relative max');
[b_(18),errormsg_{18}] = ndi.test.values_within_tolerance(doc_e.fit.naka_rushton_RB_saturation_index, doc_a.fit.naka_rushton_RB_saturation_index, tolerance.fit.naka_rushton_RB_saturation_index, 'naka rushton RB saturation index');
[b_(19),errormsg_{19}] = ndi.test.values_within_tolerance(doc_e.fit.naka_rushton_RB_sensitivity, doc_a.fit.naka_rushton_RB_sensitivity, tolerance.fit.naka_rushton_RB_sensitivity, 'naka rushton RB sensitivity');
[b_(20),errormsg_{20}] = ndi.test.values_within_tolerance(doc_e.fit.naka_rushton_RBN_parameters, doc_a.fit.naka_rushton_RBN_parameters, tolerance.fit.naka_rushton_RBN_parameters, 'naka rushton RBN parameters');
[b_(21),errormsg_{21}] = ndi.test.values_within_tolerance(doc_e.fit.naka_rushton_RBN_contrast, doc_a.fit.naka_rushton_RBN_contrast, tolerance.fit.naka_rushton_RBN_contrast, 'naka rushton RBN contrast');
[b_(22),errormsg_{22}] = ndi.test.values_within_tolerance(doc_e.fit.naka_rushton_RBN_values, doc_a.fit.naka_rushton_RBN_values, tolerance.fit.naka_rushton_RBN_values, 'naka rushton RBN values');
[b_(23),errormsg_{23}] = ndi.test.values_within_tolerance(doc_e.fit.naka_rushton_RBN_pref, doc_a.fit.naka_rushton_RBN_pref, tolerance.fit.naka_rushton_RBN_pref, 'naka rushton RBN pref');
[b_(24),errormsg_{24}] = ndi.test.values_within_tolerance(doc_e.fit.naka_rushton_RBN_empirical_c50, doc_a.fit.naka_rushton_RBN_empirical_c50, tolerance.fit.naka_rushton_RBN_empirical_c50, 'naka rushton RBN empirical c50');
%[b_(25),errormsg_{25}] = ndi.test.values_within_tolerance(doc_e.fit.naka_rushton_RBN_r2, doc_a.fit.naka_rushton_RBN_r2, tolerance.fit.naka_rushton_RBN_r2, 'naka rushton RBN r2');
%[b_(26),errormsg_{26}] = ndi.test.values_within_tolerance(doc_e.fit.naka_rushton_RBN_relative_max_gain, doc_a.fit.naka_rushton_RBN_relative_max_gain, tolerance.fit.naka_rushton_RBN_relative_max_gain, 'naka rushton RBN relative max');
[b_(27),errormsg_{27}] = ndi.test.values_within_tolerance(doc_e.fit.naka_rushton_RBN_saturation_index, doc_a.fit.naka_rushton_RBN_saturation_index, tolerance.fit.naka_rushton_RBN_saturation_index, 'naka rushton RBN saturation index');
[b_(28),errormsg_{28}] = ndi.test.values_within_tolerance(doc_e.fit.naka_rushton_RBN_sensitivity, doc_a.fit.naka_rushton_RBN_sensitivity, tolerance.fit.naka_rushton_RBN_sensitivity, 'naka rushton RBN sensitivity');
[b_(29),errormsg_{29}] = ndi.test.values_within_tolerance(doc_e.fit.naka_rushton_RBNS_parameters, doc_a.fit.naka_rushton_RBNS_parameters, tolerance.fit.naka_rushton_RBNS_parameters, 'naka rushton RBNS parameters');
[b_(30),errormsg_{30}] = ndi.test.values_within_tolerance(doc_e.fit.naka_rushton_RBNS_contrast, doc_a.fit.naka_rushton_RBNS_contrast, tolerance.fit.naka_rushton_RBNS_contrast, 'naka rushton RBNS contrast');
[b_(31),errormsg_{31}] = ndi.test.values_within_tolerance(doc_e.fit.naka_rushton_RBNS_values, doc_a.fit.naka_rushton_RBNS_values, tolerance.fit.naka_rushton_RBNS_values, 'naka rushton RBNS values');
[b_(32),errormsg_{32}] = ndi.test.values_within_tolerance(doc_e.fit.naka_rushton_RBNS_pref, doc_a.fit.naka_rushton_RBNS_pref, tolerance.fit.naka_rushton_RBNS_pref, 'naka rushton RBNS pref');
[b_(33),errormsg_{33}] = ndi.test.values_within_tolerance(doc_e.fit.naka_rushton_RBNS_empirical_c50, doc_a.fit.naka_rushton_RBNS_empirical_c50, tolerance.fit.naka_rushton_RBNS_empirical_c50, 'naka rushton RBNS empirical c50');
%[b_(34),errormsg_{34}] = ndi.test.values_within_tolerance(doc_e.fit.naka_rushton_RBNS_r2, doc_a.fit.naka_rushton_RBNS_r2, tolerance.fit.naka_rushton_RBNS_r2, 'naka rushton RBNS r2');
%[b_(35),errormsg_{35}] = ndi.test.values_within_tolerance(doc_e.fit.naka_rushton_RBNS_relative_max_gain, doc_a.fit.naka_rushton_RBNS_relative_max_gain, tolerance.fit.naka_rushton_RBNS_relative_max_gain, 'naka rushton RBNS relative max');
[b_(36),errormsg_{36}] = ndi.test.values_within_tolerance(doc_e.fit.naka_rushton_RBNS_saturation_index, doc_a.fit.naka_rushton_RBNS_saturation_index, tolerance.fit.naka_rushton_RBNS_saturation_index, 'naka rushton RBNS saturation index');
[b_(37),errormsg_{37}] = ndi.test.values_within_tolerance(doc_e.fit.naka_rushton_RBNS_sensitivity, doc_a.fit.naka_rushton_RBNS_sensitivity, tolerance.fit.naka_rushton_RBNS_sensitivity, 'naka rushton RBNS sensitivity');


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


