function [b, errormsg] = contrast_tuning_compare_docs(document_expected, document_actual, scope)
% Comparing Contrast_tuning
%
% [B, ERRORMSG] = ndi.calc.vis.test.contrast_tuning_compare_docs(DOC_EXPECTED, DOC_ACTUAL, SCOPE)
%
%


b = 1;
errormsg = '';

% start comparison

doc_e = document_expected.document_properties.contrast_tuning;
doc_a = document_actual.document_properties.contrast_tuning;

% Comparing properties
%   Response Units

properties_match = strcmpi(char(doc_e.properties.response_units), char(doc_a.properties.response_units));
if ~properties_match
   b = 0;
   errormsg = ['Expected response units in ' doc_e.properties.response_units ' but observed ' doc_a.properties.response_units];
   return;
end

%   Response Type

properties_match = strcmpi(char(doc_e.properties.response_type), char(doc_a.properties.response_type));
if ~properties_match
   b = 0;
   errormsg = ['Expected response type of ' doc_e.properties.response_units ' but observed ' doc_a.properties.response_units];
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

errormsg_list = {};
b_list = [];

[b_list(1),errormsg_list{1}] = ndi.test.values_within_tolerance(doc_e.tuning_curve.contrast, doc_a.tuning_curve.contrast, tolerance, 'contrast');
[b,errormsg] = ndi.test.values_within_tolerance(doc_e.tuning_curve.mean, doc_a.tuning_curve.mean, tolerance, 'mean');
[b,errormsg] = ndi.test.values_within_tolerance(doc_e.tuning_curve.stddev, doc_a.tuning_curve.stddev, tolerance, 'stddev');
[b,errormsg] = ndi.test.values_within_tolerance(doc_e.tuning_curve.stderr, doc_a.tuning_curve.stderr, tolerance, 'stderr');
[b,errormsg] = ndi.test.values_within_tolerance(doc_e.tuning_curve.individual, doc_a.tuning_curve.individual, tolerance, 'individual');
[b,errormsg] = ndi.test.values_within_tolerance(doc_e.tuning_curve.control_stddev, doc_a.tuning_curve.control_stddev, tolerance, 'control_stddev');
[b,errormsg] = ndi.test.values_within_tolerance(doc_e.tuning_curve.control_stderr, doc_a.tuning_curve.control_stderr, tolerance, 'control_stderr');

% Comparing Significance
%   visual_response_anova_p
%   across_stimuli_anova_p

[b,errormsg] = ndi.test.values_within_tolerance(doc_e.significance.visual_response_anova_p, doc_a.significance.visual_response_anova_p, tolerance, 'visual response anova p');
[b,errormsg] = ndi.test.values_within_tolerance(doc_e.significance.across_stimuli_anova_p, doc_a.significance.across_stimuli_anova_p, tolerance, 'across stimuli anova p');

% Comparing Fitless
%   interpolated_c50

[b,errormsg] = ndi.test.values_within_tolerance(doc_e.fitless.interpolated_c50, doc_a.fitless.interpolated_c50, tolerance, 'interpolated c50');

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
                                                                    
[b,errormsg] = ndi.test.values_within_tolerance(doc_e.fit.naka_rushton_RB_parameters, doc_a.fit.naka_rushton_RB_parameters, tolerance, 'naka rushton RB parameters');
[b,errormsg] = ndi.test.values_within_tolerance(doc_e.fit.naka_rushton_RB_contrast, doc_a.fit.naka_rushton_RB_contrast, tolerance, 'naka rushton RB contrast');
[b,errormsg] = ndi.test.values_within_tolerance(doc_e.fit.naka_rushton_RB_values, doc_a.fit.naka_rushton_RB_values, tolerance, 'naka rushton RB values');
[b,errormsg] = ndi.test.values_within_tolerance(doc_e.fit.naka_rushton_RB_pref, doc_a.fit.naka_rushton_RB_pref, tolerance, 'naka rushton RB pref');
[b,errormsg] = ndi.test.values_within_tolerance(doc_e.fit.naka_rushton_RB_empirical_c50, doc_a.fit.naka_rushton_RB_empirical_c50, tolerance, 'naka rushton RB empirical c50');
[b,errormsg] = ndi.test.values_within_tolerance(doc_e.fit.naka_rushton_RB_r2, doc_a.fit.naka_rushton_RB_r2, tolerance, 'naka rushton RB r2');
[b,errormsg] = ndi.test.values_within_tolerance(doc_e.fit.naka_rushton_RB_relative_max_gain, doc_a.fit.naka_rushton_RB_relative_max, tolerance, 'naka rushton RB relative max');
[b,errormsg] = ndi.test.values_within_tolerance(doc_e.fit.naka_rushton_RB_saturation_index, doc_a.fit.naka_rushton_RB_saturation_index, tolerance, 'naka rushton RB saturation index');
[b,errormsg] = ndi.test.values_within_tolerance(doc_e.fit.naka_rushton_RB_sensitivity, doc_a.fit.naka_rushton_RB_sensitivity, tolerance, 'naka rushton RB sensitivity');
[b,errormsg] = ndi.test.values_within_tolerance(doc_e.fit.naka_rushton_RBN_parameters, doc_a.fit.naka_rushton_RBN_parameters, tolerance, 'naka rushton RBN parameters');
[b,errormsg] = ndi.test.values_within_tolerance(doc_e.fit.naka_rushton_RBN_contrast, doc_a.fit.naka_rushton_RBN_contrast, tolerance, 'naka rushton RBN contrast');
[b,errormsg] = ndi.test.values_within_tolerance(doc_e.fit.naka_rushton_RBN_values, doc_a.fit.naka_rushton_RBN_values, tolerance, 'naka rushton RBN values');
[b,errormsg] = ndi.test.values_within_tolerance(doc_e.fit.naka_rushton_RBN_pref, doc_a.fit.naka_rushton_RBN_pref, tolerance, 'naka rushton RBN pref');
[b,errormsg] = ndi.test.values_within_tolerance(doc_e.fit.naka_rushton_RBN_empirical_c50, doc_a.fit.naka_rushton_RBN_empirical_c50, tolerance, 'naka rushton RBN empirical c50');
[b,errormsg] = ndi.test.values_within_tolerance(doc_e.fit.naka_rushton_RBN_r2, doc_a.fit.naka_rushton_RBN_r2, tolerance, 'naka rushton RBN r2');
[b,errormsg] = ndi.test.values_within_tolerance(doc_e.fit.naka_rushton_RBN_relative_max_gain, doc_a.fit.naka_rushton_RBN_relative_max, tolerance, 'naka rushton RBN relative max');
[b,errormsg] = ndi.test.values_within_tolerance(doc_e.fit.naka_rushton_RBN_saturation_index, doc_a.fit.naka_rushton_RBN_saturation_index, tolerance, 'naka rushton RBN saturation index');
[b,errormsg] = ndi.test.values_within_tolerance(doc_e.fit.naka_rushton_RBN_sensitivity, doc_a.fit.naka_rushton_RBN_sensitivity, tolerance, 'naka rushton RBN sensitivity');
[b,errormsg] = ndi.test.values_within_tolerance(doc_e.fit.naka_rushton_RBNS_parameters, doc_a.fit.naka_rushton_RBNS_parameters, tolerance, 'naka rushton RBNS parameters');
[b,errormsg] = ndi.test.values_within_tolerance(doc_e.fit.naka_rushton_RBNS_contrast, doc_a.fit.naka_rushton_RBNS_contrast, tolerance, 'naka rushton RBNS contrast');
[b,errormsg] = ndi.test.values_within_tolerance(doc_e.fit.naka_rushton_RBNS_values, doc_a.fit.naka_rushton_RBNS_values, tolerance, 'naka rushton RBNS values');
[b,errormsg] = ndi.test.values_within_tolerance(doc_e.fit.naka_rushton_RBNS_pref, doc_a.fit.naka_rushton_RBNS_pref, tolerance, 'naka rushton RBNS pref');
[b,errormsg] = ndi.test.values_within_tolerance(doc_e.fit.naka_rushton_RBNS_empirical_c50, doc_a.fit.naka_rushton_RBNS_empirical_c50, tolerance, 'naka rushton RBNS empirical c50');
[b,errormsg] = ndi.test.values_within_tolerance(doc_e.fit.naka_rushton_RBNS_r2, doc_a.fit.naka_rushton_RBNS_r2, tolerance, 'naka rushton RBNS r2');
[b,errormsg] = ndi.test.values_within_tolerance(doc_e.fit.naka_rushton_RBNS_relative_max_gain, doc_a.fit.naka_rushton_RBNS_relative_max, tolerance, 'naka rushton RBNS relative max');
[b,errormsg] = ndi.test.values_within_tolerance(doc_e.fit.naka_rushton_RBNS_saturation_index, doc_a.fit.naka_rushton_RBNS_saturation_index, tolerance, 'naka rushton RBNS saturation index');
[b,errormsg] = ndi.test.values_within_tolerance(doc_e.fit.naka_rushton_RBNS_sensitivity, doc_a.fit.naka_rushton_RBNS_sensitivity, tolerance, 'naka rushton RBNS sensitivity');


if ~any(b_list),
    b_failed_index = find(b_list==0);
    b = b_list(b_failed_index); % 0
    errormsg = errormsg_list{b_failed_index};
end;


