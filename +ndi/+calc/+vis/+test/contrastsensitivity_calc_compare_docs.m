function [b, errormsg] = contrastsensitivity_calc_compare_docs(document_expected, document_actual, scope)
% contrastsensitivity_calc_compare_docs
%
% [B, ERRORMSG] = ndi.calc.vis.test.contrastsensitivity_calc_compare_docs(DOC_EXPECTED, DOC_ACTUAL, SCOPE)
%


b = 1;
errormsg = '';

% start comparison

doc_e = document_expected.document_properties.contrastsensitivity_calc;
doc_a = document_actual.document_properties.contrastsensitivity_calc;

% Comparing
%   spatial_frequencies
%   sensitivity_RB
%   sensitivity_RBN
%   sensitivity_RBNS
%   across_stimuli_varies_p_bonferroni
%   response_varies_p_bonferroni
%   visual_response_p_bonferroni
%   is_modulated_response
%   response_type

[b,errormsg] = ndi.test.values_within_tolerance(doc_e.spatial_frequencies, doc_a.spatial_frequencies, tolerance, 'spatial frequencies');
[b,errormsg] = ndi.test.values_within_tolerance(doc_e.sensitivity_RB, doc_a.sensitivity_RB, tolerance, 'sensitivity RB');
[b,errormsg] = ndi.test.values_within_tolerance(doc_e.sensitivity_RBN, doc_a.sensitivity_RBN, tolerance, 'sensitivity_RBN');
[b,errormsg] = ndi.test.values_within_tolerance(doc_e.sensitivity_RBNS, doc_a.sensitivity_RBNS, tolerance, 'sensitivity_RBNS');
[b,errormsg] = ndi.test.values_within_tolerance(doc_e.across_stimuli_varies_p_bonferroni, doc_a.across_stimuli_varies_p_bonferroni, tolerance, 'across stimuli varies p bonferroni');
[b,errormsg] = ndi.test.values_within_tolerance(doc_e.response_varies_p_bonferroni, doc_a.response_varies_p_bonferroni, tolerance, 'response varies p bonferroni');
[b,errormsg] = ndi.test.values_within_tolerance(doc_e.visual_response_p_bonferroni, doc_a.visual_response_p_bonferroni, tolerance, 'visual response p bonferroni');

is_modulated_response_match = strcmpi(doc_e.is_modulated_response, doc_a.is_modulated_response);
if ~is_modulated_response_match
   b = 0;
   errormsg = ['Expected is modulated response is ' doc_e.is_modulated_response ' but observed ' doc_a.is_modulated_response];
   return;
end

response_type_match = strcmpi(doc_e.response_type, doc_a.response_type);
if ~response_type_match
   b = 0;
   errormsg = ['Expected response type of ' doc_e.response_type ' but observed ' doc_a.response_type];
   return;
end


end
