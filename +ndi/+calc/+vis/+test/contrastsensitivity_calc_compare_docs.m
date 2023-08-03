function [b_, errormsg_] = contrastsensitivity_calc_compare_docs(document_expected, document_actual, scope)
% contrastsensitivity_calc_compare_docs
%
% [B, ERRORMSG] = ndi.calc.vis.test.contrastsensitivity_calc_compare_docs(DOC_EXPECTED, DOC_ACTUAL, SCOPE)
%

% Initialize b_ as a row vector of ones for tracking comparison results
% Initialize errormsg_ as an empty cell array to hold error messages

b_ = ones(1,9);
errormsg_ = cell(1,9);

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

[b_(1),errormsg_{1}] = ndi.test.values_within_tolerance(doc_e.spatial_frequencies, doc_a.spatial_frequencies, tolerance, 'spatial frequencies');
[b_(2),errormsg_{2}] = ndi.test.values_within_tolerance(doc_e.sensitivity_RB, doc_a.sensitivity_RB, tolerance, 'sensitivity RB');
[b_(3),errormsg_{3}] = ndi.test.values_within_tolerance(doc_e.sensitivity_RBN, doc_a.sensitivity_RBN, tolerance, 'sensitivity_RBN');
[b_(4),errormsg_{4}] = ndi.test.values_within_tolerance(doc_e.sensitivity_RBNS, doc_a.sensitivity_RBNS, tolerance, 'sensitivity_RBNS');
[b_(5),errormsg_{5}] = ndi.test.values_within_tolerance(doc_e.across_stimuli_varies_p_bonferroni, doc_a.across_stimuli_varies_p_bonferroni, tolerance, 'across stimuli varies p bonferroni');
[b_(6),errormsg_{6}] = ndi.test.values_within_tolerance(doc_e.response_varies_p_bonferroni, doc_a.response_varies_p_bonferroni, tolerance, 'response varies p bonferroni');
[b_(7),errormsg_{7}] = ndi.test.values_within_tolerance(doc_e.visual_response_p_bonferroni, doc_a.visual_response_p_bonferroni, tolerance, 'visual response p bonferroni');

is_modulated_response_match = strcmpi(char(doc_e.is_modulated_response), char(doc_a.is_modulated_response));
if ~is_modulated_response_match
   b_(8) = 0;
   errormsg_{8} = ['Expected is modulated response is ' doc_e.is_modulated_response ' but observed ' doc_a.is_modulated_response];
   return;
end

response_type_match = strcmpi(char(doc_e.response_type), char(doc_a.response_type));
if ~response_type_match
   b_(9) = 0;
   errormsg_{9} = ['Expected response type of ' doc_e.response_type ' but observed ' doc_a.response_type];
   return;
end

% Identify the b_ values with unmatched results
% Update b_ to only include those
% Update the corresponding errormsg_ messages

if any(b_==0),
    error_indices = find(b_==0);
    b_ = b_(error_indices);
    errormsg_ = errormsg_(error_indices);
end

end
