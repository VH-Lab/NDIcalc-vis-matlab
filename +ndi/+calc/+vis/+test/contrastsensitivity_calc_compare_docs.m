function [b_, errormsg_] = contrastsensitivity_calc_compare_docs(document_expected, document_actual, scope)
% contrastsensitivity_calc_compare_docs
%
% [B, ERRORMSG] = ndi.calc.vis.test.contrastsensitivity_calc_compare_docs(DOC_EXPECTED, DOC_ACTUAL, SCOPE)
%

% Initialize b_ as a row vector of ones for tracking comparison results
% Initialize errormsg_ as an empty cell array to hold error messages

b_ = ones(1,9);
errormsg_ = cell(1,9);

% establish scope-dependent tolerances
switch(scope),
    case 'standard',

       tolerance.spatial_frequencies = 0.1;
       tolerance.sensitivity_RB = 0.1;
       tolerance.sensitivity_RBN = 0.1;
       tolerance.sensitivity_RBNS = 0.1;
       tolerance.across_stimuli_varies_p_bonferroni = 0.1;
       tolerance.response_varies_p_bonferroni = 0.1;
       tolerance.visual_response_p_bonferroni = 0.1;

    case 'low_noise',

       tolerance.spatial_frequencies = 0.1;
       tolerance.sensitivity_RB = 0.1;
       tolerance.sensitivity_RBN = 0.1;
       tolerance.sensitivity_RBNS = 0.1;
       tolerance.across_stimuli_varies_p_bonferroni = 0.1;
       tolerance.response_varies_p_bonferroni = 0.1;
       tolerance.visual_response_p_bonferroni = 0.1;

    case 'high_noise',

       tolerance.spatial_frequencies = 0.1;
       tolerance.sensitivity_RB = 0.1;
       tolerance.sensitivity_RBN = 0.1;
       tolerance.sensitivity_RBNS = 0.1;
       tolerance.across_stimuli_varies_p_bonferroni = 0.1;
       tolerance.response_varies_p_bonferroni = 0.1;
       tolerance.visual_response_p_bonferroni = 0.1;

    otherwise,
    
       error(['Unknown scope ' scope '.']);
 end;

% start comparison

doc_e = document_expected.document_properties.contrastsensitivity_calc;
doc_a = document_actual.document_properties.contrastsensitivity_calc;

% tol = 0.1

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

[b_(1),errormsg_{1}] = ndi.test.values_within_tolerance(doc_e.spatial_frequencies, doc_a.spatial_frequencies, tolerance.spatial_frequencies, 'spatial frequencies');
[b_(2),errormsg_{2}] = ndi.test.values_within_tolerance(doc_e.sensitivity_RB, doc_a.sensitivity_RB, tolerance.sensitivity_RB, 'sensitivity RB');
[b_(3),errormsg_{3}] = ndi.test.values_within_tolerance(doc_e.sensitivity_RBN, doc_a.sensitivity_RBN, tolerance.sensitivity_RBN, 'sensitivity_RBN');
[b_(4),errormsg_{4}] = ndi.test.values_within_tolerance(doc_e.sensitivity_RBNS, doc_a.sensitivity_RBNS, tolerance.sensitivity_RBNS, 'sensitivity_RBNS');
%[b_(5),errormsg_{5}] = ndi.test.values_within_tolerance(doc_e.across_stimuli_varies_p_bonferroni, doc_a.across_stimuli_varies_p_bonferroni, tolerance.across_stimuli_varies_p_bonferroni, 'across stimuli varies p bonferroni');
[b_(6),errormsg_{6}] = ndi.test.values_within_tolerance(doc_e.response_varies_p_bonferroni, doc_a.response_varies_p_bonferroni, tolerance.response_varies_p_bonferroni, 'response varies p bonferroni');
[b_(7),errormsg_{7}] = ndi.test.values_within_tolerance(doc_e.visual_response_p_bonferroni, doc_a.visual_response_p_bonferroni, tolerance.visual_response_p_bonferroni, 'visual response p bonferroni');

if doc_e.is_modulated_response ~= doc_a.is_modulated_response;
   b_(8) = 0;
   errormsg_{8} = ['Is modulated response do not match.'];
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
if all(b_) && all(cellfun('isempty', errormsg_))
    disp('No error');
end


end
