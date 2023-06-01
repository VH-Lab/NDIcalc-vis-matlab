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

%Comparing spatial_frequencies

tol = %put in appropriate tolorance;

if any(abs(doc_e.spatial_frequencies(:) - doc_a.spatial_frequencies(:)) > tol)
   b = 0;
   errormsg = ['Spatial frequencies differ by greater than ' num2str(tol) '.'];
end

%Comparing sensitivity_RB

tol = %put in appropriate tolorance;

if any(abs(doc_e.sensitivity_RB(:) - doc_a.sensitivity_RB(:)) > tol)
   b = 0;
   errormsg = ['Sensitivity RB differ by greater than ' num2str(tol) '.'];
end

%Comparing sensitivity_RBN

tol = %put in appropriate tolorance;

if any(abs(doc_e.sensitivity_RBN(:) - doc_a.sensitivity_RBN(:)) > tol)
   b = 0;
   errormsg = ['Sensitivity RBN differ by greater than ' num2str(tol) '.'];
end

%Comparing sensitivity_RBNS

tol = %put in appropriate tolorance;

if any(abs(doc_e.sensitivity_RBNS(:) - doc_a.sensitivity_RBNS(:)) > tol)
   b = 0;
   errormsg = ['Sensitivity RBNS differ by greater than ' num2str(tol) '.'];
end

% Comparing is_modulated_response

tol = %put in appropriate tolorance;

is_modulated_response_match = vlt.data.sizeeq(doc_e.is_modulated_response(:),doc_a.is_modulated_response(:));
if is_modulated_response_match
   is_modulated_response_match = max(abs(doc_e.is_modulated_response(:) - doc_a.is_modulated_response(:))) < tol;
end

if ~is_modulated_response_match
   b = 0;
   errormsg = ['Is modulated response is differed by greater than ' num2str(tol) '.'];
   return;
end

%Comparing across_stimuli_varies_p_bonferroni

tol = %put in appropriate tolorance;

if any(abs(doc_e.across_stimuli_varies_p_bonferroni(:) - doc_a.across_stimuli_varies_p_bonferroni(:)) > tol)
   b = 0;
   errormsg = ['Across Stimuli Varies p bonferroni differ by greater than ' num2str(tol) '.'];
end

%Comparing response_varies_p_bonferroni

tol = %put in appropriate tolorance;

if any(abs(doc_e.response_varies_p_bonferroni(:) - doc_a.response_varies_p_bonferroni(:)) > tol)
   b = 0;
   errormsg = ['Response varies p bonferroni differ by greater than ' num2str(tol) '.'];
end

%Comparing visual_response_p_bonferroni

tol = %put in appropriate tolorance;

if any(abs(doc_e.visual_response_p_bonferroni(:) - doc_a.visual_response_p_bonferroni(:)) > tol)
   b = 0;
   errormsg = ['Visual response p bonferroni differ by greater than ' num2str(tol) '.'];
end

%Comparing response_type
% check to make sure not NaN

if isnan(doc_a.response_type)
    disp('Response type is not a number')
end

if ~isnan(doc_a.response_type)

tol = %put in appropriate tolorance;

if any(abs(doc_e.response_type(:) - doc_a.response_type(:)) > tol)
   b = 0;
   errormsg = ['Response type differ by greater than ' num2str(tol) '.'];
end

end

end
