function [b, errormsg] = contrast_tuning_compare_docs(document_expected, document_actual, scope)
% Comparing Contrast_tuning
%
% [B, ERRORMSG] = ndi.calc.vis.test.contrast_tuning_compare_docs(DOC_EXPECTED, DOC_ACTUAL, SCOPE)
%
%


b = 1;
errormsg = '';

% start comparison for fit

doc_e = document_expected.document_properties.contrast_tuning.fit;
doc_a = document_actual.document_properties.contrast_tuning.fit;

% Comparing naka_rushton_RB_parameters

tol = %put in appropriate tolorance;

naka_rushton_RB_parameters_match = vlt.data.sizeeq(doc_e.naka_rushton_RB_parameters(:),doc_a.naka_rushton_RB_parameters(:));
if naka_rushton_RB_parameters_match
   naka_rushton_RB_parameters_match = max(abs(doc_e.naka_rushton_RB_parameters(:) - doc_a.naka_rushton_RB_parameters(:))) < tol;
end

if ~naka_rushton_RB_parameters_match
   b = 0;
   errormsg = ['Naka Rushton RB parameters differed by greater than ' num2str(tol) '.'];
   return;
end

% Comparing naka_rushton_RB_contrast

tol = %put in appropriate tolorance;

naka_rushton_RB_contrast_match = vlt.data.sizeeq(doc_e.naka_rushton_RB_contrast(:),doc_a.naka_rushton_RB_contrast(:));
if naka_rushton_RB_contrast_match
   naka_rushton_RB_contrast_match = max(abs(doc_e.naka_rushton_RB_contrast(:) - doc_a.naka_rushton_RB_contrast(:))) < tol;
end

if ~naka_rushton_RB_contrast_match
   b = 0;
   errormsg = ['Naka Rushton RB comtrast differed by greater than ' num2str(tol) '.'];
   return;
end

% Comparing naka_rushton_RB_values

tol = %put in appropriate tolorance;

naka_rushton_RB_values_match = vlt.data.sizeeq(doc_e.naka_rushton_RB_values(:),doc_a.naka_rushton_RB_values(:));
if naka_rushton_RB_values_match
   naka_rushton_RB_values_match = max(abs(doc_e.naka_rushton_RB_values(:) - doc_a.naka_rushton_RB_values(:))) < tol;
end

if ~naka_rushton_RB_values_match
   b = 0;
   errormsg = ['Naka Rushton RB values differed by greater than ' num2str(tol) '.'];
   return;
end

% Comparing naka_rushton_RB_pref

tol = %put in appropriate tolorance;

naka_rushton_RB_pref_match = vlt.data.sizeeq(doc_e.naka_rushton_RB_pref(:),doc_a.naka_rushton_RB_pref(:));
if naka_rushton_RB_pref_match
   naka_rushton_RB_pref_match = max(abs(doc_e.naka_rushton_RB_pref(:) - doc_a.naka_rushton_RB_pref(:))) < tol;
end

if ~naka_rushton_RB_pref_match
   b = 0;
   errormsg = ['Naka Rushton RB pref differed by greater than ' num2str(tol) '.'];
   return;
end

% Comparing naka_rushton_RB_empirical_c50

tol = %put in appropriate tolorance;

naka_rushton_RB_empirical_c50_match = vlt.data.sizeeq(doc_e.naka_rushton_RB_empirical_c50(:),doc_a.naka_rushton_RB_empirical_c50(:));
if naka_rushton_RB_empirical_c50_match
   naka_rushton_RB_empirical_c50_match = max(abs(doc_e.naka_rushton_RB_empirical_c50(:) - doc_a.naka_rushton_RB_empirical_c50(:))) < tol;
end

if ~naka_rushton_RB_empirical_c50_match
   b = 0;
   errormsg = ['Naka Rushton RB empirical c50 differed by greater than ' num2str(tol) '.'];
   return;
end

% Comparing naka_rushton_RB_r2

tol = %put in appropriate tolorance;

naka_rushton_RB_r2_match = vlt.data.sizeeq(doc_e.naka_rushton_RB_r2(:),doc_a.naka_rushton_RB_r2(:));
if naka_rushton_RB_r2_match
   naka_rushton_RB_r2_match = max(abs(doc_e.naka_rushton_RB_r2(:) - doc_a.naka_rushton_RB_r2(:))) < tol;
end

if ~naka_rushton_RB_r2_match
   b = 0;
   errormsg = ['Naka Rushton RB r2 differed by greater than ' num2str(tol) '.'];
   return;
end

% Comparing naka_rushton_RB_relative_max_gain

tol = %put in appropriate tolorance;

naka_rushton_RB_relative_max_gain_match = vlt.data.sizeeq(doc_e.naka_rushton_RB_relative_max_gain(:),doc_a.naka_rushton_RB_relative_max_gain(:));
if naka_rushton_RB_relative_max_gain_match
   naka_rushton_RB_relative_max_gain_match = max(abs(doc_e.naka_rushton_RB_relative_max_gain(:) - doc_a.naka_rushton_RB_relative_max_gain(:))) < tol;
end

if ~naka_rushton_RB_relative_max_gain_match
   b = 0;
   errormsg = ['Naka Rushton RB relative max gain differed by greater than ' num2str(tol) '.'];
   return;
end

% Comparing naka_rushton_RB_saturation_index

tol = %put in appropriate tolorance;

naka_rushton_RB_saturation_index_match = vlt.data.sizeeq(doc_e.naka_rushton_RB_saturation_index(:),doc_a.naka_rushton_RB_saturation_index(:));
if naka_rushton_RB_saturation_index_match
   naka_rushton_RB_saturation_index_match = max(abs(doc_e.naka_rushton_RB_saturation_index(:) - doc_a.naka_rushton_RB_saturation_index(:))) < tol;
end

if ~naka_rushton_RB_saturation_index_match
   b = 0;
   errormsg = ['Naka Rushton RB saturation index differed by greater than ' num2str(tol) '.'];
   return;
end

% Comparing naka_rushton_RB_sensitivity

tol = %put in appropriate tolorance;

naka_rushton_RB_sensitivity_match = vlt.data.sizeeq(doc_e.naka_rushton_RB_sensitivity(:),doc_a.naka_rushton_RB_sensitivity(:));
if naka_rushton_RB_sensitivity_match
   naka_rushton_RB_sensitivity_match = max(abs(doc_e.naka_rushton_RB_sensitivity(:) - doc_a.naka_rushton_RB_sensitivity(:))) < tol;
end

if ~naka_rushton_RB_sensitivity_match
   b = 0;
   errormsg = ['Naka Rushton RB sensitivity differed by greater than ' num2str(tol) '.'];
   return;
end

% Comparing naka_rushton_RBN_parameters

tol = %put in appropriate tolorance;

naka_rushton_RBN_parameters_match = vlt.data.sizeeq(doc_e.naka_rushton_RBN_parameters(:),doc_a.naka_rushton_RBN_parameters(:));
if naka_rushton_RBN_parameters_match
   naka_rushton_RBN_parameters_match = max(abs(doc_e.naka_rushton_RBN_parameters(:) - doc_a.naka_rushton_RBN_parameters(:))) < tol;
end

if ~naka_rushton_RBN_parameters_match
   b = 0;
   errormsg = ['Naka Rushton RBN parameters differed by greater than ' num2str(tol) '.'];
   return;
end

% Comparing naka_rushton_RBN_contrast

tol = %put in appropriate tolorance;

naka_rushton_RBN_contrast_match = vlt.data.sizeeq(doc_e.naka_rushton_RBN_contrast(:),doc_a.naka_rushton_RBN_contrast(:));
if naka_rushton_RBN_contrast_match
   naka_rushton_RBN_contrast_match = max(abs(doc_e.naka_rushton_RBN_contrast(:) - doc_a.naka_rushton_RBN_contrast(:))) < tol;
end

if ~naka_rushton_RBN_contrast_match
   b = 0;
   errormsg = ['Naka Rushton RBN contrast differed by greater than ' num2str(tol) '.'];
   return;
end

% Comparing naka_rushton_RBN_values

tol = %put in appropriate tolorance;

naka_rushton_RBN_values_match = vlt.data.sizeeq(doc_e.naka_rushton_RBN_values(:),doc_a.naka_rushton_RBN_values(:));
if naka_rushton_RBN_values_match
   naka_rushton_RBN_values_match = max(abs(doc_e.naka_rushton_RBN_values(:) - doc_a.naka_rushton_RBN_values(:))) < tol;
end

if ~naka_rushton_RBN_values_match
   b = 0;
   errormsg = ['Naka Rushton RBN values differed by greater than ' num2str(tol) '.'];
   return;
end

% Comparing naka_rushton_RBN_pref

tol = %put in appropriate tolorance;

naka_rushton_RBN_pref_match = vlt.data.sizeeq(doc_e.naka_rushton_RBN_pref(:),doc_a.naka_rushton_RBN_pref(:));
if naka_rushton_RBN_pref_match
   naka_rushton_RBN_pref_match = max(abs(doc_e.naka_rushton_RBN_pref(:) - doc_a.naka_rushton_RBN_pref(:))) < tol;
end

if ~naka_rushton_RBN_pref_match
   b = 0;
   errormsg = ['Naka Rushton RBN pref differed by greater than ' num2str(tol) '.'];
   return;
end

% Comparing naka_rushton_RBN_empirical_c50

tol = %put in appropriate tolorance;

naka_rushton_RBN_empirical_c50_match = vlt.data.sizeeq(doc_e.naka_rushton_RBN_empirical_c50(:),doc_a.naka_rushton_RBN_empirical_c50(:));
if naka_rushton_RBN_empirical_c50_match
   naka_rushton_RBN_empirical_c50_match = max(abs(doc_e.naka_rushton_RBN_empirical_c50(:) - doc_a.naka_rushton_RBN_empirical_c50(:))) < tol;
end

if ~naka_rushton_RBN_empirical_c50_match
   b = 0;
   errormsg = ['Naka Rushton RBN empirical c50 differed by greater than ' num2str(tol) '.'];
   return;
end

% Comparing naka_rushton_RBN_r2

tol = %put in appropriate tolorance;

naka_rushton_RBN_r2_match = vlt.data.sizeeq(doc_e.naka_rushton_RBN_r2(:),doc_a.naka_rushton_RBN_r2(:));
if naka_rushton_RBN_r2_match
   naka_rushton_RBN_r2_match = max(abs(doc_e.naka_rushton_RBN_r2(:) - doc_a.naka_rushton_RBN_r2(:))) < tol;
end

if ~naka_rushton_RBN_r2_match
   b = 0;
   errormsg = ['Naka Rushton RBN r2 differed by greater than ' num2str(tol) '.'];
   return;
end

% Comparing naka_rushton_RBN_relative_max_gain

tol = %put in appropriate tolorance;

naka_rushton_RBN_relative_max_gain_match = vlt.data.sizeeq(doc_e.naka_rushton_RBN_relative_max_gain(:),doc_a.naka_rushton_RBN_relative_max_gain(:));
if naka_rushton_RBN_relative_max_gain_match
   naka_rushton_RBN_relative_max_gain_match = max(abs(doc_e.naka_rushton_RBN_relative_max_gain(:) - doc_a.naka_rushton_RBN_relative_max_gain(:))) < tol;
end

if ~naka_rushton_RBN_relative_max_gain_match
   b = 0;
   errormsg = ['Naka Rushton RBN relative max gain differed by greater than ' num2str(tol) '.'];
   return;
end

% Comparing naka_rushton_RBN_saturation_index

tol = %put in appropriate tolorance;

naka_rushton_RBN_saturation_index_match = vlt.data.sizeeq(doc_e.naka_rushton_RBN_saturation_index(:),doc_a.naka_rushton_RBN_saturation_index(:));
if naka_rushton_RBN_saturation_index_match
   naka_rushton_RBN_saturation_index_match = max(abs(doc_e.naka_rushton_RBN_saturation_index(:) - doc_a.naka_rushton_RBN_saturation_index(:))) < tol;
end

if ~naka_rushton_RBN_saturation_index_match
   b = 0;
   errormsg = ['Naka Rushton RBN saturation index differed by greater than ' num2str(tol) '.'];
   return;
end

% Comparing naka_rushton_RBN_sensitivity

tol = %put in appropriate tolorance;

naka_rushton_RBN_sensitivity_match = vlt.data.sizeeq(doc_e.naka_rushton_RBN_sensitivity(:),doc_a.naka_rushton_RBN_sensitivity(:));
if naka_rushton_RBN_sensitivity_match
   naka_rushton_RBN_sensitivity_match = max(abs(doc_e.naka_rushton_RBN_sensitivity(:) - doc_a.naka_rushton_RBN_sensitivity(:))) < tol;
end

if ~naka_rushton_RBN_sensitivity_match
   b = 0;
   errormsg = ['Naka Rushton RBN sensitivity differed by greater than ' num2str(tol) '.'];
   return;
end

% Comparing naka_rushton_RBNS_parameters

tol = %put in appropriate tolorance;

naka_rushton_RBNS_parameters_match = vlt.data.sizeeq(doc_e.naka_rushton_RBNS_parameters(:),doc_a.naka_rushton_RBNS_parameters(:));
if naka_rushton_RBNS_parameters_match
   naka_rushton_RBNS_parameters_match = max(abs(doc_e.naka_rushton_RBNS_parameters(:) - doc_a.naka_rushton_RBNS_parameters(:))) < tol;
end

if ~naka_rushton_RBNS_parameters_match
   b = 0;
   errormsg = ['Naka Rushton RBNS parameters differed by greater than ' num2str(tol) '.'];
   return;
end

% Comparing naka_rushton_RBNS_contrast

tol = %put in appropriate tolorance;

naka_rushton_RBNS_contrast_match = vlt.data.sizeeq(doc_e.naka_rushton_RBNS_contrast(:),doc_a.naka_rushton_RBNS_contrast(:));
if naka_rushton_RBNS_contrast_match
   naka_rushton_RBNS_contrast_match = max(abs(doc_e.naka_rushton_RBNS_contrast(:) - doc_a.naka_rushton_RBNS_contrast(:))) < tol;
end

if ~naka_rushton_RBNS_contrast_match
   b = 0;
   errormsg = ['Naka Rushton RBNS contrast differed by greater than ' num2str(tol) '.'];
   return;
end

% Comparing naka_rushton_RBNS_values

tol = %put in appropriate tolorance;

naka_rushton_RBNS_values_match = vlt.data.sizeeq(doc_e.naka_rushton_RBNS_values(:),doc_a.naka_rushton_RBNS_values(:));
if naka_rushton_RBNS_values_match
   naka_rushton_RBNS_values_match = max(abs(doc_e.naka_rushton_RBNS_values(:) - doc_a.naka_rushton_RBNS_values(:))) < tol;
end

if ~naka_rushton_RBNS_values_match
   b = 0;
   errormsg = ['Naka Rushton RBNS values differed by greater than ' num2str(tol) '.'];
   return;
end

% Comparing naka_rushton_RBNS_pref

tol = %put in appropriate tolorance;

naka_rushton_RBNS_pref_match = vlt.data.sizeeq(doc_e.naka_rushton_RBNS_pref(:),doc_a.naka_rushton_RBNS_pref(:));
if naka_rushton_RBNS_pref_match
   naka_rushton_RBNS_pref_match = max(abs(doc_e.naka_rushton_RBNS_pref(:) - doc_a.naka_rushton_RBNS_pref(:))) < tol;
end

if ~naka_rushton_RBNS_pref_match
   b = 0;
   errormsg = ['Naka Rushton RBNS pref differed by greater than ' num2str(tol) '.'];
   return;
end

% Comparing naka_rushton_RBNS_empirical_c50

tol = %put in appropriate tolorance;

naka_rushton_RBNS_empirical_c50_match = vlt.data.sizeeq(doc_e.naka_rushton_RBNS_empirical_c50(:),doc_a.naka_rushton_RBNS_empirical_c50(:));
if naka_rushton_RBNS_empirical_c50_match
   naka_rushton_RBNS_empirical_c50_match = max(abs(doc_e.naka_rushton_RBNS_empirical_c50(:) - doc_a.naka_rushton_RBNS_empirical_c50(:))) < tol;
end

if ~naka_rushton_RBNS_empirical_c50_match
   b = 0;
   errormsg = ['Naka Rushton RBNS empirical c50 differed by greater than ' num2str(tol) '.'];
   return;
end

% Comparing naka_rushton_RBNS_r2

tol = %put in appropriate tolorance;

naka_rushton_RBNS_r2_match = vlt.data.sizeeq(doc_e.naka_rushton_RBNS_r2(:),doc_a.naka_rushton_RBNS_r2(:));
if naka_rushton_RBNS_r2_match
   naka_rushton_RBNS_r2_match = max(abs(doc_e.naka_rushton_RBNS_r2(:) - doc_a.naka_rushton_RBNS_r2(:))) < tol;
end

if ~naka_rushton_RBNS_r2_match
   b = 0;
   errormsg = ['Naka Rushton RBNS r2 differed by greater than ' num2str(tol) '.'];
   return;
end

% Comparing naka_rushton_RBNS_relative_max_gain

tol = %put in appropriate tolorance;

naka_rushton_RBNS_relative_max_gain_match = vlt.data.sizeeq(doc_e.naka_rushton_RBNS_relative_max_gain(:),doc_a.naka_rushton_RBNS_relative_max_gain(:));
if naka_rushton_RBNS_relative_max_gain_match
   naka_rushton_RBNS_relative_max_gain_match = max(abs(doc_e.naka_rushton_RBNS_relative_max_gain(:) - doc_a.naka_rushton_RBNS_relative_max_gain(:))) < tol;
end

if ~naka_rushton_RBNS_relative_max_gain_match
   b = 0;
   errormsg = ['Naka Rushton RBNS relative max gain differed by greater than ' num2str(tol) '.'];
   return;
end

% Comparing naka_rushton_RBNS_saturation_index

tol = %put in appropriate tolorance;

naka_rushton_RBNS_saturation_index_match = vlt.data.sizeeq(doc_e.naka_rushton_RBNS_saturation_index(:),doc_a.naka_rushton_RBNS_saturation_index(:));
if naka_rushton_RBNS_saturation_index_match
   naka_rushton_RBNS_saturation_index_match = max(abs(doc_e.naka_rushton_RBNS_saturation_index(:) - doc_a.naka_rushton_RBNS_saturation_index(:))) < tol;
end

if ~naka_rushton_RBNS_saturation_index_match
   b = 0;
   errormsg = ['Naka Rushton RBNS relative saturation index differed by greater than ' num2str(tol) '.'];
   return;
end

% Comparing naka_rushton_RBNS_sensitivity

tol = %put in appropriate tolorance;

naka_rushton_RBNS_sensitivity_match = vlt.data.sizeeq(doc_e.naka_rushton_RBNS_sensitivity(:),doc_a.naka_rushton_RBNS_sensitivity(:));
if naka_rushton_RBNS_sensitivity_match
   naka_rushton_RBNS_sensitivity_match = max(abs(doc_e.naka_rushton_RBNS_sensitivity(:) - doc_a.naka_rushton_RBNS_sensitivity(:))) < tol;
end

if ~naka_rushton_RBNS_sensitivity_match
   b = 0;
   errormsg = ['Naka Rushton RBNS sensitivity differed by greater than ' num2str(tol) '.'];
   return;
end

