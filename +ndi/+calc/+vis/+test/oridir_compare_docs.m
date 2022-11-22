function [b, errormsg] = oridir_compare_docs(document_expected, document_actual, scope)
% ORIDIR_COMPARE_DOCS
%
% [B, ERRORMSG] = ndi.calc.vis.test.oridir_compare_docs(DOC_EXPECTED, DOC_ACTUAL, SCOPE)
%


b = 1;
errormsg = '';

% start comparison

doc_e = document_expected.document_properties.orientation_direction_tuning;
doc_a = document_actual.document_properties.orientation_direction_tuning;

% Comparing Properties
% Comparing Properties.coordinates

properties_match = strcmpi(doc_e.properties.coordinates, doc_a.properties.coordinates);
if ~properties_match
   b = 0;
   errormsg = ['Expected properties.coordinates field of ' doc_e.properties.coordinates ' but observed ' doc_a.properties.coordinates];
   return;
end

% Comparing Tuning_curve
% Tuning_curve.direction match  within 1e-6

diretol = 1e-6;

direction_match = vlt.data.sizeeq(doc_e.tuning_curve.direction(:),doc_a.tuning_curve.direction(:));
if direction_match
   direction_match = max(abs(doc_e.tuning_curve.direction(:) - doc_a.tuning_curve.direction(:))) < diretol;
end

if ~direction_match
   b = 0;
   errormsg = ['Direction angles differed by greater than ' num2str(diretol) '.'];
end

% Tuning_curve.mean match within 3

meantol = 3;

mean_match = vlt.data.sizeeq(doc_e.tuning_curve.mean(:),doc_a.tuning_curve.mean(:));
if mean_match 
   mean_match = max(abs(doc_e.tuning_curve.mean(:) - doc_a.tuning_curve.mean(:))) < meantol;
end

if ~mean_match
   b = 0;
   errormsg = ['Mean angles differed by greater than ' num2str(meantol) '.'];
end

% Tuning_curve.stddev match exactly

stddev_match = vlt.data.sizeeq(doc_e.tuning_curve.stddev(:),doc_a.tuning_curve.stddev(:));
if stddev_match
   stddev_match = 1;
end

if ~stddev_match
   b = 0;
   errormsg = 'Stddev are different.';
end

% Tuning_curve.stderr match exactly

stderr_match = vlt.data.sizeeq(doc_e.tuning_curve.stderr(:),doc_a.tuning_curve.stderr(:));
if stderr_match
   stderr_match = 1;
end

if ~stderr_match
   b = 0;
   errormsg = 'Stderr are different.';
end

% Comparing Significance
% Significance.visual_response_anova_p match within 0.1

vratol = 0.1;

vra_match = vlt.data.sizeeq(doc_e.significance.visual_response_anova_p,doc_a.significance.visual_response_anova_p);
if vra_match 
   vra_match = max(abs(doc_e.significance.visual_response_anova_p - doc_a.significance.visual_response_anova_p)) < vratol;
end

if ~vra_match
   b = 0;
   errormsg = ['Visual Response Anova P differed by greater than ' num2str(vratol) '.'];
end

% Significance.across_stimuli_anova_p match within 0.1

asatol = 0.1;

asa_match = vlt.data.sizeeq(doc_e.significance.across_stimuli_anova_p,doc_a.significance.across_stimuli_anova_p);
if asa_match 
   asa_match = max(abs(doc_e.significance.across_stimuli_anova_p - doc_a.significance.across_stimuli_anova_p)) < asatol;
end

if ~asa_match
   b = 0;
   errormsg = ['Across Stimuli Anova P differed by greater than ' num2str(asatol) '.'];
end

% Comparing vector
% vector.direction_hotelling2test match within 0.1

dhttol = 0.1;

dht_match = vlt.data.sizeeq(doc_e.vector.direction_hotelling2test,doc_a.vector.direction_hotelling2test);
if dht_match 
   dht_match = max(abs(doc_e.vector.direction_hotelling2test - doc_a.vector.direction_hotelling2test)) < dhttol;
end

if ~dht_match
   b = 0;
   errormsg = ['Direction Hotelling 2 Test differed by greater than ' num2str(dhttol) '.'];
end

% vector.Hotelling2Test match within 0.1

httol = 0.1;

ht_match = vlt.data.sizeeq(doc_e.vector.Hotelling2Test,doc_a.vector.Hotelling2Test);
if ht_match 
   ht_match = max(abs(doc_e.vector.Hotelling2Test - doc_a.vector.Hotelling2Test)) < httol;
end

if ~ht_match
   b = 0;
   errormsg = ['Hotelling 2 Test differed by greater than ' num2str(httol) '.'];
end

% vector.dot_direction_significance match within 0.1

ddstol = 0.1;

dds_match = vlt.data.sizeeq(doc_e.vector.dot_direction_significance,doc_a.vector.dot_direction_significance);
if dds_match 
   dds_match = max(abs(doc_e.vector.dot_direction_significance - doc_a.vector.dot_direction_significance)) < ddstol;
end

if ~dds_match
   b = 0;
   errormsg = ['Dot direction significance differed by greater than ' num2str(ddstol) '.'];
end

% vector.circular_variance match within 0.1

cvtol = 0.1;

cv_match = vlt.data.sizeeq(doc_e.vector.circular_variance,doc_a.vector.circular_variance);
if cv_match 
   cv_match = max(abs(doc_e.vector.circular_variance - doc_a.vector.circular_variance)) < cvtol;
end

if ~cv_match
   b = 0;
   errormsg = ['Circular variance differed by greater than ' num2str(cvtol) '.'];
end

% vector.orientation_preference match within 10 degrees

optol = 10;

op_match = vlt.data.sizeeq(doc_e.vector.orientation_preference,doc_a.vector.orientation_preference);
if op_match 
   op_match = max(abs(doc_e.vector.orientation_preference - doc_a.vector.orientation_preference)) < optol;
end

if ~op_match
   b = 0;
   errormsg = ['Orientation preference differed by greater than ' num2str(optol) '.'];
end

% vector.direction_preference match within 10 degrees

 % check to make sure not NaN

if 0,

dptol = 10;

dp_match = vlt.data.sizeeq(doc_e.vector.direction_preference,doc_a.vector.direction_preference);
if dp_match 
   dp_match = max(abs(doc_e.vector.direction_preference - doc_a.vector.direction_preference)) < dptol;
end

if ~dp_match
	b = 0;
	errormsg = ['Direction preference differed by greater than ' num2str(dptol) '.'];
	return;
end

end;

% vector.direction_circular_variance match within 0.1

dcvtol = 0.1;

if 0,
dcv_match = vlt.data.sizeeq(doc_e.vector.direction_circular_variance,doc_a.vector.direction_circular_variance);
if dcv_match 
   dcv_match = max(abs(doc_e.vector.direction_circular_variance - doc_a.vector.direction_circular_variance)) < dcvtol;
end

if ~dcv_match
   b = 0;
   errormsg = ['Direction circular variance differed by greater than ' num2str(dcvtol) '.'];
end

end;

% Comparing fit
% fit.orientation_angle_preference match within 10 degrees

oaptol = 10;

oap_match = vlt.data.sizeeq(doc_e.fit.orientation_angle_preference,doc_a.fit.orientation_angle_preference);
if oap_match 
   oap_match = max(abs(doc_e.fit.orientation_angle_preference - doc_a.fit.orientation_angle_preference)) < oaptol;
end

if ~oap_match
   b = 0;
   errormsg = ['Orientation angle preference differed by greater than ' num2str(oaptol) '.'];
end

% fit.direction_angle_preference match within 10 degrees

daptol = 10;

dap_match = vlt.data.sizeeq(doc_e.fit.direction_angle_preference,doc_a.fit.direction_angle_preference);
if dap_match 
   dap_match = max(abs(doc_e.fit.direction_angle_preference - doc_a.fit.direction_angle_preference)) < daptol;
end

if ~dap_match
   b = 0;
   errormsg = ['Direction angle preference differed by greater than ' num2str(daptol) '.'];
end

% fit.hwhh match within 10 degrees

hwhhtol = 10;

hwhh_match = vlt.data.sizeeq(doc_e.fit.hwhh,doc_a.fit.hwhh);
if hwhh_match 
   hwhh_match = max(abs(doc_e.fit.hwhh - doc_a.fit.hwhh)) < hwhhtol;
end

if ~hwhh_match
   b = 0;
   errormsg = ['Half width half height differed by greater than ' num2str(hwhhtol) '.'];
end

% fit.direction_preferred_null_ratio match within 0.1

dpnrtol = 0.1;

dpnr_match = vlt.data.sizeeq(doc_e.fit.direction_preferred_null_ratio,doc_a.fit.direction_preferred_null_ratio);
if dpnr_match 
   dpnr_match = max(abs(doc_e.fit.direction_preferred_null_ratio - doc_a.fit.direction_preferred_null_ratio)) < dpnrtol;
end

if ~dpnr_match
   b = 0;
   errormsg = ['Direction preferred null ratio differed by greater than ' num2str(dpnrtol) '.'];
end

% fit.orientation_preferred_orthogonal_ratio match within 0.1

oportol = 0.1;

opor_match = vlt.data.sizeeq(doc_e.fit.orientation_preferred_orthogonal_ratio,doc_a.fit.orientation_preferred_orthogonal_ratio);
if opor_match 
   opor_match = max(abs(doc_e.fit.orientation_preferred_orthogonal_ratio - doc_a.fit.orientation_preferred_orthogonal_ratio)) < oportol;
end

if ~opor_match
   b = 0;
   errormsg = ['Orientation preferred orthogonal ratio differed by greater than ' num2str(oportol) '.'];
end

% fit.orientation_preferred_orthogonal_ratio_rectified match within 0.1

oporrtol = 0.1;

oporr_match = vlt.data.sizeeq(doc_e.fit.orientation_preferred_orthogonal_ratio_rectified,doc_a.fit.orientation_preferred_orthogonal_ratio_rectified);
if oporr_match 
   oporr_match = max(abs(doc_e.fit.orientation_preferred_orthogonal_ratio_rectified - doc_a.fit.orientation_preferred_orthogonal_ratio_rectified)) < oporrtol;
end

if ~oporr_match
   b = 0;
   errormsg = ['Orientation preferred orthogonal ratio rectified differed by greater than ' num2str(oporrtol) '.'];
end

% fit.directional_preferred_null_ratio_rectified match within 0.1

dpnrrtol = 0.1;

dpnrr_match = vlt.data.sizeeq(doc_e.fit.direction_preferred_null_ratio_rectified,doc_a.fit.direction_preferred_null_ratio_rectified);
if dpnrr_match 
   dpnrr_match = max(abs(doc_e.fit.direction_preferred_null_ratio_rectified - doc_a.fit.direction_preferred_null_ratio_rectified)) < dpnrrtol;
end

if ~dpnrr_match
   b = 0;
   errormsg = ['Directional preferred null ratio rectified differed by greater than ' num2str(dpnrrtol) '.'];
end

% Comparing fit.double_guassian_parameters
% First two parameters match within 5

dgpon = 5;

dbpo_match = vlt.data.sizeeq(doc_e.fit.double_guassian_parameters(1),doc_a.fit.double_guassian_parameters(1));
if dbpo_match 
   dbpo_match = max(abs(doc_e.fit.double_guassian_parameters(1) - doc_a.fit.double_guassian_parameters(1))) < dgpon;
end

if ~dbpo_match
   b = 0;
   errormsg = ['Parameter #1 of the double guassian parameters differed by greater than ' num2str(dgpon) '.'];
end

dbpt_match = vlt.data.sizeeq(doc_e.fit.double_guassian_parameters(2),doc_a.fit.double_guassian_parameters(2));
if dbpt_match 
   dbpt_match = max(abs(doc_e.fit.double_guassian_parameters(2) - doc_a.fit.double_guassian_parameters(2))) < dgpon;
end

if ~dbpt_match
   b = 0;
   errormsg = ['Parameter #2 of the double guassian parameters differed by greater than ' num2str(dgpon) '.'];
end

% Second two parameters match within 10

dbptf = 10;

dbpt_match = vlt.data.sizeeq(doc_e.fit.double_guassian_parameters(3),doc_a.fit.double_guassian_parameters(3));
if dbpt_match 
   dbpt_match = max(abs(doc_e.fit.double_guassian_parameters(3) - doc_a.fit.double_guassian_parameters(3))) < dbptf;
end

if ~dbpt_match
   b = 0;
   errormsg = ['Parameter #3 of the double guassian parameters differed by greater than ' num2str(dbptf) '.'];
end

dgpf_pref = 20;

dbpf_match = vlt.data.sizeeq(doc_e.fit.double_guassian_parameters(4),doc_a.fit.double_guassian_parameters(4));
if dbpf_match 
   dbpf_match = max(abs(doc_e.fit.double_guassian_parameters(4) - doc_a.fit.double_guassian_parameters(4))) < dgpf_pref;
end

if ~dbpf_match
   b = 0;
   errormsg = ['Preferred direction parameter of the double guassian fit differed by greater than ' num2str(dgpf_pref) '.'];
end

% Last parameter match within 10

dbpfi = 10;

dbpfi_match = vlt.data.sizeeq(doc_e.fit.double_guassian_parameters(5),doc_a.fit.double_guassian_parameters(5));
if dbpfi_match 
   dbpfi_match = max(abs(doc_e.fit.double_guassian_parameters(5) - doc_a.fit.double_guassian_parameters(5))) < dbpfi;
end

if ~dbpfi_match
   b = 0;
   errormsg = ['Tuning width parameter (#5) of the double guassian parameters differed by greater than ' num2str(dbpfi) '.'];
end

% fit.double_gaussian_fit_angles match exactly

dgfa_match = vlt.data.sizeeq(doc_e.fit.double_gaussian_fit_angles(:),doc_a.fit.double_gaussian_fit_angles(:));
if dgfa_match
   dgfa_match = 1;
end

if ~dgfa_match
   b = 0;
   errormsg = 'Double gaussian fit angles are different.';
end

% fit.double_gaussian_fit_values match within 5

if any(abs(doc_e.fit.double_gaussian_fit_values(:) - doc_a.fit.double_gaussian_fit_values(:)) > 5)
   b = 0;
   errormsg = 'Double gaussian fit values differed by greater than 5.';
end
