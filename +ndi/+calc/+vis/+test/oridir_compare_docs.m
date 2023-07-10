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
%   coordinates
%   response_units
%   response_type

properties_match = strcmpi(char(doc_e.properties.coordinates), char(doc_a.properties.coordinates));
if ~properties_match
   b = 0;
   errormsg = ['Expected coordinates of ' doc_e.properties.coordinates ' but observed ' doc_a.properties.coordinates];
   return;
end

properties_match = strcmpi(char(doc_e.properties.response_units), char(doc_a.properties.response_units));
if ~properties_match
   b = 0;
   errormsg = ['Expected response units of ' doc_e.properties.response_units ' but observed ' doc_a.properties.response_units];
   return;
end

properties_match = strcmpi(char(doc_e.properties.response_type), char(doc_a.properties.response_type));
if ~properties_match
   b = 0;
   errormsg = ['Expected response type of ' doc_e.properties.response_type ' but observed ' doc_a.properties.response_type];
   return;
end

% Comparing tuning_curve
%   direction
%   mean
%   stddev
%   stderr
%   individual
%   raw_individual
%   control_individual

[b,errormsg] = ndi.test.values_within_tolerance(doc_e.tuning_curve.direction, doc_a.tuning_curve.direction, tolerance, 'direction');
[b,errormsg] = ndi.test.values_within_tolerance(doc_e.tuning_curve.mean, doc_a.tuning_curve.mean, tolerance, 'mean');
[b,errormsg] = ndi.test.values_within_tolerance(doc_e.tuning_curve.stddev, doc_a.tuning_curve.stddev, tolerance, 'stddev');
[b,errormsg] = ndi.test.values_within_tolerance(doc_e.tuning_curve.stderr, doc_a.tuning_curve.stderr, tolerance, 'stderr');
[b,errormsg] = ndi.test.values_within_tolerance(doc_e.tuning_curve.individual, doc_a.tuning_curve.individual, tolerance, 'individual');
[b,errormsg] = ndi.test.values_within_tolerance(doc_e.tuning_curve.raw_individual, doc_a.tuning_curve.raw_individual, tolerance, 'raw individual');
[b,errormsg] = ndi.test.values_within_tolerance(doc_e.tuning_curve.control_individual, doc_a.tuning_curve.control_individual, tolerance, 'control individual');

% Comparing significance
%   visual_response_anova_p
%   across_stimuli_anova_p

[b,errormsg] = ndi.test.values_within_tolerance(doc_e.significance.visual_response_anova_p, doc_a.significance.visual_response_anova_p, tolerance, 'visual response anova p');
[b,errormsg] = ndi.test.values_within_tolerance(doc_e.significance.across_stimuli_anova_p, doc_a.significance.across_stimuli_anova_p, tolerance, 'across stimli anova p');

% Comparing vector
%   circular_variance
%   direction_circular_variance
%   hotelling2test
%   orientation_preference
%   direction_preference
%   direction_hotelling2test
%   dot_direction_significance

[b,errormsg] = ndi.test.values_within_tolerance(doc_e.vector.circular_variance, doc_a.vector.circular_variance, tolerance, 'circular variance');
[b,errormsg] = ndi.test.values_within_tolerance(doc_e.vector.direction_circular_variance, doc_a.vector.direction_circular_variance, tolerance, 'direction circular variance');
[b,errormsg] = ndi.test.values_within_tolerance(doc_e.vector.hotelling2test, doc_a.vector.hotelling2test, tolerance, 'hotelling2test');
[b,errormsg] = ndi.test.values_within_tolerance(doc_e.vector.orientation_preference, doc_a.vector.orientation_preference, tolerance, 'orientation preference');
[b,errormsg] = ndi.test.values_within_tolerance(doc_e.vector.direction_preference, doc_a.vector.direction_preference, tolerance, 'direction preference');
[b,errormsg] = ndi.test.values_within_tolerance(doc_e.vector.direction_hotelling2test, doc_a.vector.direction_hotelling2test, tolerance, 'direction hotelling2test');
[b,errormsg] = ndi.test.values_within_tolerance(doc_e.vector.dot_direction_significance, doc_a.vector.dot_direction_significance, tolerance, 'dot direction significance');

% Comparing fit
%   double_gaussian_parameters
%   double_gaussian_fit_angles
%   double_gaussian_fit_values		
%   orientation_angle_preference
%   direction_angle_preference
%   hwhh
%   orientation_preferred_orthogonal_ratio
%   direction_preferred_null_ratio
%   orientation_preferred_orthogonal_ratio_rectified
%   direction_preferred_null_ratio_rectified

[b,errormsg] = ndi.test.values_within_tolerance(doc_e.fit.double_gaussian_parameters, doc_a.fit.double_gaussian_parameters, tolerance, 'double gaussian parameters');
[b,errormsg] = ndi.test.values_within_tolerance(doc_e.fit.double_gaussian_fit_angles, doc_a.fit.double_gaussian_fit_angles, tolerance, 'double gaussian fit angles');
[b,errormsg] = ndi.test.values_within_tolerance(doc_e.fit.double_gaussian_fit_values, doc_a.fit.double_gaussian_fit_values, tolerance, 'double gaussian fit values');
[b,errormsg] = ndi.test.values_within_tolerance(doc_e.fit.orientation_angle_preference, doc_a.fit.orientation_angle_preference, tolerance, 'orientation angle preference');
[b,errormsg] = ndi.test.values_within_tolerance(doc_e.fit.direction_angle_preference, doc_a.fit.direction_angle_preference, tolerance, 'direction angle preference');
[b,errormsg] = ndi.test.values_within_tolerance(doc_e.fit.hwhh, doc_a.fit.hwhh, tolerance, 'hwhh');
[b,errormsg] = ndi.test.values_within_tolerance(doc_e.fit.orientation_preferred_orthogonal_ratio, doc_a.fit.orientation_preferred_orthogonal_ratio, tolerance, 'orientation preferred orthogonal ratio');
[b,errormsg] = ndi.test.values_within_tolerance(doc_e.fit.direction_preferred_null_ratio, doc_a.fit.direction_preferred_null_ratio, tolerance, 'direction preferred null ratio');
[b,errormsg] = ndi.test.values_within_tolerance(doc_e.fit.orientation_preferred_orthogonal_ratio_rectified, doc_a.fit.orientation_preferred_orthogonal_ratio_rectified, tolerance, 'orientation preferred orthogonal ratio rectified');
[b,errormsg] = ndi.test.values_within_tolerance(doc_e.fit.direction_preferred_null_ratio_rectified, doc_a.fit.direction_preferred_null_ratio_rectified, tolerance, 'direction preferred null ratio rectified');

end
