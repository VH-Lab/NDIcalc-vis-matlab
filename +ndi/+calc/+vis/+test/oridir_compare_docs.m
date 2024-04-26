function [b_, errormsg_] = oridir_compare_docs(document_expected, document_actual, scope)
% ORIDIR_COMPARE_DOCS
%
% [B, ERRORMSG] = ndi.calc.vis.test.oridir_compare_docs(DOC_EXPECTED, DOC_ACTUAL, SCOPE)
%

% Initialize b_ as a row vector of ones for tracking comparison results.
% Initialize errormsg_ as an empty cell array to hold error messages.

b_ = ones(1,35);
errormsg_ = cell(1,35);

 % establish scope-dependent tolerances
switch(scope),
    case 'standard',
    
       tolerance.tuning_curve.direction = 1e-6;
       tolerance.tuning_curve.mean = 1e-6;
       tolerance.tuning_curve.stddev = 1e-6;
       tolerance.tuning_curve.stderr = 1e-6;
       tolerance.tuning_curve.individual = 0.5;
       tolerance.tuning_curve.raw_individual = 0.5;
       tolerance.tuning_curve.control_individual = 0.5;
       tolerance.significance.visual_response_anova_p = 0.1;
       tolerance.significance.across_stimuli_anova_p = 0.1;
       tolerance.vector.circular_variance = 0.1;
       tolerance.vector.direction_circular_variance = 0.1;
       tolerance.vector.Hotelling2Test = 0.1;
       tolerance.vector.orientation_preference = 10;
       tolerance.vector.direction_preference = 0.1;
       tolerance.vector.direction_hotelling2test = 0.1;
       tolerance.vector.dot_direction_significance = 0.1;
       tolerance.fit.double_gaussian_parameters(1) = 5;
       tolerance.fit.double_gaussian_parameters(2) = 5;
       tolerance.fit.double_gaussian_parameters(3) = 10;
       tolerance.fit.double_gaussian_parameters(4) = 20;
       tolerance.fit.double_gaussian_parameters(5) = 10;
       tolerance.fit.double_gaussian_fit_angles = 1e-6;
       tolerance.fit.double_gaussian_fit_values = 5;
       tolerance.fit.orientation_angle_preference = 10;
       tolerance.fit.direction_angle_preference = 10;
       tolerance.fit.hwhh = 10;
       tolerance.fit.orientation_preferred_orthogonal_ratio = 0.1;
       tolerance.fit.direction_preferred_null_ratio = 0.1;
       tolerance.fit.orientation_preferred_orthogonal_ratio_rectified = 0.1;
       tolerance.fit.direction_preferred_null_ratio_rectified = 0.1;


    case 'low_noise',
    
       tolerance.tuning_curve.direction = 1e-6;
       tolerance.tuning_curve.mean = 2;
       tolerance.tuning_curve.stddev = 1e-6;
       tolerance.tuning_curve.stderr = 1e-6;
       tolerance.tuning_curve.individual = 0.5;
       tolerance.tuning_curve.raw_individual = 0.5;
       tolerance.tuning_curve.control_individual = 0.5;
       tolerance.significance.visual_response_anova_p = 0.1;
       tolerance.significance.across_stimuli_anova_p = 0.1;
       tolerance.vector.circular_variance = 0.1;
       tolerance.vector.direction_circular_variance = 0.1;
       tolerance.vector.Hotelling2Test = 0.1;
       tolerance.vector.orientation_preference = 10;
       tolerance.vector.direction_preference = 0.1;
       tolerance.vector.direction_hotelling2test = 0.1;
       tolerance.vector.dot_direction_significance = 0.1;
       tolerance.fit.double_gaussian_parameters(1) = 5;
       tolerance.fit.double_gaussian_parameters(2) = 5;
       tolerance.fit.double_gaussian_parameters(3) = 10;
       tolerance.fit.double_gaussian_parameters(4) = 20;
       tolerance.fit.double_gaussian_parameters(5) = 10;
       tolerance.fit.double_gaussian_fit_angles = 1e-6;
       tolerance.fit.double_gaussian_fit_values = 5;
       tolerance.fit.orientation_angle_preference = 10;
       tolerance.fit.direction_angle_preference = 10;
       tolerance.fit.hwhh = 10;
       tolerance.fit.orientation_preferred_orthogonal_ratio = 0.1;
       tolerance.fit.direction_preferred_null_ratio = 0.1;
       tolerance.fit.orientation_preferred_orthogonal_ratio_rectified = 0.1;
       tolerance.fit.direction_preferred_null_ratio_rectified = 0.1;


    case 'high_noise',
    
        tolerance.tuning_curve.direction = 1e-6;
        tolerance.tuning_curve.mean = 5;
        tolerance.tuning_curve.stddev = 1e-6;
        tolerance.tuning_curve.stddev = 1e-6;
        tolerance.tuning_curve.stderr = 1e-6;
        tolerance.tuning_curve.individual = 0.5;
        tolerance.tuning_curve.raw_individual = 0.5;
        tolerance.tuning_curve.control_individual = 0.5;
        tolerance.significance.visual_response_anova_p = 0.1;
        tolerance.significance.across_stimuli_anova_p = 0.1;
        tolerance.vector.circular_variance = 0.1;
        tolerance.vector.direction_circular_variance = 0.1;
        tolerance.vector.Hotelling2Test = 0.1;
        tolerance.vector.orientation_preference = 10;
        tolerance.vector.direction_preference = 0.1;
        tolerance.vector.direction_hotelling2test = 0.1;
        tolerance.vector.dot_direction_significance = 0.1;
        tolerance.fit.double_gaussian_parameters(1) = 5;
        tolerance.fit.double_gaussian_parameters(2) = 5;
        tolerance.fit.double_gaussian_parameters(3) = 10;
        tolerance.fit.double_gaussian_parameters(4) = 20;
        tolerance.fit.double_gaussian_parameters(5) = 10;
        tolerance.fit.double_gaussian_fit_angles = 1e-6;
        tolerance.fit.double_gaussian_fit_values = 5;
        tolerance.fit.orientation_angle_preference = 10;
        tolerance.fit.direction_angle_preference = 10;
        tolerance.fit.hwhh = 10;
        tolerance.fit.orientation_preferred_orthogonal_ratio = 0.1;
        tolerance.fit.direction_preferred_null_ratio = 0.1;
        tolerance.fit.orientation_preferred_orthogonal_ratio_rectified = 0.1;
        tolerance.fit.direction_preferred_null_ratio_rectified = 0.1;


    otherwise,
       error(['Unknown scope ' scope '.']);
end;

% start comparison

doc_e = document_expected.document_properties.orientation_direction_tuning;
doc_a = document_actual.document_properties.orientation_direction_tuning;

% Comparing Properties
%   coordinates
%   response_units
%   response_type

properties_cmatch = strcmpi(char(doc_e.properties.coordinates), char(doc_a.properties.coordinates));
if ~properties_cmatch
    b_(35) = 0;
    errormsg_{35} = ['Expected coordinates of ' doc_e.properties.coordinates ' but observed ' doc_a.properties.coordinates];
end

properties_rmatch = strcmpi(char(doc_e.properties.response_units), char(doc_a.properties.response_units));
if ~properties_rmatch
    b_(34) = 0;
    errormsg_{34} = ['Expected response units of ' doc_e.properties.response_units ' but observed ' doc_a.properties.response_units];
end

properties_match = strcmpi(char(doc_e.properties.response_type), char(doc_a.properties.response_type));
if ~properties_match
   b_(33) = 0;
   errormsg_{33} = ['Expected response type of ' doc_e.properties.response_type ' but observed ' doc_a.properties.response_type];
end

% Comparing tuning_curve
%   direction
%   mean
%   stddev
%   stderr
%   individual
%   raw_individual
%   control_individual
%b_=[]
%errormsg_={}
[b_(1),errormsg_{1}] = ndi.test.values_within_tolerance(doc_e.tuning_curve.direction, doc_a.tuning_curve.direction, tolerance.tuning_curve.direction, 'direction');
[b_(2),errormsg_{2}] = ndi.test.values_within_tolerance(doc_e.tuning_curve.mean, doc_a.tuning_curve.mean, tolerance.tuning_curve.mean, 'mean');
[b_(3),errormsg_{3}] = ndi.test.values_within_tolerance(doc_e.tuning_curve.stddev, doc_a.tuning_curve.stddev, tolerance.tuning_curve.stddev, 'stddev');
[b_(4),errormsg_{4}] = ndi.test.values_within_tolerance(doc_e.tuning_curve.stderr, doc_a.tuning_curve.stderr, tolerance.tuning_curve.stderr, 'stderr');
[b_(5),errormsg_{5}] = ndi.test.values_within_tolerance(doc_e.tuning_curve.individual, doc_a.tuning_curve.individual, tolerance.tuning_curve.individual, 'individual');
[b_(6),errormsg_{6}] = ndi.test.values_within_tolerance(doc_e.tuning_curve.raw_individual, doc_a.tuning_curve.raw_individual, tolerance.tuning_curve.raw_individual, 'raw individual');
[b_(7),errormsg_{7}] = ndi.test.values_within_tolerance(doc_e.tuning_curve.control_individual, doc_a.tuning_curve.control_individual, tolerance.tuning_curve.control_individual, 'control individual');

% Comparing significance
%   visual_response_anova_p
%   across_stimuli_anova_p

[b_(8),errormsg_{8}] = ndi.test.values_within_tolerance(doc_e.significance.visual_response_anova_p, doc_a.significance.visual_response_anova_p, tolerance.significance.visual_response_anova_p, 'visual response anova p');
[b_(9),errormsg_{9}] = ndi.test.values_within_tolerance(doc_e.significance.across_stimuli_anova_p, doc_a.significance.across_stimuli_anova_p, tolerance.significance.across_stimuli_anova_p, 'across stimli anova p');

% Comparing vector
%   circular_variance
%   direction_circular_variance
%   hotelling2test
%   orientation_preference
%   direction_preference
%   direction_hotelling2test
%   dot_direction_significance

[b_(10),errormsg_{10}] = ndi.test.values_within_tolerance(doc_e.vector.circular_variance, doc_a.vector.circular_variance, tolerance.vector.circular_variance, 'circular variance');
[b_(11),errormsg_{11}] = ndi.test.values_within_tolerance(doc_e.vector.direction_circular_variance, doc_a.vector.direction_circular_variance, tolerance.vector.direction_circular_variance, 'direction circular variance');
[b_(12),errormsg_{12}] = ndi.test.values_within_tolerance(doc_e.vector.Hotelling2Test, doc_a.vector.Hotelling2Test, tolerance.vector.Hotelling2Test, 'hotelling2test');
[b_(13),errormsg_{13}] = ndi.test.values_within_tolerance(doc_e.vector.orientation_preference, doc_a.vector.orientation_preference, tolerance.vector.orientation_preference, 'orientation preference');
[b_(14),errormsg_{14}] = ndi.test.values_within_tolerance(doc_e.vector.direction_preference, doc_a.vector.direction_preference, tolerance.vector.direction_preference, 'direction preference');
[b_(15),errormsg_{15}] = ndi.test.values_within_tolerance(doc_e.vector.direction_hotelling2test, doc_a.vector.direction_hotelling2test, tolerance.vector.direction_hotelling2test, 'direction hotelling2test');
[b_(16),errormsg_{16}] = ndi.test.values_within_tolerance(doc_e.vector.dot_direction_significance, doc_a.vector.dot_direction_significance, tolerance.vector.dot_direction_significance, 'dot direction significance');

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

[b_(17),errormsg_{17}] = ndi.test.values_within_tolerance(doc_e.fit.double_gaussian_parameters(1), doc_a.fit.double_gaussian_parameters(1), tolerance.fit.double_gaussian_parameters(1), 'double gaussian parameters #1');
[b_(18),errormsg_{18}] = ndi.test.values_within_tolerance(doc_e.fit.double_gaussian_parameters(2), doc_a.fit.double_gaussian_parameters(2), tolerance.fit.double_gaussian_parameters(2), 'double gaussian parameters #2');
[b_(19),errormsg_{19}] = ndi.test.values_within_tolerance(doc_e.fit.double_gaussian_parameters(3), doc_a.fit.double_gaussian_parameters(3), tolerance.fit.double_gaussian_parameters(3), 'double gaussian parameters #3');
[b_(20),errormsg_{20}] = ndi.test.values_within_tolerance(doc_e.fit.double_gaussian_parameters(4), doc_a.fit.double_gaussian_parameters(4), tolerance.fit.double_gaussian_parameters(4), 'double gaussian parameters #4');
[b_(21),errormsg_{21}] = ndi.test.values_within_tolerance(doc_e.fit.double_gaussian_parameters(5), doc_a.fit.double_gaussian_parameters(5), tolerance.fit.double_gaussian_parameters(5), 'double gaussian parameters #1');
[b_(22),errormsg_{22}] = ndi.test.values_within_tolerance(doc_e.fit.double_gaussian_fit_angles, doc_a.fit.double_gaussian_fit_angles, tolerance.fit.double_gaussian_fit_angles, 'double gaussian fit angles');
[b_(23),errormsg_{23}] = ndi.test.values_within_tolerance(doc_e.fit.double_gaussian_fit_values, doc_a.fit.double_gaussian_fit_values, tolerance.fit.double_gaussian_fit_values, 'double gaussian fit values');
[b_(24),errormsg_{24}] = ndi.test.values_within_tolerance(doc_e.fit.orientation_angle_preference, doc_a.fit.orientation_angle_preference, tolerance.fit.orientation_angle_preference, 'orientation angle preference');
[b_(25),errormsg_{25}] = ndi.test.values_within_tolerance(doc_e.fit.direction_angle_preference, doc_a.fit.direction_angle_preference, tolerance.fit.direction_angle_preference, 'direction angle preference');
[b_(26),errormsg_{26}] = ndi.test.values_within_tolerance(doc_e.fit.hwhh, doc_a.fit.hwhh, tolerance.fit.hwhh, 'hwhh');
[b_(27),errormsg_{27}] = ndi.test.values_within_tolerance(doc_e.fit.orientation_preferred_orthogonal_ratio, doc_a.fit.orientation_preferred_orthogonal_ratio, tolerance.fit.orientation_preferred_orthogonal_ratio, 'orientation preferred orthogonal ratio');
[b_(28),errormsg_{28}] = ndi.test.values_within_tolerance(doc_e.fit.direction_preferred_null_ratio, doc_a.fit.direction_preferred_null_ratio, tolerance.fit.direction_preferred_null_ratio, 'direction preferred null ratio');
[b_(29),errormsg_{29}] = ndi.test.values_within_tolerance(doc_e.fit.orientation_preferred_orthogonal_ratio_rectified, doc_a.fit.orientation_preferred_orthogonal_ratio_rectified, tolerance.fit.orientation_preferred_orthogonal_ratio_rectified, 'orientation preferred orthogonal ratio rectified');
[b_(30),errormsg_{30}] = ndi.test.values_within_tolerance(doc_e.fit.direction_preferred_null_ratio_rectified, doc_a.fit.direction_preferred_null_ratio_rectified, tolerance.fit.direction_preferred_null_ratio_rectified, 'direction preferred null ratio rectified');

% The following code does three things:
% 1. Identify the b_ values with unmatched results
% 2. Update b_ to only include those
% 3. Update the corresponding errormsg_ messages

if any(b_==0),
    error_indices = find(b_==0);
    b_ = b_(error_indices);
    errormsg_ = errormsg_(error_indices);
end


end
