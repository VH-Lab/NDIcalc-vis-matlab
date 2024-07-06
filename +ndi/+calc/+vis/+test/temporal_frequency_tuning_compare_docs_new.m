function [b_, errormsg_] = temporal_frequency_tuning_compare_docs(document_expected, document_actual, scope)
% Comparing temporal_frequency_tuning
%
% [B, ERRORMSG] = ndi.calc.vis.test.temporal_frequency_tuning_compare_docs(DOC_EXPECTED, DOC_ACTUAL, SCOPE)
%
%

% Initialize b_ as a row vector of ones for tracking comparison results
% Initialize errormsg_ as an empty cell array to hold error messages

b_ = ones(1,56);
errormsg_ = cell(1,56);

% establish scope-dependent tolerances
switch(scope),
    case 'standard',
    
      tolerance.tuning_curve.temporal_frequency = 0.5;    
      tolerance.tuning_curve.mean = 0.5;
      tolerance.tuning_curve.stddev = 0.5;
      tolerance.tuning_curve.stderr = 0.5;
      tolerance.tuning_curve.individual = 0.5;
      tolerance.tuning_curve.control_mean = 0.5;
      tolerance.tuning_curve.control_stddev = 0.5;
      tolerance.tuning_curve.control_stderr = 0.5;
      tolerance.tuning_curve.control_mean_stddev = 0.5;
      tolerance.tuning_curve.control_mean_stderr = 0.5;
      tolerance.significance.visual_response_anova_p = 0.5;
      tolerance.significance.across_stimuli_anova_p = 0.5;
      tolerance.fitless.L50 = 0.5;
      tolerance.fitless.Pref = 0.5;
      tolerance.fitless.H50 = 0.5;
      tolerance.fitless.bandwidth = 0.5;
      tolerance.fitless.low_pass_index = 0.5;
      tolerance.fitless.high_pass_index = 0.5;
      tolerance.fit_dog.parameters, = 0.5;
      tolerance.fit_dog.values = 0.5;
      tolerance.fit_dog.fit = 0.5;
      tolerance.fit_dog.L50 = 0.5;
      tolerance.fit_dog.Pref = 0.5;
      tolerance.fit_dog.H50 = 0.5;
      tolerance.fit_dog.bandwidth = 0.5;
      tolerance.fit_movshon.parameters = 0.5;
      tolerance.fit_movshon.values = 0.5;
      tolerance.fit_movshon.fit = 0.5;
      tolerance.fit_movshon.L50 = 0.5;
      tolerance.fit_movshon.Pref = 0.5;
      tolerance.fit_movshon.H50 = 0.5;
      tolerance.fit_movshon.R2 = 0.5;
      tolerance.fit_movshon.bandwidth = 0.5;
      tolerance.fit_movshon_c.parameters = 0.5;
      tolerance.fit_movshon_c.values = 0.5;
      tolerance.fit_movshon_c.fit = 0.5;
      tolerance.fit_movshon_c.L50 = 0.5;
      tolerance.fit_movshon_c.Pref = 0.5;
      tolerance.fit_movshon_c.H50 = 0.5;
      tolerance.fit_movshon_c.R2 = 0.5;
      tolerance.fit_movshon_c.bandwidth = 0.5;
      tolerance.fit_spline.values = 0.5;
      tolerance.fit_spline.fit = 0.5;
      tolerance.fit_spline.L50 = 0.5;
      tolerance.fit_spline.Pref = 0.5;
      tolerance.fit_spline.H50 = 0.5;
      tolerance.fit_spline.bandwidth = 0.5;
      tolerance.fit_gausslog.parameters = 0.5;
      tolerance.fit_gausslog.values = 0.5;
      tolerance.fit_gausslog.fit = 0.5;
      tolerance.fit_gausslog.L50 = 0.5;
      tolerance.fit_gausslog.Pref = 0.5;
      tolerance.fit_gausslog.H50 = 0.5;
      tolerance.fit_gausslog.bandwidth = 0.5;

    case 'low_noise',

      tolerance.tuning_curve.temporal_frequency = 0.5;    
      tolerance.tuning_curve.mean = 0.5;
      tolerance.tuning_curve.stddev = 0.5;
      tolerance.tuning_curve.stderr = 0.5;
      tolerance.tuning_curve.individual = 0.5;
      tolerance.tuning_curve.control_mean = 0.5;
      tolerance.tuning_curve.control_stddev = 0.5;
      tolerance.tuning_curve.control_stderr = 0.5;
      tolerance.tuning_curve.control_mean_stddev = 0.5;
      tolerance.tuning_curve.control_mean_stderr = 0.5;
      tolerance.significance.visual_response_anova_p = 0.5;
      tolerance.significance.across_stimuli_anova_p = 0.5;
      tolerance.fitless.L50 = 0.5;
      tolerance.fitless.Pref = 0.5;
      tolerance.fitless.H50 = 0.5;
      tolerance.fitless.bandwidth = 0.5;
      tolerance.fitless.low_pass_index = 0.5;
      tolerance.fitless.high_pass_index = 0.5;
      tolerance.fit_dog.parameters, = 0.5;
      tolerance.fit_dog.values = 0.5;
      tolerance.fit_dog.fit = 0.5;
      tolerance.fit_dog.L50 = 0.5;
      tolerance.fit_dog.Pref = 0.5;
      tolerance.fit_dog.H50 = 0.5;
      tolerance.fit_dog.bandwidth = 0.5;
      tolerance.fit_movshon.parameters = 0.5;
      tolerance.fit_movshon.values = 0.5;
      tolerance.fit_movshon.fit = 0.5;
      tolerance.fit_movshon.L50 = 0.5;
      tolerance.fit_movshon.Pref = 0.5;
      tolerance.fit_movshon.H50 = 0.5;
      tolerance.fit_movshon.R2 = 0.5;
      tolerance.fit_movshon.bandwidth = 0.5;
      tolerance.fit_movshon_c.parameters = 0.5;
      tolerance.fit_movshon_c.values = 0.5;
      tolerance.fit_movshon_c.fit = 0.5;
      tolerance.fit_movshon_c.L50 = 0.5;
      tolerance.fit_movshon_c.Pref = 0.5;
      tolerance.fit_movshon_c.H50 = 0.5;
      tolerance.fit_movshon_c.R2 = 0.5;
      tolerance.fit_movshon_c.bandwidth = 0.5;
      tolerance.fit_spline.values = 0.5;
      tolerance.fit_spline.fit = 0.5;
      tolerance.fit_spline.L50 = 0.5;
      tolerance.fit_spline.Pref = 0.5;
      tolerance.fit_spline.H50 = 0.5;
      tolerance.fit_spline.bandwidth = 0.5;
      tolerance.fit_gausslog.parameters = 0.5;
      tolerance.fit_gausslog.values = 0.5;
      tolerance.fit_gausslog.fit = 0.5;
      tolerance.fit_gausslog.L50 = 0.5;
      tolerance.fit_gausslog.Pref = 0.5;
      tolerance.fit_gausslog.H50 = 0.5;
      tolerance.fit_gausslog.bandwidth = 0.5;

    case 'high_noise',

      tolerance.tuning_curve.temporal_frequency = 0.5;    
      tolerance.tuning_curve.mean = 0.5;
      tolerance.tuning_curve.stddev = 0.5;
      tolerance.tuning_curve.stderr = 0.5;
      tolerance.tuning_curve.individual = 0.5;
      tolerance.tuning_curve.control_mean = 0.5;
      tolerance.tuning_curve.control_stddev = 0.5;
      tolerance.tuning_curve.control_stderr = 0.5;
      tolerance.tuning_curve.control_mean_stddev = 0.5;
      tolerance.tuning_curve.control_mean_stderr = 0.5;
      tolerance.significance.visual_response_anova_p = 0.5;
      tolerance.significance.across_stimuli_anova_p = 0.5;
      tolerance.fitless.L50 = 0.5;
      tolerance.fitless.Pref = 0.5;
      tolerance.fitless.H50 = 0.5;
      tolerance.fitless.bandwidth = 0.5;
      tolerance.fitless.low_pass_index = 0.5;
      tolerance.fitless.high_pass_index = 0.5;
      tolerance.fit_dog.parameters, = 0.5;
      tolerance.fit_dog.values = 0.5;
      tolerance.fit_dog.fit = 0.5;
      tolerance.fit_dog.L50 = 0.5;
      tolerance.fit_dog.Pref = 0.5;
      tolerance.fit_dog.H50 = 0.5;
      tolerance.fit_dog.bandwidth = 0.5;
      tolerance.fit_movshon.parameters = 0.5;
      tolerance.fit_movshon.values = 0.5;
      tolerance.fit_movshon.fit = 0.5;
      tolerance.fit_movshon.L50 = 0.5;
      tolerance.fit_movshon.Pref = 0.5;
      tolerance.fit_movshon.H50 = 0.5;
      tolerance.fit_movshon.R2 = 0.5;
      tolerance.fit_movshon.bandwidth = 0.5;
      tolerance.fit_movshon_c.parameters = 0.5;
      tolerance.fit_movshon_c.values = 0.5;
      tolerance.fit_movshon_c.fit = 0.5;
      tolerance.fit_movshon_c.L50 = 0.5;
      tolerance.fit_movshon_c.Pref = 0.5;
      tolerance.fit_movshon_c.H50 = 0.5;
      tolerance.fit_movshon_c.R2 = 0.5;
      tolerance.fit_movshon_c.bandwidth = 0.5;
      tolerance.fit_spline.values = 0.5;
      tolerance.fit_spline.fit = 0.5;
      tolerance.fit_spline.L50 = 0.5;
      tolerance.fit_spline.Pref = 0.5;
      tolerance.fit_spline.H50 = 0.5;
      tolerance.fit_spline.bandwidth = 0.5;
      tolerance.fit_gausslog.parameters = 0.5;
      tolerance.fit_gausslog.values = 0.5;
      tolerance.fit_gausslog.fit = 0.5;
      tolerance.fit_gausslog.L50 = 0.5;
      tolerance.fit_gausslog.Pref = 0.5;
      tolerance.fit_gausslog.H50 = 0.5;
      tolerance.fit_gausslog.bandwidth = 0.5;

    otherwise,
       error(['Unknown scope ' scope '.']);
end;

% start comparison

doc_e = document_expected.document_properties.temporal_frequency_tuning;
doc_a = document_actual.document_properties.temporal_frequency_tuning;

% Comparing properties
%   Response Units
%   Response Type

properties_rmatch = strcmpi(char(doc_e.properties.response_units), char(doc_a.properties.response_units));
if ~properties_rmatch
    b_(1) = 0;
    errormsg_{1} = ['Expected response units of ' doc_e.properties.response_units ' but observed ' doc_a.properties.response_units];
end

properties_match = strcmpi(char(doc_e.properties.response_type), char(doc_a.properties.response_type));
if ~properties_match
   b_(2) = 0;
   errormsg_{2} = ['Expected response type of ' doc_e.properties.response_type ' but observed ' doc_a.properties.response_type];
end

% Comparing tuning_curve
%    temporal_frequency
%    mean
%    stddev
%    stderr
%    individual
%    control_mean
%    control_stddev
%    control_stderr
%    control_mean_stddev
%    control_mean_stderr

[b_(3),errormsg_{3}] = ndi.test.values_within_tolerance(doc_e.tuning_curve.temporal_frequency, doc_a.tuning_curve.temporal_frequency, tolerance.tuning_curve.temporal_frequency, 'tuning_curve temporal_frequency');
[b_(4),errormsg_{4}] = ndi.test.values_within_tolerance(doc_e.tuning_curve.mean, doc_a.tuning_curve.mean, tolerance.tuning_curve.mean, 'tuning_curve mean');
[b_(5),errormsg_{5}] = ndi.test.values_within_tolerance(doc_e.tuning_curve.stddev, doc_a.tuning_curve.stddev, tolerance.tuning_curve.stddev, 'tuning_curve stddev');
[b_(6),errormsg_{6}] = ndi.test.values_within_tolerance(doc_e.tuning_curve.stderr, doc_a.tuning_curve.stderr, tolerance.tuning_curve.stderr, 'tuning_curve stderr');
[b_(7),errormsg_{7}] = ndi.test.values_within_tolerance(doc_e.tuning_curve.individual, doc_a.tuning_curve.individual, tolerance.tuning_curve.individual, 'tuning_curve individual');
[b_(8),errormsg_{8}] = ndi.test.values_within_tolerance(doc_e.tuning_curve.control_mean, doc_a.tuning_curve.control_mean, tolerance.tuning_curve.control_mean, 'tuning_curve control_mean');
[b_(9),errormsg_{9}] = ndi.test.values_within_tolerance(doc_e.tuning_curve.control_stddev, doc_a.tuning_curve.control_stddev, tolerance.tuning_curve.control_stddev, 'tuning_curve control_stddev');
[b_(10),errormsg_{10}] = ndi.test.values_within_tolerance(doc_e.tuning_curve.control_stderr, doc_a.tuning_curve.control_stderr, tolerance.tuning_curve.control_stderr, 'tuning_curve control_stderr');
[b_(11),errormsg_{11}] = ndi.test.values_within_tolerance(doc_e.tuning_curve.control_mean_stddev, doc_a.tuning_curve.control_mean_stddev, tolerance.tuning_curve.control_mean_stddev, 'tuning_curve control_mean_stddev');
[b_(12),errormsg_{12}] = ndi.test.values_within_tolerance(doc_e.tuning_curve.control_mean_stderr, doc_a.tuning_curve.control_mean_stderr, tolerance.tuning_curve.control_mean_stderr, 'tuning_curve control_mean_stderr');

% Comparing significance
%    visual_response_anova_p
%    across_stimuli_anova_p

[b_(13),errormsg_{13}] = ndi.test.values_within_tolerance(doc_e.significance.visual_response_anova_p, doc_a.significance.visual_response_anova_p, tolerance.significance.visual_response_anova_p, 'visual_response_anova_p');
[b_(14),errormsg_{14}] = ndi.test.values_within_tolerance(doc_e.significance.across_stimuli_anova_p, doc_a.significance.across_stimuli_anova_p, tolerance.significance.across_stimuli_anova_p, 'across_stimuli_anova_p');

% Comparing fitless
%    L50
%    Pref
%    H50
%    bandwidth
%    low_pass_index
%    high_pass_index

[b_(15),errormsg_{15}] = ndi.test.values_within_tolerance(doc_e.fitless.L50, doc_a.fitless.L50, tolerance.fitless.L50, 'fitless L50');
[b_(16),errormsg_{16}] = ndi.test.values_within_tolerance(doc_e.fitless.Pref, doc_a.fitless.Pref, tolerance.fitless.Pref, 'fitless Pref');
[b_(17),errormsg_{17}] = ndi.test.values_within_tolerance(doc_e.fitless.H50, doc_a.fitless.H50, tolerance.fitless.H50, 'fitless H50');
[b_(18),errormsg_{18}] = ndi.test.values_within_tolerance(doc_e.fitless.bandwidth, doc_a.fitless.bandwidth, tolerance.fitless.bandwidth, 'fitless bandwidth');
[b_(19),errormsg_{19}] = ndi.test.values_within_tolerance(doc_e.fitless.low_pass_index, doc_a.fitless.low_pass_index, tolerance.fitless.low_pass_index, 'fitless low_pass_index');
[b_(20),errormsg_{20}] = ndi.test.values_within_tolerance(doc_e.fitless.high_pass_index, doc_a.fitless.high_pass_index, tolerance.fitless.high_pass_index, 'fitless high_pass_index');

% Comparing fit_dog
%    parameters
%    values
%    fit
%    L50
%    Pref
%    H50
%    bandwidth

[b_(21),errormsg_{21}] = ndi.test.values_within_tolerance(doc_e.fit_dog.parameters, doc_a.fit_dog.parameters, tolerance.fit_dog.parameters, 'fit_dog parameters');
[b_(22),errormsg_{22}] = ndi.test.values_within_tolerance(doc_e.fit_dog.values, doc_a.fit_dog.values, tolerance.fit_dog.values, 'fit_dog values');
[b_(23),errormsg_{23}] = ndi.test.values_within_tolerance(doc_e.fit_dog.fit, doc_a.fit_dog.fit, tolerance.fit_dog.fit, 'fit_dog fit');
[b_(24),errormsg_{24}] = ndi.test.values_within_tolerance(doc_e.fit_dog.L50, doc_a.fit_dog.L50, tolerance.fit_dog.L50, 'fit_dog L50');
[b_(25),errormsg_{25}] = ndi.test.values_within_tolerance(doc_e.fit_dog.Pref, doc_a.fit_dog.Pref, tolerance.fit_dog.Pref, 'fit_dog Pref');
[b_(26),errormsg_{26}] = ndi.test.values_within_tolerance(doc_e.fit_dog.H50, doc_a.fit_dog.H50, tolerance.fit_dog.H50, 'fit_dog H50');
[b_(27),errormsg_{27}] = ndi.test.values_within_tolerance(doc_e.fit_dog.bandwidth, doc_a.fit_dog.bandwidth, tolerance.fit_dog.bandwidth, 'fit_dog bandwidth');

% Comparing fit_movshon
%    parameters
%    values
%    fit
%    L50
%    Pref
%    H50
%    R2
%    bandwidth

[b_(28),errormsg_{28}] = ndi.test.values_within_tolerance(doc_e.fit_movshon.parameters, doc_a.fit_movshon.parameters, tolerance.fit_movshon.parameters, 'fit_movshon parameters');
[b_(29),errormsg_{29}] = ndi.test.values_within_tolerance(doc_e.fit_movshon.values, doc_a.fit_movshon.values, tolerance.fit_movshon.values, 'fit_movshon values');
[b_(30),errormsg_{30}] = ndi.test.values_within_tolerance(doc_e.fit_movshon.fit, doc_a.fit_movshon.fit, tolerance.fit_movshon.fit, 'fit_movshon fit');
[b_(31),errormsg_{31}] = ndi.test.values_within_tolerance(doc_e.fit_movshon.L50, doc_a.fit_movshon.L50, tolerance.fit_movshon.L50, 'fit_movshon L50');
[b_(32),errormsg_{32}] = ndi.test.values_within_tolerance(doc_e.fit_movshon.Pref, doc_a.fit_movshon.Pref, tolerance.fit_movshon.Pref, 'fit_movshon Pref');
[b_(33),errormsg_{33}] = ndi.test.values_within_tolerance(doc_e.fit_movshon.H50, doc_a.fit_movshon.H50, tolerance.fit_movshon.H50, 'fit_movshon H50');
[b_(34),errormsg_{34}] = ndi.test.values_within_tolerance(doc_e.fit_movshon.R2, doc_a.fit_movshon.R2, tolerance.fit_movshon.R2, 'fit_movshon R2');
[b_(35),errormsg_{35}] = ndi.test.values_within_tolerance(doc_e.fit_movshon.bandwidth, doc_a.fit_movshon.bandwidth, tolerance.fit_movshon.bandwidth, 'fit_movshon bandwidth');

% Comparing fit_movshon_c
%    parameters
%    values
%    fit
%    L50
%    Pref
%    H50
%    R2
%    bandwidth

[b_(36),errormsg_{36}] = ndi.test.values_within_tolerance(doc_e.fit_movshon_c.parameters, doc_a.fit_movshon_c.parameters, tolerance.fit_movshon_c.parameters, 'fit_movshon_c parameters');
[b_(37),errormsg_{37}] = ndi.test.values_within_tolerance(doc_e.fit_movshon_c.values, doc_a.fit_movshon_c.values, tolerance.fit_movshon_c.values, 'fit_movshon_c values');
[b_(38),errormsg_{38}] = ndi.test.values_within_tolerance(doc_e.fit_movshon_c.fit, doc_a.fit_movshon_c.fit, tolerance.fit_movshon_c.fit, 'fit_movshon_c fit');
[b_(39),errormsg_{39}] = ndi.test.values_within_tolerance(doc_e.fit_movshon_c.L50, doc_a.fit_movshon_c.L50, tolerance.fit_movshon_c.L50, 'fit_movshon_c L50');
[b_(40),errormsg_{40}] = ndi.test.values_within_tolerance(doc_e.fit_movshon_c.Pref, doc_a.fit_movshon_c.Pref, tolerance.fit_movshon_c.Pref, 'fit_movshon_c Pref');
[b_(41),errormsg_{41}] = ndi.test.values_within_tolerance(doc_e.fit_movshon_c.H50, doc_a.fit_movshon_c.H50, tolerance.fit_movshon_c.H50, 'fit_movshon_c H50');
[b_(42),errormsg_{42}] = ndi.test.values_within_tolerance(doc_e.fit_movshon_c.R2, doc_a.fit_movshon_c.R2, tolerance.fit_movshon_c.R2, 'fit_movshon_c R2');
[b_(43),errormsg_{43}] = ndi.test.values_within_tolerance(doc_e.fit_movshon_c.bandwidth, doc_a.fit_movshon_c.bandwidth, tolerance.fit_movshon_c.bandwidth, 'fit_movshon_c bandwidth');

% Comparing fit_spline
%    values
%    fit
%    L50
%    Pref
%    H50
%    bandwidth

[b_(44),errormsg_{44}] = ndi.test.values_within_tolerance(doc_e.fit_spline.values, doc_a.fit_spline.values, tolerance.fit_spline.values, 'fit_spline values');
[b_(45),errormsg_{45}] = ndi.test.values_within_tolerance(doc_e.fit_spline.fit, doc_a.fit_spline.fit, tolerance.fit_spline.fit, 'fit_spline fit');
[b_(46),errormsg_{46}] = ndi.test.values_within_tolerance(doc_e.fit_spline.L50, doc_a.fit_spline.L50, tolerance.fit_spline.L50, 'fit_spline L50');
[b_(47),errormsg_{47}] = ndi.test.values_within_tolerance(doc_e.fit_spline.Pref, doc_a.fit_spline.Pref, tolerance.fit_spline.Pref, 'fit_spline Pref');
[b_(48),errormsg_{48}] = ndi.test.values_within_tolerance(doc_e.fit_spline.H50, doc_a.fit_spline.H50, tolerance.fit_spline.H50, 'fit_spline H50');
[b_(49),errormsg_{49}] = ndi.test.values_within_tolerance(doc_e.fit_spline.bandwidth, doc_a.fit_spline.bandwidth, tolerance.fit_spline.bandwidth, 'fit_spline bandwidth');

% Comparing fit_gausslog
%    parameters
%    values
%    fit
%    L50
%    Pref
%    H50
%    bandwidth

[b_(50),errormsg_{50}] = ndi.test.values_within_tolerance(doc_e.fit_gausslog.parameters, doc_a.fit_gausslog.parameters, tolerance.fit_gausslog.parameters, 'fit_gausslog parameters');
[b_(51),errormsg_{51}] = ndi.test.values_within_tolerance(doc_e.fit_gausslog.values, doc_a.fit_gausslog.values, tolerance.fit_gausslog.values, 'fit_gausslog values');
[b_(52),errormsg_{52}] = ndi.test.values_within_tolerance(doc_e.fit_gausslog.fit, doc_a.fit_gausslog.fit, tolerance.fit_gausslog.fit, 'fit_gausslog fit');
[b_(53),errormsg_{53}] = ndi.test.values_within_tolerance(doc_e.fit_gausslog.L50, doc_a.fit_gausslog.L50, tolerance.fit_gausslog.L50, 'fit_gausslog L50');
[b_(54),errormsg_{54}] = ndi.test.values_within_tolerance(doc_e.fit_gausslog.Pref, doc_a.fit_gausslog.Pref, tolerance.fit_gausslog.Pref, 'fit_gausslog Pref');
[b_(55),errormsg_{55}] = ndi.test.values_within_tolerance(doc_e.fit_gausslog.H50, doc_a.fit_gausslog.H50, tolerance.fit_gausslog.H50, 'fit_gausslog H50');
[b_(56),errormsg_{56}] = ndi.test.values_within_tolerance(doc_e.fit_gausslog.bandwidth, doc_a.fit_gausslog.bandwidth, tolerance.fit_gausslog.bandwidth, 'fit_gausslog bandwidth');

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
