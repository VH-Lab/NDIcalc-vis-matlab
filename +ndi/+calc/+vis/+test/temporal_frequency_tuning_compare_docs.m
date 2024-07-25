function [b_, errormsg_] = temporal_frequency_tuning_compare_docs(document_expected, document_actual, scope)
% Comparing temporal_frequency_tuning
%
% [B, ERRORMSG] = ndi.calc.vis.test.temporal_frequency_tuning_compare_docs(DOC_EXPECTED, DOC_ACTUAL, SCOPE)
%
%

% Initialize b_ as a row vector of ones for tracking comparison results
% Initialize errormsg_ as an empty cell array to hold error messages

b_ = ones(1,98);
errormsg_ = cell(1,98);

% establish scope-dependent tolerances
switch(scope),
    case 'standard',
    
      tolerance.tuning_curve.temporal_frequency = 1e-6; %tuning curve fields should be equivalent to the golden file when there's no noised added   
      tolerance.tuning_curve.mean = 1e-6;
      tolerance.tuning_curve.stddev = 1e-6;
      tolerance.tuning_curve.stderr = 1e-6;
      tolerance.tuning_curve.individual = 1e-6;
      tolerance.tuning_curve.control_mean = 1e-6;
      tolerance.tuning_curve.control_stddev = 1e-6;
      tolerance.tuning_curve.control_stderr = 1e-6;
      tolerance.tuning_curve.control_mean_stddev = 1e-6;
      tolerance.tuning_curve.control_mean_stderr = 1e-6;
      tolerance.significance.visual_response_anova_p = 1e-6;
      tolerance.significance.across_stimuli_anova_p = 1e-6;
      tolerance.fitless.L50 = 1e-6; %fitless should have no variability either
      tolerance.fitless.Pref = 1e-6;
      tolerance.fitless.H50 = 1e-6;
      tolerance.fitless.bandwidth = 1e-6;
      tolerance.fitless.low_pass_index = 1e-6;
      tolerance.fitless.high_pass_index = 1e-6;
      % DOG function: y = A1*exp(-x.^2/(2*B1^2)) - A2*exp(-x.^2/(2*B2^2))
      tolerance.fit_dog.parameters(1) = 0.5; %A1
      tolerance.fit_dog.parameters(2) = 0.5; %B1
      tolerance.fit_dog.parameters(3) = 0.5; %A2
      tolerance.fit_dog.parameters(4) = 0.5; %B2
      tolerance.fit_dog.values = 1e-6; 
      tolerance.fit_dog.fit = 0.5; %0 to ~50
      tolerance.fit_dog.L50 = 0.5; %0.02 to 60
      tolerance.fit_dog.Pref = 0.5; %0.02 to 60
      tolerance.fit_dog.H50 = 0.5; %0.02 to 60
      tolerance.fit_dog.bandwidth = 0.5; %0.02 to 60
      %Movshon function: R(f) = k * exp(-(f./fc).^2) ./ (1+(fh./f).^B)
      tolerance.fit_movshon.parameters(1) = 0.5;%   k - P(1) - scaling factor
      tolerance.fit_movshon.parameters(2) = 0.5;%   fc - P(2) - characteristic temporal frequency
      tolerance.fit_movshon.parameters(3) = 0.5;%   fh - P(3) - corner frequency of the low-frequency limb
      tolerance.fit_movshon.parameters(4) = 1;%   B - P(4) - slope of the low-frequency limb (no limit)
      tolerance.fit_movshon.values = 1e-6;
      tolerance.fit_movshon.fit = 0.5; %0 to ~50
      tolerance.fit_movshon.L50 = 0.5; %0.02 to 60
      tolerance.fit_movshon.Pref = 0.5; %0.02 to 60
      tolerance.fit_movshon.H50 = 0.5; %0.02 to 60
      tolerance.fit_movshon.R2 = 0.05; %0 to 1
      tolerance.fit_movshon.bandwidth = 0.5; %0.02 to 60
      %Movshon_c function: R(f) = k * exp(-(f./fc).^2) ./ (1+(fh./f).^B) + C
      tolerance.fit_movshon_c.parameters(1) = 0.5;%   k - P(1) - scaling factor
      tolerance.fit_movshon_c.parameters(2) = 0.5;%   fc - P(2) - characteristic temporal frequency
      tolerance.fit_movshon_c.parameters(3) = 0.5;%   fh - P(3) - corner frequency of the low-frequency limb
      tolerance.fit_movshon_c.parameters(4) = 1;%   B - P(4) - slope of the low-frequency limb (no limit)
      tolerance.fit_movshon_c.parameters(5) = 0.2;%   c - P(5) - constant term (0 to ~10)
      tolerance.fit_movshon_c.values = 1e-6;
      tolerance.fit_movshon_c.fit = 0.5; %0 to ~50
      tolerance.fit_movshon_c.L50 = 0.5; %0.02 to 60
      tolerance.fit_movshon_c.Pref = 0.5; %0.02 to 60
      tolerance.fit_movshon_c.H50 = 0.5; %0.02 to 60
      tolerance.fit_movshon_c.R2 = 0.05; %0 to 1
      tolerance.fit_movshon_c.bandwidth = 0.5; %0.02 to 60
      tolerance.fit_spline.values = 1e-6;
      tolerance.fit_spline.fit = 0.5; %0 to ~50
      tolerance.fit_spline.L50 = 0.5; %0.02 to 60
      tolerance.fit_spline.Pref = 0.5; %0.02 to 60
      tolerance.fit_spline.H50 = 0.5; %0.02 to 60
      tolerance.fit_spline.bandwidth = 0.5; %0.02 to 60
      %Gausslog function: Y = a+b*exp(-((log10(x)-log10(c)).^2/(2*d^2)))
      tolerance.fit_gausslog.parameters(1) = 0.5; %a is an offset parameter
      tolerance.fit_gausslog.parameters(2) = 0.5; %b is a height parameter above the offset
      tolerance.fit_gausslog.parameters(3) = 0.5; %c is the peak location
      tolerance.fit_gausslog.parameters(4) = 0.5; %d is the width
      tolerance.fit_gausslog.values = 1e-6;
      tolerance.fit_gausslog.fit = 0.5; 
      tolerance.fit_gausslog.L50 = 0.5; 
      tolerance.fit_gausslog.Pref = 0.5; 
      tolerance.fit_gausslog.H50 = 0.5; 
      tolerance.fit_gausslog.bandwidth = 0.5; 
      tolerance.abs.fitless.L50 = 1e-6; %fitless should have no variability either
      tolerance.abs.fitless.Pref = 1e-6;
      tolerance.abs.fitless.H50 = 1e-6;
      tolerance.abs.fitless.bandwidth = 1e-6;
      tolerance.abs.fitless.low_pass_index = 1e-6;
      tolerance.abs.fitless.high_pass_index = 1e-6;
      % DOG function: y = A1*exp(-x.^2/(2*B1^2)) - A2*exp(-x.^2/(2*B2^2))
      tolerance.abs.fit_dog.parameters(1) = 0.5; %A1
      tolerance.abs.fit_dog.parameters(2) = 0.5; %B1
      tolerance.abs.fit_dog.parameters(3) = 0.5; %A2
      tolerance.abs.fit_dog.parameters(4) = 0.5; %B2
      tolerance.abs.fit_dog.values = 1e-6; 
      tolerance.abs.fit_dog.fit = 0.5; %0 to ~50
      tolerance.abs.fit_dog.L50 = 0.5; %0.02 to 60
      tolerance.abs.fit_dog.Pref = 0.5; %0.02 to 60
      tolerance.abs.fit_dog.H50 = 0.5; %0.02 to 60
      tolerance.abs.fit_dog.bandwidth = 0.5; %0.02 to 60
      %Movshon function: R(f) = k * exp(-(f./fc).^2) ./ (1+(fh./f).^B)
      tolerance.abs.fit_movshon.parameters(1) = 0.5;%   k - P(1) - scaling factor
      tolerance.abs.fit_movshon.parameters(2) = 0.5;%   fc - P(2) - characteristic temporal frequency
      tolerance.abs.fit_movshon.parameters(3) = 0.5;%   fh - P(3) - corner frequency of the low-frequency limb
      tolerance.abs.fit_movshon.parameters(4) = 1;%   B - P(4) - slope of the low-frequency limb (no limit)
      tolerance.abs.fit_movshon.values = 1e-6;
      tolerance.abs.fit_movshon.fit = 0.5; %0 to ~50
      tolerance.abs.fit_movshon.L50 = 0.5; %0.02 to 60
      tolerance.abs.fit_movshon.Pref = 0.5; %0.02 to 60
      tolerance.abs.fit_movshon.H50 = 0.5; %0.02 to 60
      tolerance.abs.fit_movshon.R2 = 0.05; %0 to 1
      tolerance.abs.fit_movshon.bandwidth = 0.5; %0.02 to 60
      %Movshon_c function: R(f) = k * exp(-(f./fc).^2) ./ (1+(fh./f).^B) + C
      tolerance.abs.fit_movshon_c.parameters(1) = 0.5;%   k - P(1) - scaling factor
      tolerance.abs.fit_movshon_c.parameters(2) = 0.5;%   fc - P(2) - characteristic temporal frequency
      tolerance.abs.fit_movshon_c.parameters(3) = 0.5;%   fh - P(3) - corner frequency of the low-frequency limb
      tolerance.abs.fit_movshon_c.parameters(4) = 1;%   B - P(4) - slope of the low-frequency limb (no limit)
      tolerance.abs.fit_movshon_c.parameters(5) = 0.2;%   c - P(5) - constant term (0 to ~10)
      tolerance.abs.fit_movshon_c.values = 1e-6;
      tolerance.abs.fit_movshon_c.fit = 0.5; %0 to ~50
      tolerance.abs.fit_movshon_c.L50 = 0.5; %0.02 to 60
      tolerance.abs.fit_movshon_c.Pref = 0.5; %0.02 to 60
      tolerance.abs.fit_movshon_c.H50 = 0.5; %0.02 to 60
      tolerance.abs.fit_movshon_c.R2 = 0.05; %0 to 1
      tolerance.abs.fit_movshon_c.bandwidth = 0.5; %0.02 to 60
      tolerance.abs.fit_spline.values = 1e-6;
      tolerance.abs.fit_spline.fit = 0.5; %0 to ~50
      tolerance.abs.fit_spline.L50 = 0.5; %0.02 to 60
      tolerance.abs.fit_spline.Pref = 0.5; %0.02 to 60
      tolerance.abs.fit_spline.H50 = 0.5; %0.02 to 60
      tolerance.abs.fit_spline.bandwidth = 0.5; %0.02 to 60
      %Gausslog function: Y = a+b*exp(-((log10(x)-log10(c)).^2/(2*d^2)))
      tolerance.abs.fit_gausslog.parameters(1) = 0.5; %a is an offset parameter
      tolerance.abs.fit_gausslog.parameters(2) = 0.5; %b is a height parameter above the offset
      tolerance.abs.fit_gausslog.parameters(3) = 0.5; %c is the peak location
      tolerance.abs.fit_gausslog.parameters(4) = 0.5; %d is the width
      tolerance.abs.fit_gausslog.values = 1e-6;
      tolerance.abs.fit_gausslog.fit = 0.5; 
      tolerance.abs.fit_gausslog.L50 = 0.5; 
      tolerance.abs.fit_gausslog.Pref = 0.5; 
      tolerance.abs.fit_gausslog.H50 = 0.5; 
      tolerance.abs.fit_gausslog.bandwidth = 0.5; 

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
      tolerance.fit_dog.parameters = 0.5;
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
      tolerance.abs.fitless.L50 = 0.5;
      tolerance.abs.fitless.Pref = 0.5;
      tolerance.abs.fitless.H50 = 0.5;
      tolerance.abs.fitless.bandwidth = 0.5;
      tolerance.abs.fitless.low_pass_index = 0.5;
      tolerance.abs.fitless.high_pass_index = 0.5;
      tolerance.abs.fit_dog.parameters = 0.5;
      tolerance.abs.fit_dog.values = 0.5;
      tolerance.abs.fit_dog.fit = 0.5;
      tolerance.abs.fit_dog.L50 = 0.5;
      tolerance.abs.fit_dog.Pref = 0.5;
      tolerance.abs.fit_dog.H50 = 0.5;
      tolerance.abs.fit_dog.bandwidth = 0.5;
      tolerance.abs.fit_movshon.parameters = 0.5;
      tolerance.abs.fit_movshon.values = 0.5;
      tolerance.abs.fit_movshon.fit = 0.5;
      tolerance.abs.fit_movshon.L50 = 0.5;
      tolerance.abs.fit_movshon.Pref = 0.5;
      tolerance.abs.fit_movshon.H50 = 0.5;
      tolerance.abs.fit_movshon.R2 = 0.5;
      tolerance.abs.fit_movshon.bandwidth = 0.5;
      tolerance.abs.fit_movshon_c.parameters = 0.5;
      tolerance.abs.fit_movshon_c.values = 0.5;
      tolerance.abs.fit_movshon_c.fit = 0.5;
      tolerance.abs.fit_movshon_c.L50 = 0.5;
      tolerance.abs.fit_movshon_c.Pref = 0.5;
      tolerance.abs.fit_movshon_c.H50 = 0.5;
      tolerance.abs.fit_movshon_c.R2 = 0.5;
      tolerance.abs.fit_movshon_c.bandwidth = 0.5;
      tolerance.abs.fit_spline.values = 0.5;
      tolerance.abs.fit_spline.fit = 0.5;
      tolerance.abs.fit_spline.L50 = 0.5;
      tolerance.abs.fit_spline.Pref = 0.5;
      tolerance.abs.fit_spline.H50 = 0.5;
      tolerance.abs.fit_spline.bandwidth = 0.5;
      tolerance.abs.fit_gausslog.parameters = 0.5;
      tolerance.abs.fit_gausslog.values = 0.5;
      tolerance.abs.fit_gausslog.fit = 0.5;
      tolerance.abs.fit_gausslog.L50 = 0.5;
      tolerance.abs.fit_gausslog.Pref = 0.5;
      tolerance.abs.fit_gausslog.H50 = 0.5;
      tolerance.abs.fit_gausslog.bandwidth = 0.5;

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
      tolerance.fit_dog.parameters = 0.5;
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
      tolerance.abs.fitless.L50 = 0.5;
      tolerance.abs.fitless.Pref = 0.5;
      tolerance.abs.fitless.H50 = 0.5;
      tolerance.abs.fitless.bandwidth = 0.5;
      tolerance.abs.fitless.low_pass_index = 0.5;
      tolerance.abs.fitless.high_pass_index = 0.5;
      tolerance.abs.fit_dog.parameters = 0.5;
      tolerance.abs.fit_dog.values = 0.5;
      tolerance.abs.fit_dog.fit = 0.5;
      tolerance.abs.fit_dog.L50 = 0.5;
      tolerance.abs.fit_dog.Pref = 0.5;
      tolerance.abs.fit_dog.H50 = 0.5;
      tolerance.abs.fit_dog.bandwidth = 0.5;
      tolerance.abs.fit_movshon.parameters = 0.5;
      tolerance.abs.fit_movshon.values = 0.5;
      tolerance.abs.fit_movshon.fit = 0.5;
      tolerance.abs.fit_movshon.L50 = 0.5;
      tolerance.abs.fit_movshon.Pref = 0.5;
      tolerance.abs.fit_movshon.H50 = 0.5;
      tolerance.abs.fit_movshon.R2 = 0.5;
      tolerance.abs.fit_movshon.bandwidth = 0.5;
      tolerance.abs.fit_movshon_c.parameters = 0.5;
      tolerance.abs.fit_movshon_c.values = 0.5;
      tolerance.abs.fit_movshon_c.fit = 0.5;
      tolerance.abs.fit_movshon_c.L50 = 0.5;
      tolerance.abs.fit_movshon_c.Pref = 0.5;
      tolerance.abs.fit_movshon_c.H50 = 0.5;
      tolerance.abs.fit_movshon_c.R2 = 0.5;
      tolerance.abs.fit_movshon_c.bandwidth = 0.5;
      tolerance.abs.fit_spline.values = 0.5;
      tolerance.abs.fit_spline.fit = 0.5;
      tolerance.abs.fit_spline.L50 = 0.5;
      tolerance.abs.fit_spline.Pref = 0.5;
      tolerance.abs.fit_spline.H50 = 0.5;
      tolerance.abs.fit_spline.bandwidth = 0.5;
      tolerance.abs.fit_gausslog.parameters = 0.5;
      tolerance.abs.fit_gausslog.values = 0.5;
      tolerance.abs.fit_gausslog.fit = 0.5;
      tolerance.abs.fit_gausslog.L50 = 0.5;
      tolerance.abs.fit_gausslog.Pref = 0.5;
      tolerance.abs.fit_gausslog.H50 = 0.5;
      tolerance.abs.fit_gausslog.bandwidth = 0.5;

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
%    parameters 1-4
%    values
%    fit
%    L50
%    Pref
%    H50
%    bandwidth

[b_(21),errormsg_{21}] = ndi.test.values_within_tolerance(doc_e.fit_dog.parameters(1), doc_a.fit_dog.parameters(1), tolerance.fit_dog.parameters(1), 'fit_dog parameter 1');
[b_(end+1),errormsg_{end+1}] = ndi.test.values_within_tolerance(doc_e.fit_dog.parameters(2), doc_a.fit_dog.parameters(2), tolerance.fit_dog.parameters(2), 'fit_dog parameter 2');
[b_(end+1),errormsg_{end+1}] = ndi.test.values_within_tolerance(doc_e.fit_dog.parameters(3), doc_a.fit_dog.parameters(3), tolerance.fit_dog.parameters(3), 'fit_dog parameter 3');
[b_(end+1),errormsg_{end+1}] = ndi.test.values_within_tolerance(doc_e.fit_dog.parameters(4), doc_a.fit_dog.parameters(4), tolerance.fit_dog.parameters(4), 'fit_dog parameter 4');
[b_(22),errormsg_{22}] = ndi.test.values_within_tolerance(doc_e.fit_dog.values, doc_a.fit_dog.values, tolerance.fit_dog.values, 'fit_dog values');
[b_(23),errormsg_{23}] = ndi.test.values_within_tolerance(doc_e.fit_dog.fit, doc_a.fit_dog.fit, tolerance.fit_dog.fit, 'fit_dog fit');
[b_(24),errormsg_{24}] = ndi.test.values_within_tolerance(doc_e.fit_dog.L50, doc_a.fit_dog.L50, tolerance.fit_dog.L50, 'fit_dog L50');
[b_(25),errormsg_{25}] = ndi.test.values_within_tolerance(doc_e.fit_dog.Pref, doc_a.fit_dog.Pref, tolerance.fit_dog.Pref, 'fit_dog Pref');
[b_(26),errormsg_{26}] = ndi.test.values_within_tolerance(doc_e.fit_dog.H50, doc_a.fit_dog.H50, tolerance.fit_dog.H50, 'fit_dog H50');
[b_(27),errormsg_{27}] = ndi.test.values_within_tolerance(doc_e.fit_dog.bandwidth, doc_a.fit_dog.bandwidth, tolerance.fit_dog.bandwidth, 'fit_dog bandwidth');

% Comparing fit_movshon
%    parameters 1-4
%    values
%    fit
%    L50
%    Pref
%    H50
%    R2
%    bandwidth

[b_(28),errormsg_{28}] = ndi.test.values_within_tolerance(doc_e.fit_movshon.parameters(1), doc_a.fit_movshon.parameters(1), tolerance.fit_movshon.parameters(1), 'fit_movshon parameter 1');
[b_(end+1),errormsg_{end+1}] = ndi.test.values_within_tolerance(doc_e.fit_movshon.parameters(2), doc_a.fit_movshon.parameters(2), tolerance.fit_movshon.parameters(2), 'fit_movshon parameter 2');
[b_(end+1),errormsg_{end+1}] = ndi.test.values_within_tolerance(doc_e.fit_movshon.parameters(3), doc_a.fit_movshon.parameters(3), tolerance.fit_movshon.parameters(3), 'fit_movshon parameter 3');
[b_(end+1),errormsg_{end+1}] = ndi.test.values_within_tolerance(doc_e.fit_movshon.parameters(4), doc_a.fit_movshon.parameters(4), tolerance.fit_movshon.parameters(4), 'fit_movshon parameter 4');
[b_(29),errormsg_{29}] = ndi.test.values_within_tolerance(doc_e.fit_movshon.values, doc_a.fit_movshon.values, tolerance.fit_movshon.values, 'fit_movshon values');
[b_(30),errormsg_{30}] = ndi.test.values_within_tolerance(doc_e.fit_movshon.fit, doc_a.fit_movshon.fit, tolerance.fit_movshon.fit, 'fit_movshon fit');
[b_(31),errormsg_{31}] = ndi.test.values_within_tolerance(doc_e.fit_movshon.L50, doc_a.fit_movshon.L50, tolerance.fit_movshon.L50, 'fit_movshon L50');
[b_(32),errormsg_{32}] = ndi.test.values_within_tolerance(doc_e.fit_movshon.Pref, doc_a.fit_movshon.Pref, tolerance.fit_movshon.Pref, 'fit_movshon Pref');
[b_(33),errormsg_{33}] = ndi.test.values_within_tolerance(doc_e.fit_movshon.H50, doc_a.fit_movshon.H50, tolerance.fit_movshon.H50, 'fit_movshon H50');
[b_(34),errormsg_{34}] = ndi.test.values_within_tolerance(doc_e.fit_movshon.R2, doc_a.fit_movshon.R2, tolerance.fit_movshon.R2, 'fit_movshon R2');
[b_(35),errormsg_{35}] = ndi.test.values_within_tolerance(doc_e.fit_movshon.bandwidth, doc_a.fit_movshon.bandwidth, tolerance.fit_movshon.bandwidth, 'fit_movshon bandwidth');

% Comparing fit_movshon_c
%    parameters 1-5
%    values
%    fit
%    L50
%    Pref
%    H50
%    R2
%    bandwidth

[b_(36),errormsg_{36}] = ndi.test.values_within_tolerance(doc_e.fit_movshon_c.parameters(1), doc_a.fit_movshon_c.parameters(1), tolerance.fit_movshon_c.parameters(1), 'fit_movshon_c parameter 1');
[b_(end+1),errormsg_{end+1}] = ndi.test.values_within_tolerance(doc_e.fit_movshon_c.parameters(2), doc_a.fit_movshon_c.parameters(2), tolerance.fit_movshon_c.parameters(2), 'fit_movshon_c parameter 2');
[b_(end+1),errormsg_{end+1}] = ndi.test.values_within_tolerance(doc_e.fit_movshon_c.parameters(3), doc_a.fit_movshon_c.parameters(3), tolerance.fit_movshon_c.parameters(3), 'fit_movshon_c parameter 3');
[b_(end+1),errormsg_{end+1}] = ndi.test.values_within_tolerance(doc_e.fit_movshon_c.parameters(4), doc_a.fit_movshon_c.parameters(4), tolerance.fit_movshon_c.parameters(4), 'fit_movshon_c parameter 4');
[b_(end+1),errormsg_{end+1}] = ndi.test.values_within_tolerance(doc_e.fit_movshon_c.parameters(5), doc_a.fit_movshon_c.parameters(5), tolerance.fit_movshon_c.parameters(5), 'fit_movshon_c parameter 5');
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
%    parameters 1-4
%    values
%    fit
%    L50
%    Pref
%    H50
%    bandwidth

[b_(50),errormsg_{50}] = ndi.test.values_within_tolerance(doc_e.fit_gausslog.parameters(1), doc_a.fit_gausslog.parameters(1), tolerance.fit_gausslog.parameters(1), 'fit_gausslog parameter 1');
[b_(end+1),errormsg_{end+1}] = ndi.test.values_within_tolerance(doc_e.fit_gausslog.parameters(2), doc_a.fit_gausslog.parameters(2), tolerance.fit_gausslog.parameters(2), 'fit_gausslog parameter 2');
[b_(end+1),errormsg_{end+1}] = ndi.test.values_within_tolerance(doc_e.fit_gausslog.parameters(3), doc_a.fit_gausslog.parameters(3), tolerance.fit_gausslog.parameters(3), 'fit_gausslog parameter 3');
[b_(end+1),errormsg_{end+1}] = ndi.test.values_within_tolerance(doc_e.fit_gausslog.parameters(4), doc_a.fit_gausslog.parameters(4), tolerance.fit_gausslog.parameters(4), 'fit_gausslog parameter 4');
[b_(51),errormsg_{51}] = ndi.test.values_within_tolerance(doc_e.fit_gausslog.values, doc_a.fit_gausslog.values, tolerance.fit_gausslog.values, 'fit_gausslog values');
[b_(52),errormsg_{52}] = ndi.test.values_within_tolerance(doc_e.fit_gausslog.fit, doc_a.fit_gausslog.fit, tolerance.fit_gausslog.fit, 'fit_gausslog fit');
[b_(53),errormsg_{53}] = ndi.test.values_within_tolerance(doc_e.fit_gausslog.L50, doc_a.fit_gausslog.L50, tolerance.fit_gausslog.L50, 'fit_gausslog L50');
[b_(54),errormsg_{54}] = ndi.test.values_within_tolerance(doc_e.fit_gausslog.Pref, doc_a.fit_gausslog.Pref, tolerance.fit_gausslog.Pref, 'fit_gausslog Pref');
[b_(55),errormsg_{55}] = ndi.test.values_within_tolerance(doc_e.fit_gausslog.H50, doc_a.fit_gausslog.H50, tolerance.fit_gausslog.H50, 'fit_gausslog H50');
[b_(56),errormsg_{56}] = ndi.test.values_within_tolerance(doc_e.fit_gausslog.bandwidth, doc_a.fit_gausslog.bandwidth, tolerance.fit_gausslog.bandwidth, 'fit_gausslog bandwidth');

% abs
%    Comparing fitless
%        L50
%        Pref
%        H50
%        bandwidth
%        low_pass_index
%        high_pass_index

[b_(57),errormsg_{57}] = ndi.test.values_within_tolerance(doc_e.abs.fitless.L50, doc_a.abs.fitless.L50, tolerance.abs.fitless.L50, 'abs fitless L50');
[b_(58),errormsg_{58}] = ndi.test.values_within_tolerance(doc_e.abs.fitless.Pref, doc_a.abs.fitless.Pref, tolerance.abs.fitless.Pref, 'abs fitless Pref');
[b_(59),errormsg_{59}] = ndi.test.values_within_tolerance(doc_e.abs.fitless.H50, doc_a.abs.fitless.H50, tolerance.abs.fitless.H50, 'abs fitless H50');
[b_(60),errormsg_{60}] = ndi.test.values_within_tolerance(doc_e.abs.fitless.bandwidth, doc_a.abs.fitless.bandwidth, tolerance.abs.fitless.bandwidth, 'abs fitless bandwidth');
[b_(61),errormsg_{61}] = ndi.test.values_within_tolerance(doc_e.abs.fitless.low_pass_index, doc_a.abs.fitless.low_pass_index, tolerance.abs.fitless.low_pass_index, 'abs fitless low_pass_index');
[b_(62),errormsg_{62}] = ndi.test.values_within_tolerance(doc_e.abs.fitless.high_pass_index, doc_a.abs.fitless.high_pass_index, tolerance.abs.fitless.high_pass_index, 'abs fitless high_pass_index');

%    Comparing fit_dog
%        parameters 1-4
%        values
%        fit
%        L50
%        Pref
%        H50
%        bandwidth

[b_(63),errormsg_{63}] = ndi.test.values_within_tolerance(doc_e.abs.fit_dog.parameters(1), doc_a.abs.fit_dog.parameters(1), tolerance.abs.fit_dog.parameters(1), 'abs fit_dog parameter 1');
[b_(end+1),errormsg_{end+1}] = ndi.test.values_within_tolerance(doc_e.abs.fit_dog.parameters(2), doc_a.abs.fit_dog.parameters(2), tolerance.abs.fit_dog.parameters(2), 'abs fit_dog parameter 2');
[b_(end+1),errormsg_{end+1}] = ndi.test.values_within_tolerance(doc_e.abs.fit_dog.parameters(3), doc_a.abs.fit_dog.parameters(3), tolerance.abs.fit_dog.parameters(3), 'abs fit_dog parameter 3');
[b_(end+1),errormsg_{end+1}] = ndi.test.values_within_tolerance(doc_e.abs.fit_dog.parameters(4), doc_a.abs.fit_dog.parameters(4), tolerance.abs.fit_dog.parameters(4), 'abs fit_dog parameter 4');
[b_(64),errormsg_{64}] = ndi.test.values_within_tolerance(doc_e.abs.fit_dog.values, doc_a.abs.fit_dog.values, tolerance.abs.fit_dog.values, 'abs fit_dog values');
[b_(65),errormsg_{65}] = ndi.test.values_within_tolerance(doc_e.abs.fit_dog.fit, doc_a.abs.fit_dog.fit, tolerance.abs.fit_dog.fit, 'abs fit_dog fit');
[b_(66),errormsg_{66}] = ndi.test.values_within_tolerance(doc_e.abs.fit_dog.L50, doc_a.abs.fit_dog.L50, tolerance.abs.fit_dog.L50, 'abs fit_dog L50');
[b_(67),errormsg_{67}] = ndi.test.values_within_tolerance(doc_e.abs.fit_dog.Pref, doc_a.abs.fit_dog.Pref, tolerance.abs.fit_dog.Pref, 'abs fit_dog Pref');
[b_(68),errormsg_{68}] = ndi.test.values_within_tolerance(doc_e.abs.fit_dog.H50, doc_a.abs.fit_dog.H50, tolerance.abs.fit_dog.H50, 'abs fit_dog H50');
[b_(69),errormsg_{69}] = ndi.test.values_within_tolerance(doc_e.abs.fit_dog.bandwidth, doc_a.abs.fit_dog.bandwidth, tolerance.abs.fit_dog.bandwidth, 'abs fit_dog bandwidth');

%    Comparing fit_movshon
%        parameters
%        values
%        fit
%        L50
%        Pref
%        H50
%        R2
%        bandwidth

[b_(70),errormsg_{70}] = ndi.test.values_within_tolerance(doc_e.abs.fit_movshon.parameters(1), doc_a.abs.fit_movshon.parameters(1), tolerance.abs.fit_movshon.parameters(1), 'abs fit_movshon parameter 1');
[b_(end+1),errormsg_{end+1}] = ndi.test.values_within_tolerance(doc_e.abs.fit_movshon.parameters(2), doc_a.abs.fit_movshon.parameters(2), tolerance.abs.fit_movshon.parameters(2), 'abs fit_movshon parameter 2');
[b_(end+1),errormsg_{end+1}] = ndi.test.values_within_tolerance(doc_e.abs.fit_movshon.parameters(3), doc_a.abs.fit_movshon.parameters(3), tolerance.abs.fit_movshon.parameters(3), 'abs fit_movshon parameter 3');
[b_(end+1),errormsg_{end+1}] = ndi.test.values_within_tolerance(doc_e.abs.fit_movshon.parameters(4), doc_a.abs.fit_movshon.parameters(4), tolerance.abs.fit_movshon.parameters(4), 'abs fit_movshon parameter 4');
[b_(71),errormsg_{71}] = ndi.test.values_within_tolerance(doc_e.abs.fit_movshon.values, doc_a.abs.fit_movshon.values, tolerance.abs.fit_movshon.values, 'abs fit_movshon values');
[b_(72),errormsg_{72}] = ndi.test.values_within_tolerance(doc_e.abs.fit_movshon.fit, doc_a.abs.fit_movshon.fit, tolerance.abs.fit_movshon.fit, 'abs fit_movshon fit');
[b_(73),errormsg_{73}] = ndi.test.values_within_tolerance(doc_e.abs.fit_movshon.L50, doc_a.abs.fit_movshon.L50, tolerance.abs.fit_movshon.L50, 'abs fit_movshon L50');
[b_(74),errormsg_{74}] = ndi.test.values_within_tolerance(doc_e.abs.fit_movshon.Pref, doc_a.abs.fit_movshon.Pref, tolerance.abs.fit_movshon.Pref, 'abs fit_movshon Pref');
[b_(75),errormsg_{75}] = ndi.test.values_within_tolerance(doc_e.abs.fit_movshon.H50, doc_a.abs.fit_movshon.H50, tolerance.abs.fit_movshon.H50, 'abs fit_movshon H50');
[b_(76),errormsg_{76}] = ndi.test.values_within_tolerance(doc_e.abs.fit_movshon.R2, doc_a.abs.fit_movshon.R2, tolerance.abs.fit_movshon.R2, 'abs fit_movshon R2');
[b_(77),errormsg_{77}] = ndi.test.values_within_tolerance(doc_e.abs.fit_movshon.bandwidth, doc_a.abs.fit_movshon.bandwidth, tolerance.abs.fit_movshon.bandwidth, 'abs fit_movshon bandwidth');

%    Comparing fit_movshon_c
%        parameters
%        values
%        fit
%        L50
%        Pref
%        H50
%        R2
%        bandwidth

[b_(78),errormsg_{78}] = ndi.test.values_within_tolerance(doc_e.abs.fit_movshon_c.parameters(1), doc_a.abs.fit_movshon_c.parameters(1), tolerance.abs.fit_movshon_c.parameters(1), 'abs fit_movshon_c parameter 1');
[b_(end+1),errormsg_{end+1}] = ndi.test.values_within_tolerance(doc_e.abs.fit_movshon_c.parameters(2), doc_a.abs.fit_movshon_c.parameters(2), tolerance.abs.fit_movshon_c.parameters(2), 'abs fit_movshon_c parameter 2');
[b_(end+1),errormsg_{end+1}] = ndi.test.values_within_tolerance(doc_e.abs.fit_movshon_c.parameters(3), doc_a.abs.fit_movshon_c.parameters(3), tolerance.abs.fit_movshon_c.parameters(3), 'abs fit_movshon_c parameter 3');
[b_(end+1),errormsg_{end+1}] = ndi.test.values_within_tolerance(doc_e.abs.fit_movshon_c.parameters(4), doc_a.abs.fit_movshon_c.parameters(4), tolerance.abs.fit_movshon_c.parameters(4), 'abs fit_movshon_c parameter 4');
[b_(end+1),errormsg_{end+1}] = ndi.test.values_within_tolerance(doc_e.abs.fit_movshon_c.parameters(5), doc_a.abs.fit_movshon_c.parameters(5), tolerance.abs.fit_movshon_c.parameters(5), 'abs fit_movshon_c parameter 5');
[b_(79),errormsg_{79}] = ndi.test.values_within_tolerance(doc_e.abs.fit_movshon_c.values, doc_a.abs.fit_movshon_c.values, tolerance.abs.fit_movshon_c.values, 'abs fit_movshon_c values');
[b_(80),errormsg_{80}] = ndi.test.values_within_tolerance(doc_e.abs.fit_movshon_c.fit, doc_a.abs.fit_movshon_c.fit, tolerance.abs.fit_movshon_c.fit, 'abs fit_movshon_c fit');
[b_(81),errormsg_{81}] = ndi.test.values_within_tolerance(doc_e.abs.fit_movshon_c.L50, doc_a.abs.fit_movshon_c.L50, tolerance.abs.fit_movshon_c.L50, 'abs fit_movshon_c L50');
[b_(82),errormsg_{82}] = ndi.test.values_within_tolerance(doc_e.abs.fit_movshon_c.Pref, doc_a.abs.fit_movshon_c.Pref, tolerance.abs.fit_movshon_c.Pref, 'abs fit_movshon_c Pref');
[b_(83),errormsg_{83}] = ndi.test.values_within_tolerance(doc_e.abs.fit_movshon_c.H50, doc_a.abs.fit_movshon_c.H50, tolerance.abs.fit_movshon_c.H50, 'abs fit_movshon_c H50');
[b_(84),errormsg_{84}] = ndi.test.values_within_tolerance(doc_e.abs.fit_movshon_c.R2, doc_a.abs.fit_movshon_c.R2, tolerance.abs.fit_movshon_c.R2, 'abs fit_movshon_c R2');
[b_(85),errormsg_{85}] = ndi.test.values_within_tolerance(doc_e.abs.fit_movshon_c.bandwidth, doc_a.abs.fit_movshon_c.bandwidth, tolerance.abs.fit_movshon_c.bandwidth, 'abs fit_movshon_c bandwidth');

%    Comparing fit_spline
%        values
%        fit
%        L50
%        Pref
%        H50
%        bandwidth

[b_(86),errormsg_{86}] = ndi.test.values_within_tolerance(doc_e.abs.fit_spline.values, doc_a.abs.fit_spline.values, tolerance.abs.fit_spline.values, 'abs fit_spline values');
[b_(87),errormsg_{87}] = ndi.test.values_within_tolerance(doc_e.abs.fit_spline.fit, doc_a.abs.fit_spline.fit, tolerance.abs.fit_spline.fit, 'abs fit_spline fit');
[b_(88),errormsg_{88}] = ndi.test.values_within_tolerance(doc_e.abs.fit_spline.L50, doc_a.abs.fit_spline.L50, tolerance.abs.fit_spline.L50, 'abs fit_spline L50');
[b_(89),errormsg_{89}] = ndi.test.values_within_tolerance(doc_e.abs.fit_spline.Pref, doc_a.abs.fit_spline.Pref, tolerance.abs.fit_spline.Pref, 'abs fit_spline Pref');
[b_(90),errormsg_{90}] = ndi.test.values_within_tolerance(doc_e.abs.fit_spline.H50, doc_a.abs.fit_spline.H50, tolerance.abs.fit_spline.H50, 'abs fit_spline H50');
[b_(91),errormsg_{91}] = ndi.test.values_within_tolerance(doc_e.abs.fit_spline.bandwidth, doc_a.abs.fit_spline.bandwidth, tolerance.abs.fit_spline.bandwidth, 'abs fit_spline bandwidth');

%    Comparing fit_gausslog
%        parameters
%        values
%        fit
%        L50
%        Pref
%        H50
%        bandwidth

[b_(92),errormsg_{92}] = ndi.test.values_within_tolerance(doc_e.abs.fit_gausslog.parameters(1), doc_a.abs.fit_gausslog.parameters(1), tolerance.abs.fit_gausslog.parameters(1), 'abs fit_gausslog parameter 1');
[b_(end+1),errormsg_{end+1}] = ndi.test.values_within_tolerance(doc_e.abs.fit_gausslog.parameters(2), doc_a.abs.fit_gausslog.parameters(2), tolerance.abs.fit_gausslog.parameters(2), 'abs fit_gausslog parameter 2');
[b_(end+1),errormsg_{end+1}] = ndi.test.values_within_tolerance(doc_e.abs.fit_gausslog.parameters(3), doc_a.abs.fit_gausslog.parameters(3), tolerance.abs.fit_gausslog.parameters(3), 'abs fit_gausslog parameter 3');
[b_(end+1),errormsg_{end+1}] = ndi.test.values_within_tolerance(doc_e.abs.fit_gausslog.parameters(4), doc_a.abs.fit_gausslog.parameters(4), tolerance.abs.fit_gausslog.parameters(4), 'abs fit_gausslog parameter 4');
[b_(93),errormsg_{93}] = ndi.test.values_within_tolerance(doc_e.abs.fit_gausslog.values, doc_a.abs.fit_gausslog.values, tolerance.abs.fit_gausslog.values, 'abs fit_gausslog values');
[b_(94),errormsg_{94}] = ndi.test.values_within_tolerance(doc_e.abs.fit_gausslog.fit, doc_a.abs.fit_gausslog.fit, tolerance.abs.fit_gausslog.fit, 'abs fit_gausslog fit');
[b_(95),errormsg_{95}] = ndi.test.values_within_tolerance(doc_e.abs.fit_gausslog.L50, doc_a.abs.fit_gausslog.L50, tolerance.abs.fit_gausslog.L50, 'abs fit_gausslog L50');
[b_(96),errormsg_{96}] = ndi.test.values_within_tolerance(doc_e.abs.fit_gausslog.Pref, doc_a.abs.fit_gausslog.Pref, tolerance.abs.fit_gausslog.Pref, 'abs fit_gausslog Pref');
[b_(97),errormsg_{97}] = ndi.test.values_within_tolerance(doc_e.abs.fit_gausslog.H50, doc_a.abs.fit_gausslog.H50, tolerance.abs.fit_gausslog.H50, 'abs fit_gausslog H50');
[b_(98),errormsg_{98}] = ndi.test.values_within_tolerance(doc_e.abs.fit_gausslog.bandwidth, doc_a.abs.fit_gausslog.bandwidth, tolerance.abs.fit_gausslog.bandwidth, 'abs fit_gausslog bandwidth');


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
