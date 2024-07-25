function [b_, errormsg_] = spatial_frequency_tuning_compare_docs(document_expected, document_actual, scope)
% Comparing spatial_frequency_tuning
%
% [B, ERRORMSG] = ndi.calc.vis.test.spatial_frequency_tuning_compare_docs(DOC_EXPECTED, DOC_ACTUAL, SCOPE)
%
%

% Initialize b_ as a row vector of ones for tracking comparison results
% Initialize errormsg_ as an empty cell array to hold error messages

b_ = ones(1,100);
errormsg_ = cell(1,100);

 % establish scope-dependent tolerances
switch(scope),
    case 'standard',
    
       tol_tuning_curve.spatial_frequency = 1e-6;
       tol_tuning_curve.mean = 1e-6;
       tol_tuning_curve.stddev = 1e-6;
       tol_tuning_curve.stderr = 1e-6;
       tol_tuning_curve.individual = .5; %responses in units of Hz
       tol_tuning_curve.control_mean = 1e-6;
       tol_tuning_curve.control_stddev = 1e-6;
       tol_tuning_curve.control_stderr = 1e-6;
       tol_tuning_curve.control_mean_stddev = 1e-6;
       tol_tuning_curve.control_mean_stderr = 1e-6;
       tol_significance.visual_response_anova_p = .1;
       tol_significance.across_stimuli_anova_p = .1;
       tol_fitless.H50 = 1; %spatial freq range: 10^-2 to 60
       tol_fitless.Pref = 1;
       tol_fitless.L50 = 1;
       tol_fitless.bandwidth = 1;
       tol_fitless.low_pass_index = 1; % 1 to 100 technically, but in practice should be 1
       tol_fitless.high_pass_index = 1; % 1 to 100 technically, but in practice is 100
       tol_fit_spline.fit = 1;
       tol_fit_spline.H50 = 1;
       tol_fit_spline.Pref = 1;
       tol_fit_spline.L50 = 1;
       tol_fit_spline.values = 1;
       tol_fit_spline.bandwidth = 1;
       tol_fit_dog.fit = 2; %double the tolerance for fits because fits can look different even for the same underlying data
       tol_fit_dog.H50 = 2;
       tol_fit_dog.Pref = 2;
       tol_fit_dog.L50 = 2;
       tol_fit_dog.values = 2;
       tol_fit_dog.bandwidth = 2;
       tol_fit_dog.parameters(1) = 2;
       tol_fit_dog.parameters(2) = 2;
       tol_fit_dog.parameters(3) = 2;
       tol_fit_dog.parameters(4) = 2;
       tol_fit_gausslog.fit = 2;
       tol_fit_gausslog.H50 = 2;
       tol_fit_gausslog.Pref = 2;
       tol_fit_gausslog.L50 = 2;
       tol_fit_gausslog.values = 2;
       tol_fit_gausslog.bandwidth = 2;
       tol_fit_gausslog.parameters(1) = 1; %offset: about -5 to 5
       tol_fit_gausslog.parameters(2) = 2; %height above offset: about 10 to 100
       tol_fit_gausslog.parameters(3) = 2; %peak location: 0 to 60, but shouldn't vary much
       tol_fit_gausslog.parameters(4) = 2; %width: 0 to 60
       %tol_fit_gausslog.parameters(5) = 0.1; %not included in output
       %document for some reason
       tol_fit_movshon.fit = 2;
       tol_fit_movshon.H50 = 2;
       tol_fit_movshon.Pref = 2;
       tol_fit_movshon.L50 = 2;
       tol_fit_movshon.values = 2;
       tol_fit_movshon.parameters(1) = 2; %scaling factor
       tol_fit_movshon.parameters(2) = 2; %characteristic spatial frequency
       tol_fit_movshon.parameters(3) = 2; %corner frequency of low-frequency limb
       tol_fit_movshon.parameters(4) = 2; %slope of low-frequency limb
       tol_fit_movshon.R2 = .1;
       tol_fit_movshon.bandwidth = .6; % about 0 to 6
       
       tol_fit_movshon_c.fit = 2;
       tol_fit_movshon_c.H50 = 2;
       tol_fit_movshon_c.Pref = 2;
       tol_fit_movshon_c.L50 = 2;
       tol_fit_movshon_c.values = 2;
       tol_fit_movshon_c.parameters(1) = 2;
       tol_fit_movshon_c.parameters(2) = 2;
       tol_fit_movshon_c.parameters(3) = 2;
       tol_fit_movshon_c.parameters(4) = 2;
       tol_fit_movshon_c.parameters(5) = 1;
       tol_fit_movshon_c.R2 = 0.1;
       tol_fit_movshon_c.bandwidth = .6;
       
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


    case 'low_noise',
    
       tol_tuning_curve.spatial_frequency = 1e-6;
       tol_tuning_curve.mean = 3;
       tol_tuning_curve.stddev = 3;
       tol_tuning_curve.stderr = 1;
       tol_tuning_curve.individual = 9; %responses in units of Hz
       tol_tuning_curve.control_mean = 3;
       tol_tuning_curve.control_stddev = 3;
       tol_tuning_curve.control_stderr = 1;
       tol_tuning_curve.control_mean_stddev = 1;
       tol_tuning_curve.control_mean_stderr = .3;
       tol_significance.visual_response_anova_p = .1;
       tol_significance.across_stimuli_anova_p = .1;
       tol_fitless.H50 = 2; %spatial freq range: 10^-2 to 60
       tol_fitless.Pref = 2;
       tol_fitless.L50 = 2;
       tol_fitless.bandwidth = 2;
       tol_fitless.low_pass_index = 1;
       tol_fitless.high_pass_index = 1;
       tol_fit_spline.fit = 2;
       tol_fit_spline.H50 = 2;
       tol_fit_spline.Pref = 2;
       tol_fit_spline.L50 = 2;
       tol_fit_spline.values = 2;
       tol_fit_spline.bandwidth = 2;
       tol_fit_dog.fit = 4; %double the tolerance for fits because fits can look different even for the same underlying data
       tol_fit_dog.H50 = 4;
       tol_fit_dog.Pref = 4;
       tol_fit_dog.L50 = 4;
       tol_fit_dog.values = 4;
       tol_fit_dog.bandwidth = 4;
       tol_fit_dog.parameters(1) = 4;
       tol_fit_dog.parameters(2) = 4;
       tol_fit_dog.parameters(3) = 4;
       tol_fit_dog.parameters(4) = 4;
       tol_fit_gausslog.fit = 4;
       tol_fit_gausslog.H50 = 4;
       tol_fit_gausslog.Pref = 4;
       tol_fit_gausslog.L50 = 4;
       tol_fit_gausslog.values = 4;
       tol_fit_gausslog.bandwidth = 4;
       tol_fit_gausslog.parameters(1) = 2; %offset: about -5 to 5
       tol_fit_gausslog.parameters(2) = 4; %height above offset: about 10 to 100
       tol_fit_gausslog.parameters(3) = 4; %peak location: 0 to 60, but shouldn't vary much
       tol_fit_gausslog.parameters(4) = 4; %width: 0 to 60
       %tol_fit_gausslog.parameters(5) = 0.1; %not included in output
       %document for some reason
       tol_fit_movshon.fit = 4;
       tol_fit_movshon.H50 = 4;
       tol_fit_movshon.Pref = 4;
       tol_fit_movshon.L50 = 4;
       tol_fit_movshon.values = 4;
       tol_fit_movshon.parameters(1) = 4; %scaling factor
       tol_fit_movshon.parameters(2) = 4; %characteristic spatial frequency
       tol_fit_movshon.parameters(3) = 4; %corner frequency of low-frequency limb
       tol_fit_movshon.parameters(4) = 4; %slope of low-frequency limb
       tol_fit_movshon.R2 = .2;
       tol_fit_movshon.bandwidth = 1.2; % about 0 to 6
       
       tol_fit_movshon_c.fit = 4;
       tol_fit_movshon_c.H50 = 4;
       tol_fit_movshon_c.Pref = 4;
       tol_fit_movshon_c.L50 = 4;
       tol_fit_movshon_c.values = 4;
       tol_fit_movshon_c.parameters(1) = 4;
       tol_fit_movshon_c.parameters(2) = 4;
       tol_fit_movshon_c.parameters(3) = 4;
       tol_fit_movshon_c.parameters(4) = 4;
       tol_fit_movshon_c.parameters(5) = 3;
       tol_fit_movshon_c.R2 = 0.2;
       tol_fit_movshon_c.bandwidth = 1.2;

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
    
       tol_tuning_curve.spatial_frequency = 1e-6;
       tol_tuning_curve.mean = 6;
       tol_tuning_curve.stddev = 6;
       tol_tuning_curve.stderr = 2;
       tol_tuning_curve.individual = 18; %responses in units of Hz
       tol_tuning_curve.control_mean = 3;
       tol_tuning_curve.control_stddev = 3;
       tol_tuning_curve.control_stderr = 1;
       tol_tuning_curve.control_mean_stddev = 1;
       tol_tuning_curve.control_mean_stderr = .3;
       tol_significance.visual_response_anova_p = .1;
       tol_significance.across_stimuli_anova_p = .1;
       tol_fitless.H50 = 4; %spatial freq range: 10^-2 to 60
       tol_fitless.Pref = 4;
       tol_fitless.L50 = 4;
       tol_fitless.bandwidth = 4;
       tol_fitless.low_pass_index = 4;
       tol_fitless.high_pass_index = 4;
       tol_fit_spline.fit = 4;
       tol_fit_spline.H50 = 4;
       tol_fit_spline.Pref = 4;
       tol_fit_spline.L50 = 4;
       tol_fit_spline.values = 4;
       tol_fit_spline.bandwidth = 4;
       tol_fit_dog.fit = 8; %double the tolerance for fits because fits can look different even for the same underlying data
       tol_fit_dog.H50 = 8;
       tol_fit_dog.Pref = 8;
       tol_fit_dog.L50 = 8;
       tol_fit_dog.values = 8;
       tol_fit_dog.bandwidth = 8;
       tol_fit_dog.parameters(1) = 8;
       tol_fit_dog.parameters(2) = 8;
       tol_fit_dog.parameters(3) = 8;
       tol_fit_dog.parameters(4) = 8;
       tol_fit_gausslog.fit = 8;
       tol_fit_gausslog.H50 = 8;
       tol_fit_gausslog.Pref = 8;
       tol_fit_gausslog.L50 = 8;
       tol_fit_gausslog.values = 8;
       tol_fit_gausslog.bandwidth = 8;
       tol_fit_gausslog.parameters(1) = 4; %offset: about -5 to 5
       tol_fit_gausslog.parameters(2) = 8; %height above offset: about 10 to 100
       tol_fit_gausslog.parameters(3) = 8; %peak location: 0 to 60, but shouldn't vary much
       tol_fit_gausslog.parameters(4) = 8; %width: 0 to 60
       %tol_fit_gausslog.parameters(5) = 0.1; %not included in output
       %document for some reason
       tol_fit_movshon.fit = 8;
       tol_fit_movshon.H50 = 8;
       tol_fit_movshon.Pref = 8;
       tol_fit_movshon.L50 = 8;
       tol_fit_movshon.values = 8;
       tol_fit_movshon.parameters(1) = 8; %scaling factor
       tol_fit_movshon.parameters(2) = 8; %characteristic spatial frequency
       tol_fit_movshon.parameters(3) = 8; %corner frequency of low-frequency limb
       tol_fit_movshon.parameters(4) = 8; %slope of low-frequency limb
       tol_fit_movshon.R2 = .4;
       tol_fit_movshon.bandwidth = 2.4; % about 0 to 6
       
       tol_fit_movshon_c.fit = 8;
       tol_fit_movshon_c.H50 = 8;
       tol_fit_movshon_c.Pref = 8;
       tol_fit_movshon_c.L50 = 8;
       tol_fit_movshon_c.values = 8;
       tol_fit_movshon_c.parameters(1) = 8;
       tol_fit_movshon_c.parameters(2) = 8;
       tol_fit_movshon_c.parameters(3) = 8;
       tol_fit_movshon_c.parameters(4) = 8;
       tol_fit_movshon_c.parameters(5) = 4;
       tol_fit_movshon_c.R2 = 0.4;
       tol_fit_movshon_c.bandwidth = 2.4;

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

doc_e = document_expected.document_properties.spatial_frequency_tuning;
doc_a = document_actual.document_properties.spatial_frequency_tuning;

% Comparing properties
%   Response Units

properties_match = strcmpi(char(doc_e.properties.response_units), char(doc_a.properties.response_units));
if ~properties_match
   b_(80) = 0;
   errormsg_{80} = ['Expected response units in ' doc_e.properties.response_units ' but observed ' doc_a.properties.response_units];
   return;
end

%   Response Type

properties_match = strcmpi(char(doc_e.properties.response_type), char(doc_a.properties.response_type));
if ~properties_match
   b_(81) = 0;
   errormsg_{81} = ['Expected response type of ' doc_e.properties.response_type ' but observed ' doc_a.properties.response_type];
   return;
end

% Comparing Tuning_curve
%	spatial frequency                               
%	mean                                   
%	stddev                                 
%	stderr                                
%	individual
%   control_mean
%	control_stddev                        
%	control_stderr
%   control_mean_stddev
%   control_mean_stderr

[b_(1),errormsg_{1}] = ndi.test.values_within_tolerance(doc_e.tuning_curve.spatial_frequency, doc_a.tuning_curve.spatial_frequency, tol_tuning_curve.spatial_frequency, 'spatial_frequency');
[b_(2),errormsg_{2}] = ndi.test.values_within_tolerance(doc_e.tuning_curve.mean, doc_a.tuning_curve.mean, tol_tuning_curve.mean, 'mean');
[b_(3),errormsg_{3}] = ndi.test.values_within_tolerance(doc_e.tuning_curve.stddev, doc_a.tuning_curve.stddev, tol_tuning_curve.stddev, 'stddev');
[b_(4),errormsg_{4}] = ndi.test.values_within_tolerance(doc_e.tuning_curve.stderr, doc_a.tuning_curve.stderr, tol_tuning_curve.stderr, 'stderr');
%[b_(5),errormsg_{5}] = ndi.test.values_within_tolerance(doc_e.tuning_curve.individual, doc_a.tuning_curve.individual, tol_tuning_curve.individual, 'individual');
%compare each rep in doc_a with the first rep of doc_e
for i = 1:size(doc_a.tuning_curve.individual,1)
    [b_(end+1),errormsg_{end+1}] = ndi.test.values_within_tolerance(doc_e.tuning_curve.individual(1,:), doc_a.tuning_curve.individual(i,:), tol_tuning_curve.individual, ['individual ',num2str(i)]);
end
[b_(end+1),errormsg_{end+1}] = ndi.test.values_within_tolerance(doc_e.tuning_curve.control_mean, doc_a.tuning_curve.control_mean, tol_tuning_curve.control_mean, 'control_mean');
[b_(6),errormsg_{6}] = ndi.test.values_within_tolerance(doc_e.tuning_curve.control_stddev, doc_a.tuning_curve.control_stddev, tol_tuning_curve.control_stddev, 'control_stddev');
[b_(7),errormsg_{7}] = ndi.test.values_within_tolerance(doc_e.tuning_curve.control_stderr, doc_a.tuning_curve.control_stderr, tol_tuning_curve.control_stderr, 'control_stderr');
[b_(end+1),errormsg_{end+1}] = ndi.test.values_within_tolerance(doc_e.tuning_curve.control_mean_stddev, doc_a.tuning_curve.control_mean_stddev, tol_tuning_curve.control_mean_stddev, 'control_mean_stddev');
[b_(end+1),errormsg_{end+1}] = ndi.test.values_within_tolerance(doc_e.tuning_curve.control_mean_stderr, doc_a.tuning_curve.control_mean_stderr, tol_tuning_curve.control_mean_stderr, 'control_mean_stderr');

% Comparing Significance
%   visual_response_anova_p
%   across_stimuli_anova_p

[b_(8),errormsg_{8}] = ndi.test.values_within_tolerance(doc_e.significance.visual_response_anova_p, doc_a.significance.visual_response_anova_p, tol_significance.visual_response_anova_p, 'visual response anova p');
[b_(9),errormsg_{9}] = ndi.test.values_within_tolerance(doc_e.significance.across_stimuli_anova_p, doc_a.significance.across_stimuli_anova_p, tol_significance.across_stimuli_anova_p, 'across stimuli anova p');

% Comparing Fitless
%   H50
%   Pref
%   L50
%   bandwidth
%   low_pass_index
%   high_pass_index

[b_(10),errormsg_{10}] = ndi.test.values_within_tolerance(doc_e.fitless.H50, doc_a.fitless.H50, tol_fitless.H50, 'fitless H50');
[b_(11),errormsg_{11}] = ndi.test.values_within_tolerance(doc_e.fitless.Pref, doc_a.fitless.Pref, tol_fitless.Pref, 'fitless Pref');
[b_(12),errormsg_{12}] = ndi.test.values_within_tolerance(doc_e.fitless.L50, doc_a.fitless.L50, tol_fitless.L50, 'fitless L50');

% Comparing fit_spline
%   fit
%   H50
%   Pref
%   L50
%   values
%   bandwidth

[b_(13),errormsg_{13}] = ndi.test.values_within_tolerance(doc_e.fit_spline.fit, doc_a.fit_spline.fit, tol_fit_spline.fit, 'fit_spline fit');
[b_(14),errormsg_{14}] = ndi.test.values_within_tolerance(doc_e.fit_spline.H50, doc_a.fit_spline.H50, tol_fit_spline.H50, 'fit_spline H50');
[b_(15),errormsg_{15}] = ndi.test.values_within_tolerance(doc_e.fit_spline.Pref, doc_a.fit_spline.Pref, tol_fit_spline.Pref, 'fit_spline Pref');
[b_(16),errormsg_{16}] = ndi.test.values_within_tolerance(doc_e.fit_spline.L50, doc_a.fit_spline.L50, tol_fit_spline.L50, 'fit_spline L50');
[b_(17),errormsg_{17}] = ndi.test.values_within_tolerance(doc_e.fit_spline.values, doc_a.fit_spline.values, tol_fit_spline.values, 'fit_spline values');

% Comparing fit_dog
%   fit
%   H50
%   Pref
%   L50
%   values
%   bandwidth
%   parameters

[b_(18),errormsg_{18}] = ndi.test.values_within_tolerance(doc_e.fit_dog.fit, doc_a.fit_dog.fit, tol_fit_dog.fit, 'fit_dog fit');
[b_(19),errormsg_{19}] = ndi.test.values_within_tolerance(doc_e.fit_dog.H50, doc_a.fit_dog.H50, tol_fit_dog.H50, 'fit_dog H50');
[b_(20),errormsg_{20}] = ndi.test.values_within_tolerance(doc_e.fit_dog.Pref, doc_a.fit_dog.Pref, tol_fit_dog.Pref, 'fit_dog Pref');
[b_(21),errormsg_{21}] = ndi.test.values_within_tolerance(doc_e.fit_dog.L50, doc_a.fit_dog.L50, tol_fit_dog.L50, 'fit_dog L50');
[b_(22),errormsg_{22}] = ndi.test.values_within_tolerance(doc_e.fit_dog.values, doc_a.fit_dog.values, tol_fit_dog.values, 'fit_dog values');
[b_(23),errormsg_{23] = ndi.test.values_within_tolerance(doc_e.fit_dog.parameters(1), doc_a.fit_dog.parameters(1), tol_fit_dog.parameters(1), 'fit_dog parameter 1');
[b_(75),errormsg_{75}] = ndi.test.values_within_tolerance(doc_e.fit_dog.parameters(2), doc_a.fit_dog.parameters(2), tol_fit_dog.parameters(2), 'fit_dog parameter 2');
[b_(76),errormsg_{76}] = ndi.test.values_within_tolerance(doc_e.fit_dog.parameters(3), doc_a.fit_dog.parameters(3), tol_fit_dog.parameters(3), 'fit_dog parameter 3');
[b_(77),errormsg_{77}] = ndi.test.values_within_tolerance(doc_e.fit_dog.parameters(4), doc_a.fit_dog.parameters(4), tol_fit_dog.parameters(4), 'fit_dog parameter 4');


% Comparing fit_gausslog
%   fit
%   H50
%   Pref
%   L50
%   values
%   bandwidth
%   parameters

[b_(24),errormsg_{24}] = ndi.test.values_within_tolerance(doc_e.fit_dog.fit, doc_a.fit_dog.fit, tol_fit_gausslog.fit, 'fit_gausslog fit');
[b_(25),errormsg_{25}] = ndi.test.values_within_tolerance(doc_e.fit_dog.H50, doc_a.fit_dog.H50, tol_fit_gausslog.H50, 'fit_gausslog H50');
[b_(26),errormsg_{26}] = ndi.test.values_within_tolerance(doc_e.fit_dog.Pref, doc_a.fit_dog.Pref, tol_fit_gausslog.Pref, 'fit_gausslog Pref');
[b_(27),errormsg_{27}] = ndi.test.values_within_tolerance(doc_e.fit_dog.L50, doc_a.fit_dog.L50, tol_fit_gausslog.L50, 'fit_gausslog L50');
[b_(28),errormsg_{28}] = ndi.test.values_within_tolerance(doc_e.fit_dog.values, doc_a.fit_dog.values, tol_fit_gausslog.values, 'fit_gausslog values');
[b_(29),errormsg_{29}] = ndi.test.values_within_tolerance(doc_e.fit_dog.parameters(1), doc_a.fit_dog.parameters(1), tol_fit_gausslog.parameters(1), 'fit_gausslog parameter 1');
[b_(30),errormsg_{30}] = ndi.test.values_within_tolerance(doc_e.fit_dog.parameters(2), doc_a.fit_dog.parameters(2), tol_fit_gausslog.parameters(2), 'fit_gausslog parameter 2');
[b_(31),errormsg_{31}] = ndi.test.values_within_tolerance(doc_e.fit_dog.parameters(3), doc_a.fit_dog.parameters(3), tol_fit_gausslog.parameters(3), 'fit_gausslog parameter 3');
[b_(32),errormsg_{32}] = ndi.test.values_within_tolerance(doc_e.fit_dog.parameters(4), doc_a.fit_dog.parameters(4), tol_fit_gausslog.parameters(4), 'fit_gausslog parameter 4');

% abs values
%  Comparing fitless
%    L50
%    Pref
%    H50
%    bandwidth
%    low_pass_index
%    high_pass_index

[b_(33),errormsg_{33}] = ndi.test.values_within_tolerance(doc_e.abs.fitless.L50, doc_a.abs.fitless.L50, tolerance.abs.fitless.L50, 'abs fitless L50');
[b_(34),errormsg_{34}] = ndi.test.values_within_tolerance(doc_e.abs.fitless.Pref, doc_a.abs.fitless.Pref, tolerance.abs.fitless.Pref, 'abs fitless Pref');
[b_(35),errormsg_{35}] = ndi.test.values_within_tolerance(doc_e.abs.fitless.H50, doc_a.abs.fitless.H50, tolerance.abs.fitless.H50, 'abs fitless H50');
[b_(36),errormsg_{36}] = ndi.test.values_within_tolerance(doc_e.abs.fitless.bandwidth, doc_a.abs.fitless.bandwidth, tolerance.abs.fitless.bandwidth, 'abs fitless bandwidth');
[b_(37),errormsg_{37}] = ndi.test.values_within_tolerance(doc_e.abs.fitless.low_pass_index, doc_a.abs.fitless.low_pass_index, tolerance.abs.fitless.low_pass_index, 'abs fitless low_pass_index');
[b_(38),errormsg_{38}] = ndi.test.values_within_tolerance(doc_e.abs.fitless.high_pass_index, doc_a.abs.fitless.high_pass_index, tolerance.abs.fitless.high_pass_index, 'abs fitless high_pass_index');

%  Comparing fit_dog
%    parameters
%    values
%    fit
%    L50
%    Pref
%    H50
%    bandwidth

[b_(39),errormsg_{39}] = ndi.test.values_within_tolerance(doc_e.abs.fit_dog.parameters, doc_a.abs.fit_dog.parameters, tolerance.abs.fit_dog.parameters, 'abs fit_dog parameters');
[b_(40),errormsg_{40}] = ndi.test.values_within_tolerance(doc_e.abs.fit_dog.values, doc_a.abs.fit_dog.values, tolerance.abs.fit_dog.values, 'abs fit_dog values');
[b_(41),errormsg_{41}] = ndi.test.values_within_tolerance(doc_e.abs.fit_dog.fit, doc_a.abs.fit_dog.fit, tolerance.abs.fit_dog.fit, 'abs fit_dog fit');
[b_(42),errormsg_{42}] = ndi.test.values_within_tolerance(doc_e.abs.fit_dog.L50, doc_a.abs.fit_dog.L50, tolerance.abs.fit_dog.L50, 'abs fit_dog L50');
[b_(43),errormsg_{43}] = ndi.test.values_within_tolerance(doc_e.abs.fit_dog.Pref, doc_a.abs.fit_dog.Pref, tolerance.abs.fit_dog.Pref, 'abs fit_dog Pref');
[b_(44),errormsg_{44}] = ndi.test.values_within_tolerance(doc_e.abs.fit_dog.H50, doc_a.abs.fit_dog.H50, tolerance.abs.fit_dog.H50, 'abs fit_dog H50');
[b_(45),errormsg_{45}] = ndi.test.values_within_tolerance(doc_e.abs.fit_dog.bandwidth, doc_a.abs.fit_dog.bandwidth, tolerance.abs.fit_dog.bandwidth, 'abs fit_dog bandwidth');

%  Comparing fit_movshon
%    parameters
%    values
%    fit
%    L50
%    Pref
%    H50
%    R2
%    bandwidth

[b_(46),errormsg_{46}] = ndi.test.values_within_tolerance(doc_e.abs.fit_movshon.parameters, doc_a.abs.fit_movshon.parameters, tolerance.abs.fit_movshon.parameters, 'abs fit_movshon parameters');
[b_(47),errormsg_{47}] = ndi.test.values_within_tolerance(doc_e.abs.fit_movshon.values, doc_a.abs.fit_movshon.values, tolerance.abs.fit_movshon.values, 'abs fit_movshon values');
[b_(48),errormsg_{48}] = ndi.test.values_within_tolerance(doc_e.abs.fit_movshon.fit, doc_a.abs.fit_movshon.fit, tolerance.abs.fit_movshon.fit, 'abs fit_movshon fit');
[b_(49),errormsg_{49}] = ndi.test.values_within_tolerance(doc_e.abs.fit_movshon.L50, doc_a.abs.fit_movshon.L50, tolerance.abs.fit_movshon.L50, 'abs fit_movshon L50');
[b_(50),errormsg_{50}] = ndi.test.values_within_tolerance(doc_e.abs.fit_movshon.Pref, doc_a.abs.fit_movshon.Pref, tolerance.abs.fit_movshon.Pref, 'abs fit_movshon Pref');
[b_(51),errormsg_{51}] = ndi.test.values_within_tolerance(doc_e.abs.fit_movshon.H50, doc_a.abs.fit_movshon.H50, tolerance.abs.fit_movshon.H50, 'abs fit_movshon H50');
[b_(52),errormsg_{52}] = ndi.test.values_within_tolerance(doc_e.abs.fit_movshon.R2, doc_a.abs.fit_movshon.R2, tolerance.abs.fit_movshon.R2, 'abs fit_movshon R2');
[b_(53),errormsg_{53}] = ndi.test.values_within_tolerance(doc_e.abs.fit_movshon.bandwidth, doc_a.abs.fit_movshon.bandwidth, tolerance.abs.fit_movshon.bandwidth, 'abs fit_movshon bandwidth');

%  Comparing fit_movshon_c
%    parameters
%    values
%    fit
%    L50
%    Pref
%    H50
%    R2
%    bandwidth

[b_(54),errormsg_{54}] = ndi.test.values_within_tolerance(doc_e.abs.fit_movshon_c.parameters, doc_a.abs.fit_movshon_c.parameters, tolerance.abs.fit_movshon_c.parameters, 'abs fit_movshon_c parameters');
[b_(55),errormsg_{55}] = ndi.test.values_within_tolerance(doc_e.abs.fit_movshon_c.values, doc_a.abs.fit_movshon_c.values, tolerance.abs.fit_movshon_c.values, 'abs fit_movshon_c values');
[b_(56),errormsg_{56}] = ndi.test.values_within_tolerance(doc_e.abs.fit_movshon_c.fit, doc_a.abs.fit_movshon_c.fit, tolerance.abs.fit_movshon_c.fit, 'abs fit_movshon_c fit');
[b_(57),errormsg_{57}] = ndi.test.values_within_tolerance(doc_e.abs.fit_movshon_c.L50, doc_a.abs.fit_movshon_c.L50, tolerance.abs.fit_movshon_c.L50, 'abs fit_movshon_c L50');
[b_(58),errormsg_{58}] = ndi.test.values_within_tolerance(doc_e.abs.fit_movshon_c.Pref, doc_a.abs.fit_movshon_c.Pref, tolerance.abs.fit_movshon_c.Pref, 'abs fit_movshon_c Pref');
[b_(59),errormsg_{59}] = ndi.test.values_within_tolerance(doc_e.abs.fit_movshon_c.H50, doc_a.abs.fit_movshon_c.H50, tolerance.abs.fit_movshon_c.H50, 'abs fit_movshon_c H50');
[b_(60),errormsg_{60}] = ndi.test.values_within_tolerance(doc_e.abs.fit_movshon_c.R2, doc_a.abs.fit_movshon_c.R2, tolerance.abs.fit_movshon_c.R2, 'abs fit_movshon_c R2');
[b_(61),errormsg_{61}] = ndi.test.values_within_tolerance(doc_e.abs.fit_movshon_c.bandwidth, doc_a.abs.fit_movshon_c.bandwidth, tolerance.abs.fit_movshon_c.bandwidth, 'abs fit_movshon_c bandwidth');

%  Comparing fit_spline
%    values
%    fit
%    L50
%    Pref
%    H50
%    bandwidth

[b_(62),errormsg_{62}] = ndi.test.values_within_tolerance(doc_e.abs.fit_spline.values, doc_a.abs.fit_spline.values, tolerance.abs.fit_spline.values, 'abs fit_spline values');
[b_(63),errormsg_{63}] = ndi.test.values_within_tolerance(doc_e.abs.fit_spline.fit, doc_a.abs.fit_spline.fit, tolerance.abs.fit_spline.fit, 'abs fit_spline fit');
[b_(64),errormsg_{64}] = ndi.test.values_within_tolerance(doc_e.abs.fit_spline.L50, doc_a.abs.fit_spline.L50, tolerance.abs.fit_spline.L50, 'abs fit_spline L50');
[b_(65),errormsg_{65}] = ndi.test.values_within_tolerance(doc_e.abs.fit_spline.Pref, doc_a.abs.fit_spline.Pref, tolerance.abs.fit_spline.Pref, 'abs fit_spline Pref');
[b_(66),errormsg_{66}] = ndi.test.values_within_tolerance(doc_e.abs.fit_spline.H50, doc_a.abs.fit_spline.H50, tolerance.abs.fit_spline.H50, 'abs fit_spline H50');
[b_(67),errormsg_{67}] = ndi.test.values_within_tolerance(doc_e.abs.fit_spline.bandwidth, doc_a.abs.fit_spline.bandwidth, tolerance.abs.fit_spline.bandwidth, 'abs fit_spline bandwidth');

%  Comparing fit_gausslog
%    parameters
%    values
%    fit
%    L50
%    Pref
%    H50
%    bandwidth

[b_(68),errormsg_{68}] = ndi.test.values_within_tolerance(doc_e.abs.fit_gausslog.parameters, doc_a.abs.fit_gausslog.parameters, tolerance.abs.fit_gausslog.parameters, 'abs fit_gausslog parameters');
[b_(69),errormsg_{69}] = ndi.test.values_within_tolerance(doc_e.abs.fit_gausslog.values, doc_a.abs.fit_gausslog.values, tolerance.abs.fit_gausslog.values, 'abs fit_gausslog values');
[b_(70),errormsg_{70}] = ndi.test.values_within_tolerance(doc_e.abs.fit_gausslog.fit, doc_a.abs.fit_gausslog.fit, tolerance.abs.fit_gausslog.fit, 'abs fit_gausslog fit');
[b_(71),errormsg_{71}] = ndi.test.values_within_tolerance(doc_e.abs.fit_gausslog.L50, doc_a.abs.fit_gausslog.L50, tolerance.abs.fit_gausslog.L50, 'abs fit_gausslog L50');
[b_(72),errormsg_{72}] = ndi.test.values_within_tolerance(doc_e.abs.fit_gausslog.Pref, doc_a.abs.fit_gausslog.Pref, tolerance.abs.fit_gausslog.Pref, 'abs fit_gausslog Pref');
[b_(73),errormsg_{73}] = ndi.test.values_within_tolerance(doc_e.abs.fit_gausslog.H50, doc_a.abs.fit_gausslog.H50, tolerance.abs.fit_gausslog.H50, 'abs fit_gausslog H50');
[b_(74),errormsg_{74}] = ndi.test.values_within_tolerance(doc_e.abs.fit_gausslog.bandwidth, doc_a.abs.fit_gausslog.bandwidth, tolerance.abs.fit_gausslog.bandwidth, 'abs fit_gausslog bandwidth');
                                                                    
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
