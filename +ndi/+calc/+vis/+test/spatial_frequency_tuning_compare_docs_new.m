function [b_, errormsg_] = spatial_frequency_tuning_compare_docs(document_expected, document_actual, scope)
% Comparing spatial_frequency_tuning
%
% [B, ERRORMSG] = ndi.calc.vis.test.spatial_frequency_tuning_compare_docs(DOC_EXPECTED, DOC_ACTUAL, SCOPE)
%
%

% Initialize b_ as a row vector of ones for tracking comparison results
% Initialize errormsg_ as an empty cell array to hold error messages

b_ = ones(1,39);
errormsg_ = cell(1,39);

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
       tol_fitless.low_pass_index = 1e-6;
       tol_fitless.high_pass_index = 1e-6;
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
       
       %tol_fit_abs still needs to be added

    case 'low_noise',
    
       tol_tuning_curve.spatial_frequency = 1e-6;
       tol_tuning_curve.mean = 1e-6;
       tol_tuning_curve.stddev = 1e-6;
       tol_tuning_curve.stderr = 1e-6;
       tol_tuning_curve.individual = 1e-6;
       tol_tuning_curve.control_stddev = 1e-6;
       tol_tuning_curve.control_stderr = 1e-6;
       tol_significance.visual_response_anova_p = 1e-6;
       tol_significance.across_stimuli_anova_p = 1e-6;
       tol_fitless.H50 = 0.1;
       tol_fitless.Pref = 0.1;
       tol_fitless.L50 = 0.1;
       tol_fit_spline.fit = 0.1;
       tol_fit_spline.H50 = 0.1;
       tol_fit_spline.Pref = 0.1;
       tol_fit_spline.L50 = 0.1;
       tol_fit_spline.values = 0.1;
       tol_fit_dog.fit = 0.1;
       tol_fit_dog.H50 = 0.1;
       tol_fit_dog.Pref = 0.1;
       tol_fit_dog.L50 = 0.1;
       tol_fit_dog.values = 0.1;
       tol_fit_dog.parameters = 0.1;
       tol_fit_gausslog.fit = 0.1;
       tol_fit_gausslog.H50 = 0.1;
       tol_fit_gausslog.Pref = 0.1;
       tol_fit_gausslog.L50 = 0.1;
       tol_fit_gausslog.values = 0.1;
       tol_fit_gausslog.parameters = 0.1;

    case 'high_noise',
    
       tol_tuning_curve.spatial_frequency = 1e-6;
       tol_tuning_curve.mean = 1e-6;
       tol_tuning_curve.stddev = 1e-6;
       tol_tuning_curve.stderr = 1e-6;
       tol_tuning_curve.individual = 1e-6;
       tol_tuning_curve.control_stddev = 1e-6;
       tol_tuning_curve.control_stderr = 1e-6;
       tol_significance.visual_response_anova_p = 1e-6;
       tol_significance.across_stimuli_anova_p = 1e-6;
       tol_fitless.H50 = 0.1;
       tol_fitless.Pref = 0.1;
       tol_fitless.L50 = 0.1;
       tol_fit_spline.fit = 0.1;
       tol_fit_spline.H50 = 0.1;
       tol_fit_spline.Pref = 0.1;
       tol_fit_spline.L50 = 0.1;
       tol_fit_spline.values = 0.1;
       tol_fit_dog.fit = 0.1;
       tol_fit_dog.H50 = 0.1;
       tol_fit_dog.Pref = 0.1;
       tol_fit_dog.L50 = 0.1;
       tol_fit_dog.values = 0.1;
       tol_fit_dog.parameters = 0.1;
       tol_fit_gausslog.fit = 0.1;
       tol_fit_gausslog.H50 = 0.1;
       tol_fit_gausslog.Pref = 0.1;
       tol_fit_gausslog.L50 = 0.1;
       tol_fit_gausslog.values = 0.1;
       tol_fit_gausslog.parameters = 0.1;

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
   b_(38) = 0;
   errormsg_{38} = ['Expected response units in ' doc_e.properties.response_units ' but observed ' doc_a.properties.response_units];
   return;
end

%   Response Type

properties_match = strcmpi(char(doc_e.properties.response_type), char(doc_a.properties.response_type));
if ~properties_match
   b_(39) = 0;
   errormsg_{39} = ['Expected response type of ' doc_e.properties.response_type ' but observed ' doc_a.properties.response_type];
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
[b_(23),errormsg_{23}] = ndi.test.values_within_tolerance(doc_e.fit_dog.parameters, doc_a.fit_dog.parameters, tol_fit_dog.parameters, 'fit_dog parameters');

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
[b_(29),errormsg_{29}] = ndi.test.values_within_tolerance(doc_e.fit_dog.parameters, doc_a.fit_dog.parameters, tol_fit_gausslog.parameters, 'fit_gausslog parameters');
                                                                    
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
