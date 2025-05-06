function [b_, errormsg_] = hartley_reverse_correlation_compare_docs(document_expected, document_actual, scope)
% hartley_reverse_correlation_compare_docs
%
% [B, ERRORMSG] = ndi.calc.vis.test.hartley_reverse_correlation_compare_docs(DOC_EXPECTED, DOC_ACTUAL, SCOPE)
%

% Initialize b_ as a row vector of ones for tracking comparison results
% Initialize errormsg_ as an empty cell array to hold error messages

b_ = ones(1,11);
errormsg_ = cell(1,11);

 % establish scope-dependent tolerances
switch(scope),
    case 'standard',

        tolerance.stimulus_properties.M = 1e-6;
        tolerance.stimulus_properties.L_max = 1e-6;
        tolerance.stimulus_properties.K_max = 1e-6;
        tolerance.stimulus_properties.sf_max = 1e-6;
        tolerance.stimulus_properties.fps = 1e-6;
        tolerance.stimulus_properties.color_high = 1e-6;
        tolerance.stimulus_properties.color_low = 1e-6;
        tolerance.stimulus_properties.rect = 1e-6;
        tolerance.reconstruction_properties.T_coords = 1e-6;
        tolerance.reconstruction_properties.X_coords = 1e-6;
        tolerance.reconstruction_properties.Y_coords = 1e-6;

    case 'low_noise',

        tolerance.stimulus_properties.M = 1e-6;
        tolerance.stimulus_properties.L_max = 1e-6;
        tolerance.stimulus_properties.K_max = 1e-6;
        tolerance.stimulus_properties.sf_max = 1e-6;
        tolerance.stimulus_properties.fps = 1e-6;
        tolerance.stimulus_properties.color_high = 1e-6;
        tolerance.stimulus_properties.color_low = 1e-6;
        tolerance.stimulus_properties.rect = 1e-6;
        tolerance.reconstruction_properties.T_coords = 1e-6;
        tolerance.reconstruction_properties.X_coords = 1e-6;
        tolerance.reconstruction_properties.Y_coords = 1e-6;

    case 'high_noise',

        tolerance.stimulus_properties.M = 1e-6;
        tolerance.stimulus_properties.L_max = 1e-6;
        tolerance.stimulus_properties.K_max = 1e-6;
        tolerance.stimulus_properties.sf_max = 1e-6;
        tolerance.stimulus_properties.fps = 1e-6;
        tolerance.stimulus_properties.color_high = 1e-6;
        tolerance.stimulus_properties.color_low = 1e-6;
        tolerance.stimulus_properties.rect = 1e-6;
        tolerance.reconstruction_properties.T_coords = 1e-6;
        tolerance.reconstruction_properties.X_coords = 1e-6;
        tolerance.reconstruction_properties.Y_coords = 1e-6;

    otherwise,
       error(['Unknown scope ' scope '.']);
end;

% start comparison

doc_e = document_expected.document_properties.hartley_reverse_correlation;
doc_a = document_actual.document_properties.hartley_reverse_correlation;


% tol =  1e-6

% Comparing stimulus_properties
%   M					
%	L_max
%	K_max
%	sf_max
%	fps
%	color_high
%	color_low
%	rect

[b_(1),errormsg_{1}] = ndi.test.values_within_tolerance(doc_e.stimulus_properties.M, doc_a.stimulus_properties.M, tolerance.stimulus_properties.M, 'M');
[b_(2),errormsg_{2}] = ndi.test.values_within_tolerance(doc_e.stimulus_properties.L_max, doc_a.stimulus_properties.L_max, tolerance.stimulus_properties.L_max, 'L max');
[b_(3),errormsg_{3}] = ndi.test.values_within_tolerance(doc_e.stimulus_properties.K_max, doc_a.stimulus_properties.K_max, tolerance.stimulus_properties.K_max, 'K max');
[b_(4),errormsg_{4}] = ndi.test.values_within_tolerance(doc_e.stimulus_properties.sf_max, doc_a.stimulus_properties.sf_max, tolerance.stimulus_properties.sf_max, 'sf max');
[b_(5),errormsg_{5}] = ndi.test.values_within_tolerance(doc_e.stimulus_properties.fps, doc_a.stimulus_properties.fps, tolerance.stimulus_properties.fps, 'fps');
[b_(6),errormsg_{6}] = ndi.test.values_within_tolerance(doc_e.stimulus_properties.color_high, doc_a.stimulus_properties.color_high, tolerance.stimulus_properties.color_high, 'color high');
[b_(7),errormsg_{7}] = ndi.test.values_within_tolerance(doc_e.stimulus_properties.color_low, doc_a.stimulus_properties.color_low, tolerance.stimulus_properties.color_low, 'color low');
[b_(8),errormsg_{8}] = ndi.test.values_within_tolerance(doc_e.stimulus_properties.rect, doc_a.stimulus_properties.rect, tolerance.stimulus_properties.rect, 'rect');




% Comparing reconstruction_properties
%   T_coords
%	X_coords
%	Y_coords

% tol = 1e-6

[b_(9),errormsg_{9}] = ndi.test.values_within_tolerance(doc_e.reconstruction_properties.T_coords, doc_a.reconstruction_properties.T_coords, tolerance.reconstruction_properties.T_coords, 'T coords');
[b_(10),errormsg_{10}] = ndi.test.values_within_tolerance(doc_e.reconstruction_properties.X_coords, doc_a.reconstruction_properties.X_coords, tolerance.reconstruction_properties.X_coords, 'X coords');
[b_(11),errormsg_{11}] = ndi.test.values_within_tolerance(doc_e.reconstruction_properties.Y_coords, doc_a.reconstruction_properties.Y_coords, tolerance.reconstruction_properties.Y_coords, 'Y coords');

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

