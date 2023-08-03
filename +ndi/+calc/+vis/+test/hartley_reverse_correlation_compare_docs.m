function [b_, errormsg_] = hartley_reverse_correlation_compare_docs(document_expected, document_actual, scope)
% hartley_reverse_correlation_compare_docs
%
% [B, ERRORMSG] = ndi.calc.vis.test.hartley_reverse_correlation_compare_docs(DOC_EXPECTED, DOC_ACTUAL, SCOPE)
%

% Initialize b_ as a row vector of ones for tracking comparison results
% Initialize errormsg_ as an empty cell array to hold error messages

b_ = ones(1,11);
errormsg_ = cell(1,11);

% start comparison

doc_e = document_expected.document_properties.hartley_reverse_correlation;
doc_a = document_actual.document_properties.hartley_reverse_correlation;

% Comparing stimulus_properties
%  M					
%	L_max
%	K_max
%	sf_max
%	fps
%	color_high
%	color_low
%	rect

[b_(1),errormsg_{1}] = ndi.test.values_within_tolerance(doc_e.stimulus_properties.M, doc_a.stimulus_properties.M, tolerance, 'M');
[b_(2),errormsg_{2}] = ndi.test.values_within_tolerance(doc_e.stimulus_properties.L_max, doc_a.stimulus_properties.L_max, tolerance, 'L max');
[b_(3),errormsg_{3}] = ndi.test.values_within_tolerance(doc_e.stimulus_properties.K_max, doc_a.stimulus_properties.K_max, tolerance, 'K max');
[b_(4),errormsg_{4}] = ndi.test.values_within_tolerance(doc_e.stimulus_properties.sf_max, doc_a.stimulus_properties.sf_max, tolerance, 'sf max');
[b_(5),errormsg_{5}] = ndi.test.values_within_tolerance(doc_e.stimulus_properties.fps, doc_a.stimulus_properties.fps, tolerance, 'fps');
[b_(6),errormsg_{6}] = ndi.test.values_within_tolerance(doc_e.stimulus_properties.color_high, doc_a.stimulus_properties.color_high, tolerance, 'color high');
[b_(7),errormsg_{7}] = ndi.test.values_within_tolerance(doc_e.stimulus_properties.color_low, doc_a.stimulus_properties.color_low, tolerance, 'color low');
[b_(8),errormsg_{8}] = ndi.test.values_within_tolerance(doc_e.stimulus_properties.rect, doc_a.stimulus_properties.rect, tolerance, 'rect');

% Comparing reconstruction_properties
%  T_coords
%	X_coords
%	Y_coords

[b_(9),errormsg_{9}] = ndi.test.values_within_tolerance(doc_e.reconstruction_properties.T_coords, doc_a.reconstruction_properties.T_coords, tolerance, 'T coords');
[b_(10),errormsg_{10}] = ndi.test.values_within_tolerance(doc_e.reconstruction_properties.X_coords, doc_a.reconstruction_properties.X_coords, tolerance, 'X coords');
[b_(11),errormsg_{11}] = ndi.test.values_within_tolerance(doc_e.reconstruction_properties.Y_coords, doc_a.reconstruction_properties.Y_coords, tolerance, 'Y coords');

% Identify the b_ values with unmatched results
% Update b_ to only include those
% Update the corresponding errormsg_ messages

if any(b_==0),
    error_indices = find(b_==0);
    b_ = b_(error_indices);
    errormsg_ = errormsg_(error_indices);
end

end

