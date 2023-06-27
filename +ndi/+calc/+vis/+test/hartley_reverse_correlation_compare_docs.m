function [b, errormsg] = hartley_reverse_correlation_compare_docs(document_expected, document_actual, scope)
% hartley_reverse_correlation_compare_docs
%
% [B, ERRORMSG] = ndi.calc.vis.test.hartley_reverse_correlation_compare_docs(DOC_EXPECTED, DOC_ACTUAL, SCOPE)
%


b = 1;
errormsg = '';

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

[b,errormsg] = ndi.test.values_within_tolerance(doc_e.stimulus_properties.M, doc_a.stimulus_properties.M, tolerance, 'M');
[b,errormsg] = ndi.test.values_within_tolerance(doc_e.stimulus_properties.L_max, doc_a.stimulus_properties.L_max, tolerance, 'L max');
[b,errormsg] = ndi.test.values_within_tolerance(doc_e.stimulus_properties.K_max, doc_a.stimulus_properties.K_max, tolerance, 'K max');
[b,errormsg] = ndi.test.values_within_tolerance(doc_e.stimulus_properties.sf_max, doc_a.stimulus_properties.sf_max, tolerance, 'sf max');
[b,errormsg] = ndi.test.values_within_tolerance(doc_e.stimulus_properties.fps, doc_a.stimulus_properties.fps, tolerance, 'fps');
[b,errormsg] = ndi.test.values_within_tolerance(doc_e.stimulus_properties.color_high, doc_a.stimulus_properties.color_high, tolerance, 'color high');
[b,errormsg] = ndi.test.values_within_tolerance(doc_e.stimulus_properties.color_low, doc_a.stimulus_properties.color_low, tolerance, 'color low');
[b,errormsg] = ndi.test.values_within_tolerance(doc_e.stimulus_properties.rect, doc_a.stimulus_properties.rect, tolerance, 'rect');

% Comparing reconstruction_properties
%  T_coords
%	X_coords
%	Y_coords

[b,errormsg] = ndi.test.values_within_tolerance(doc_e.reconstruction_properties.T_coords, doc_a.reconstruction_properties.T_coords, tolerance, 'T coords');
[b,errormsg] = ndi.test.values_within_tolerance(doc_e.reconstruction_properties.X_coords, doc_a.reconstruction_properties.X_coords, tolerance, 'X coords');
[b,errormsg] = ndi.test.values_within_tolerance(doc_e.reconstruction_properties.Y_coords, doc_a.reconstruction_properties.Y_coords, tolerance, 'Y coords');

end

