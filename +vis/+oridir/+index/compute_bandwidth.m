function [bandwidth]= compute_bandwidth( x, y )

% Calculating bandwidth
% The tuning curve is calculated before this bandwidth can be found
% x is tc_doc_to_examine.document_properties.stimulus_tuningcurve.independent_variable_value
% y is tc_doc_to_examine.document_properties.stimulus_tuningcurve.response_mean
%

[low, ~, high] = compute_halfwidth(x, y);

if isnan(low) || isnan(high)
    bandwidth = Inf;
else
    bandwidth = log2(high/low);
end