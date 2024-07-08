function big_sf_table = tuning_curve_sf(S, element_name,element_ref, condition_name, tuning_curve_search_string, element_type, acceptable_response_types);

fit_doc_type = 'spatial_frequency_tuning_calc';
fit_type_property = 'spatial_frequency_tuning';
property_name = 'spatial_frequency_tuning';
x_axis_name = 'spatial_frequency';
test_prefix = 'sf';
table_func = 'ndi.viz.bigtable.tc_table';

if nargin<5,
	tuning_curve_search_string = 'Spatial Frequency (best DIR)';
end;

if nargin<6,
	element_type = 'spikes';
end;

if nargin<7,
    acceptable_response_types = {'mean','F1'};
end;

big_sf_table = ndi.viz.bigtable.tuning_curve2(S,element_name,element_ref, table_func,...
	'tuning_curve_search_string',tuning_curve_search_string, ...
	'fit_type_property', fit_type_property, ...
	'fit_doc_type', fit_doc_type, ...
	'property_name',property_name, ...
	'element_type',element_type, ...
	'x_axis_name', x_axis_name, ...
	'test_prefix','sf');

for i=1:size(big_sf_table,1),
	condition_struct(i).condition_name = condition_name;
end;

big_sf_table = [struct2table(condition_struct) big_sf_table];

