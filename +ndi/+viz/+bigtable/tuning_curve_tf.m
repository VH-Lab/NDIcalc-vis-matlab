function big_tf_table = tuning_curve_tf(S, element_name,element_ref, condition_name, tuning_curve_search_string, element_type);

fit_doc_type = 'temporal_frequency_tuning_calc';
fit_type_property = 'temporal_frequency_tuning';
property_name = 'temporal_frequency_tuning';
x_axis_name = 'temporal_frequency';
test_prefix = 'tf';
table_func = 'ndi.viz.bigtable.freq_struct';

if nargin<5,
	tuning_curve_search_string = 'Temporal_Frequency';
end;

if nargin<6,
	element_type = 'spikes';
end;

big_tf_table = ndi.viz.bigtable.tuning_curve(S,element_name,element_ref, table_func,...
	'tuning_curve_search_string',tuning_curve_search_string,...
	'fit_type_property', fit_type_property,...
	'fit_doc_type',fit_doc_type,...
	'property_name',property_name,...
	'x_axis_name', x_axis_name, ...
	'test_prefix',test_prefix,...
	'element_type',element_type);

for i=1:size(big_tf_table,1),
	condition_struct(i).condition_name = condition_name;
end;

big_tf_table = [struct2table(condition_struct) big_tf_table];

