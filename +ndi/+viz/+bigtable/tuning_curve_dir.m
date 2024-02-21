function big_dir_table = tuning_curve_dir(S, element_name,element_ref, condition_name, tuning_curve_search_string, element_type);

fit_doc_type = 'oridirtuning_calc';
fit_type_property = 'orientation_direction_tuning';
property_name = 'orientation_direction_tuning';
x_axis_name = 'direction';
test_prefix = 'DIR';
table_func = 'ndi.viz.bigtable.dir_struct';

if nargin<5,
	tuning_curve_search_string = 'angle';
end;

if nargin<6,
	element_type = 'spikes';
end;

big_dir_table = ndi.viz.bigtable.tuning_curve(S,element_name,element_ref, table_func,...
	'fit_doc_type',fit_doc_type,...
	'tuning_curve_search_string',tuning_curve_search_string,...
	'fit_type_property', fit_type_property,...
	'property_name',property_name,...
	'x_axis_name', x_axis_name, ...
	'test_prefix',test_prefix,...
	'element_type',element_type);

for i=1:size(big_dir_table,1),
	condition_struct(i).condition_name = condition_name;
end;

big_dir_table = [struct2table(condition_struct) big_dir_table];


