function [thestruct] = dir_struct(fit_doc, stimulus_response_doc, response_type, f1f0_struct, varargin)


prefix = '';
property_name = 'orientation_direction_tuning';
x_axis_name = 'direction';

did.datastructures.assign(varargin{:});

tune_info = getfield(fit_doc.document_properties, property_name);

thestruct.tuning_curve_doc = fit_doc.document_properties.base.id;
thestruct.stimulus_response_doc = stimulus_response_doc.document_properties.base.id; 
thestruct.response_type = response_type;
thestruct.visual_response_anova_p = tune_info.significance.visual_response_anova_p;
thestruct.across_stimuli_anova_p = tune_info.significance.across_stimuli_anova_p;
x_values = getfield(tune_info.tuning_curve,x_axis_name);
[minx,minx_loc] = min(x_values);
[maxx,maxx_loc] = max(x_values);
thestruct.low_x_value = minx;
thestruct.high_x_value = maxx;
thestruct.x_values = mat2str(x_values);
thestruct.mean_values = mat2str(tune_info.tuning_curve.mean);
thestruct.stddev_values = mat2str(tune_info.tuning_curve.stddev);
thestruct.stderr_values = mat2str(tune_info.tuning_curve.stderr);
[maxvalue,maxvalue_loc] = max(tune_info.tuning_curve.mean);
thestruct.empirical_max_response_value = maxvalue;
thestruct.empirical_max_response_location = x_values(maxvalue_loc);
[minvalue,minvalue_loc] = min(tune_info.tuning_curve.mean);
thestruct.empirical_min_response_value = minvalue;
thestruct.empirical_min_response_location = x_values(minvalue_loc);
thestruct.empirical_control_individual = mat2str(tune_info.tuning_curve.control_individual);

if ~isempty(f1f0_struct),
	thestruct.f0_empirical = f1f0_struct.f0;
	thestruct.f1_empirical = f1f0_struct.f1;
	thestruct.f1f0_ratio = f1f0_struct.f1f0_ratio;
	thestruct.f1f0_2f1overf1f0 = f1f0_struct.f1f0_2f1overf1f0;
else,
	thestruct.f0_empirical = [];
	thestruct.f1_empirical = [];
	thestruct.f1f0_ratio = [];
	thestruct.f1f0_2f1overf1f0 = [];
end;

fn = fieldnames(tune_info.vector);
for i=1:numel(fn), 
	thestruct = setfield(thestruct,['vector_' fn{i}], getfield(tune_info.vector,fn{i}));
end;

fn = fieldnames(tune_info.fit);
for i=1:numel(fn),
	if any(strcmp(fn{i},{'double_gaussian_parameters','double_gaussian_fit_angles','double_gaussian_fit_values'})),
		thestruct = setfield(thestruct,['fit_' fn{i}],mat2str(getfield(tune_info.fit,fn{i})));
	else,
		thestruct = setfield(thestruct,['fit_' fn{i}],getfield(tune_info.fit,fn{i}));
	end;
end;

thestruct = cell2struct(struct2cell(thestruct), strcat(prefix, fieldnames(thestruct)));


