function [thestruct] = freq_struct(tuning_doc, stimulus_response_doc, response_type, f1f0_struct, varargin)


prefix = '';
property_name = 'temporal_frequency_tuning';
x_axis_name = 'temporal_frequency';

did.datastructures.assign(varargin{:});

tune_info = getfield(tuning_doc.document_properties, property_name);

thestruct.tuning_curve_doc = tuning_doc.document_properties.base.id;
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
thestruct.empirical_control_mean = mean(tune_info.tuning_curve.control_mean);
thestruct.f0_empirical = f1f0_struct.f0;
thestruct.f1_empirical = f1f0_struct.f1;
thestruct.f1f0_ratio = f1f0_struct.f1f0_ratio;
thestruct.f1f0_2f1overf1f0 = f1f0_struct.f1f0_2f1overf1f0;

thestruct.empirical_low_pass_index = tune_info.tuning_curve.mean(minx_loc)/maxvalue;
thestruct.empirical_high_pass_index = tune_info.tuning_curve.mean(maxx_loc)/maxvalue;

thestruct.fitless_L50 = tune_info.fitless.L50;
thestruct.fitless_Pref = tune_info.fitless.Pref;
thestruct.fitless_H50 = tune_info.fitless.H50;
thestruct.fitless_bandwidth = vlt.math.bandwidth(thestruct.fitless_L50,thestruct.fitless_H50);


thestruct.fit_dog_L50 = tune_info.fit_dog.L50;
thestruct.fit_dog_Pref = tune_info.fit_dog.Pref;
thestruct.fit_dog_H50 = tune_info.fit_dog.H50;
thestruct.fit_dog_bandwidth = vlt.math.bandwidth(thestruct.fit_dog_L50,thestruct.fit_dog_H50);

thestruct.fit_dog_offset = tune_info.fit_dog.parameters(1);
thestruct.fit_dog_a1 = tune_info.fit_dog.parameters(2);
thestruct.fit_dog_b1 = tune_info.fit_dog.parameters(3);
thestruct.fit_dog_a2 = tune_info.fit_dog.parameters(4);
thestruct.fit_dog_b2 = tune_info.fit_dog.parameters(5);

if thestruct.fit_dog_Pref > thestruct.high_x_value | thestruct.fit_dog_Pref < thestruct.low_x_value,
	% if the fit is bad, use the empirical
	thestruct.ultimate_L50 = thestruct.fitless_L50;
	thestruct.ultimate_Pref = thestruct.fitless_Pref;
	thestruct.ultimate_H50 = thestruct.fitless_H50;
	thestruct.ultimate_bandwidth = thestruct.fitless_bandwidth;
	thestruct.ultimate_method = 'empirical';
else, % use the fit
	thestruct.ultimate_L50 = thestruct.fit_dog_L50;
	thestruct.ultimate_Pref = thestruct.fit_dog_Pref;
	thestruct.ultimate_H50 = thestruct.fit_dog_H50;
	thestruct.ultimate_bandwidth = thestruct.fit_dog_bandwidth;
	thestruct.ultimate_method = 'dog_fit';
end;

thestruct = cell2struct(struct2cell(thestruct), strcat(prefix, fieldnames(thestruct)));
