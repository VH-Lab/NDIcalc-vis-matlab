function [stats] = plot_tf(bigtable, condition_name, reference_group, group_name, prefix)
% PLOT_TF - visualize temporal frequency results across conditions
%
% 

prefix = 'best';

colors = [0 0 0 ; 1 0 0 ]; 
group_line_color = [ 1 1 0];


I = find(bigtable.best_tf_visual_response_anova_p<0.05);
bigtable_tf = bigtable(I,:);

Y_values = {[prefix '_tf_empirical_low_pass_index'],...
	[prefix '_tf_ultimate_Pref'],...
	[prefix '_tf_fitless_bandwidth'],...
	[prefix '_tf_empirical_high_pass_index'],...
	[prefix '_tf_empirical_max_response_value']};
Y_labels = {'TF low pass index','TF Pref (Hz)','TF Bandwidth','TF high pass index','TF Max Response'};

log_type = [ 0 1 0 0 0];
plot_type = [1 1 2 1 1];
stat_type = [ 1 1 2 1 1];

figlist_exist = get(0,'children');

stats.condition_name = condition_name;
stats.group_name = group_name;
stats.prefix = prefix;
stats.Y_values = Y_values;
stats.Y_labels = Y_labels;

[stats.lme,stats.lme_]=vlt.stats.plot_lme_array(bigtable_tf, condition_name, Y_values, Y_labels, ...
	{'Y','Y','vlt.math.clip(Y,[0 8])','Y','Y'}, reference_group, group_name, log_type, plot_type, stat_type,...
	'colors',colors,'category_mean_color',[0.5 0.5 0.5],'group_mean_color',group_line_color,'point_marker_size',2,...
	'within_category_space',2,'across_category_space',4);

new_figs = setdiff(get(0,'children'),figlist_exist);
for i=1:numel(new_figs),
    set(new_figs(i),'tag',['TF_' int2str(i)]);
end;
figlist_exist = get(0,'children');

