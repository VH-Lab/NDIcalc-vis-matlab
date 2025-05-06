function [stats] = plot_tf(bigtable, condition_name, reference_group, group_name, options)
% PLOT_TF - visualize temporal frequency results across conditions
%
% STATS = PLOT_TF(TBL, CONDITION_NAME, REFERENCE_GROUP, GROUP_NAME, ...)
%
% Given a table TBL and table column name that describes the different experimental
% conditions CONDITION_NAME and the name of the REFERENCE_GROUP (the condition
% that is the control) and the name of a column that determines random factors
% (GROUP_NAME), plots many features relevant for TEMPORAL FREQUENCY tuning.
%
% The fields that are examined are as follows:
%   [prefix '_tf_empirical_low_pass_index']
%   [prefix '_tf_ultimate_Pref']
%   [prefix '_tf_fitless_bandwidth']
%   [prefix '_tf_empirical_high_pass_index']
%   [prefix '_tf_empirical_max_response_value']
%   
%
% The function takes options as name/value pairs:
% ----------------------------------------------------------------------------------
% | Parameter (default)          | Description                                     |
% |------------------------------|-------------------------------------------------|
% |prefix ('best')               | The prefix to be applied in front of temporal   |
% |                              |    frequency tuning fields.                     |
% | colors (vlt.plot.colorlist())| The colors that should be assigned to the       |
% |                              |    conditions.                                  |
% | group_line_color([1 1 0])    | The group mean line color.                      |
% |------------------------------|-------------------------------------------------|
%

arguments
	bigtable (:,:)	table
	condition_name (1,:) char
	reference_group (1,:) char
	group_name (1,:) 
	options.prefix (1,:) char = 'best'
	options.colors (:,3) double = vlt.plot.colorlist
	options.group_line_color (1,3) double = [ 1 1 0];
end

prefix = options.prefix;
colors = options.colors;
group_line_color = options.group_line_color;

sig_values = [prefix '_tf_TF_tuning.sig.visual_response_anova_p'];
I = find(bigtable.(sig_values) <0.05);
bigtable_tf = bigtable(I,:);

Y_values = { [prefix '_tf_TF_tuning.TC.mean'],...
        [prefix '_tf_TF_tuning.TC.mean'],...
        [prefix '_tf_TF_tuning.TC.control_mean']};
Y_labels = {'Max response','Suppression','Control resp'};

log_type = [ 1 1 1 ];
plot_type = [1 1 1 ];
stat_type = [2 2 2 ];

figlist_exist = get(0,'children');

stats.condition_name = condition_name;
stats.group_name = group_name;
stats.prefix = prefix;
stats.Y_values = Y_values;
stats.Y_labels = Y_labels;

max_code = ['0.1+vlt.data.cell2matn(cellfun(@rectify,cellfun(@max,Y,''UniformOutput'',false),''UniformOutput'',false));'];
min_code = ['0.1-vlt.data.cell2matn(cellfun(@min,cellfun(@(x) -rectify(-x),Y,''UniformOutput'',false),''UniformOutput'',false));'];
median_code = ['0.1+vlt.data.cell2matn(cellfun(@median,Y,''UniformOutput'',false));'];

[stats.lme,stats.lme_]=vlt.stats.plot_lme_array(bigtable_tf, condition_name, Y_values, Y_labels, ...
	{max_code, min_code, median_code}, reference_group, group_name, log_type, plot_type, stat_type,...
	'colors',colors,'category_mean_color',[0.5 0.5 0.5],'group_mean_color',group_line_color,'point_marker_size',2,...
	'within_category_space',2,'across_category_space',4);

new_figs = setdiff(get(0,'children'),figlist_exist);
for i=1:numel(new_figs),
    set(new_figs(i),'tag',['respmag_' int2str(i)]);
end;
figlist_exist = get(0,'children');

