function [stats] = plot_oridir(bigtable, condition_name, reference_group, group_name, options)
% PLOT_ORIDIR - visualize temporal frequency results across conditions
%
% STATS = PLOT_ORIDIR(TBL, CONDITION_NAME, REFERENCE_GROUP, GROUP_NAME, ...)
%
% Given a table TBL and table column name that describes the different experimental
% conditions CONDITION_NAME and the name of the REFERENCE_GROUP (the condition
% that is the control) and the name of a column that determines random factors
% (GROUP_NAME), plots many features relevant for ORIENTATION/DIRECTION tuning.
%
% The fields that are examined are as follows:
%   [prefix '_DIR_empirical_low_pass_index']
%   [prefix '_DIR_ultimate_Pref']
%   [prefix '_DIR_fitless_bandwidth']
%   [prefix '_DIR_empirical_high_pass_index']
%   [prefix '_DIR_empirical_max_response_value']
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


I = find(bigtable.best_DIR_visual_response_anova_p<0.05);
bigtable_dir = bigtable(I,:);

Y_values = {...
	[prefix '_DIR_vector_circular_variance'],...
	[prefix '_DIR_vector_direction_circular_variance'],...
	[prefix '_DIR_empirical_max_response_value']};
Y_labels = {'1 - CircVar','1 -DirCircVar',...
	'DIR Max Response'};

log_type = [ 0 0 0];
plot_type = [1 1 1];
stat_type = [ 1 1 1];


figlist_exist = get(0,'children');

stats.condition_name = condition_name;
stats.group_name = group_name;
stats.prefix = prefix;
stats.Y_values = Y_values;
stats.Y_labels = Y_labels;

[stats.lme,stats.lme_]=vlt.stats.plot_lme_array(bigtable_dir, condition_name,...
	 Y_values, Y_labels, {'1-Y','1-Y','Y'}, ...
	reference_group, group_name, log_type, plot_type, stat_type,...
	'colors',colors,'category_mean_color',[0.5 0.5 0.5],...
	'group_mean_color',group_line_color,'point_marker_size',2,...
	'within_category_space',2,'across_category_space',4);

new_figs = setdiff(get(0,'children'),figlist_exist);
for i=1:numel(new_figs),
    set(new_figs(i),'tag',['TF_' int2str(i)]);
end;
figlist_exist = get(0,'children');

