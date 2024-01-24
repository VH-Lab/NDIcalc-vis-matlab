function plot_bayes_sftf(sf_tf_list, bayes_curves)
% viz.plot_bayes_sftf - plot an array of bayesian theta_pref likelihoods for a range of SF and TF preferences
%
% Plots in current axes
%
%
%
%



prop_name = '.document_properties.orientation_direction_bayes.marginal_likelihoods.theta_pref';
prop_base_name = '.document_properties.stimulus_tuningcurve.';
prop_name_ind = [prop_base_name 'independent_variable_value'];
prop_name_values = [prop_base_name 'response_mean'];


angles = eval(['bayes_curves{1,1}' prop_name_ind ]);

dx = mean(diff(angles));
X_total = angles(end) - angles(1);
dy = 0.1;

% step 1: find the maximum value of the data

M = -Inf;

for i=1:size(bayes_curves,1),
	for j=1:size(bayes_curves,2),
		M = max(M,max(eval(['bayes_curves{i,j}' prop_name_values])));
		angles = eval(['bayes_curves{i,j}' prop_name_values]);
		if ~isequal(angles,angles),
			error(['All angles are not equal.']);
		end;
	end;
end;

for i=1:size(bayes_curves,1),
	for j=1:size(bayes_curves,2),
    		angles = eval(['bayes_curves{i,j}' prop_name_ind]);
		v = eval(['bayes_curves{i,j}' prop_name_values]);
		hold on;
                plot( (i-1)*(dx+X_total)+angles, (j-1)*(dy+1)+0*v/M,'k--');
		plot( (i-1)*(dx+X_total)+angles, (j-1)*(dy+1)+v/M);
	end;
end;

title(['Max value: ' num2str(M)]);

