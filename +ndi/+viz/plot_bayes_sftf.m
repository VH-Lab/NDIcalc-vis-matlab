function plot_bayes_sftf(sf_tf_list, bayes_curves, bayes_prop)
% viz.plot_bayes_sftf - plot an array of bayesian theta_pref likelihoods for a range of SF and TF preferences
%
% Plots in current axes
%
%
%
%

if nargin<3,
	bayes_prop = 'theta_pref';
end;

prop_name = ['.document_properties.orientation_direction_bayes.marginal_likelihoods.' bayes_prop];
likelihood_angles = eval(['bayes_curves{1,1}' prop_name '.values']);

dx = mean(diff(likelihood_angles));

if strcmp(bayes_prop,'theta_pref'),
    dx = 45;
end;

X_total = likelihood_angles(end) - likelihood_angles(1);
dy = 0.1;

% step 1: find the maximum value of the data

M = -Inf;

for i=1:size(bayes_curves,1),
	for j=1:size(bayes_curves,2),
		M = max(M,max(eval(['bayes_curves{i,j}' prop_name '.likelihoods'])));
		angles = eval(['bayes_curves{i,j}' prop_name '.values']);
		if ~isequal(angles,likelihood_angles),
			error(['All angles are not equal.']);
		end;
	end;
end;

for i=1:size(bayes_curves,1),
	for j=1:size(bayes_curves,2),
    		likelihood_angles = eval(['bayes_curves{i,j}' prop_name '.values']);
		l = eval(['bayes_curves{i,j}' prop_name '.likelihoods']);
		hold on;
		plot( (i-1)*(dx+X_total)+likelihood_angles, (j-1)*(dy+1)+0*l/M,'k--');
		plot( (i-1)*(dx+X_total)+likelihood_angles, (j-1)*(dy+1)+l/M);
	end;
end;


title(['Max value: ' num2str(M)]);
