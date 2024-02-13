function plot(output_struct)
% PLOT - plot the output of the Bayesian Double Gaussian grid model
%
% vis.bayes.naka_rushton.plot(OUTPUT_STRUCT)
%
% Plot the output of the Bayesian Naka_Rushton parameter estimation
% in the current figure.
%
% See also: vis.bayes.naka_rushton.grid_proportional_noise()
%

subplot(3,3,4);
plot(output_struct.marginal_likelihoods.r100.values,output_struct.marginal_likelihoods.r100.likelihoods,'b-o'),
set(gca,'xscale','log')
xlabel("R100 (Hz)"),
ylabel("probability of R100");
box off;

subplot(3,3,5);
plot(output_struct.marginal_likelihoods.c50.values,output_struct.marginal_likelihoods.c50.likelihoods,'b-o'),
xlabel("c50"),
ylabel("probability of c50");
box off;

subplot(3,3,6);
plot(output_struct.marginal_likelihoods.n.values,output_struct.marginal_likelihoods.n.likelihoods,'b-o'),
xlabel("n"),
ylabel("probability of n");
box off;

subplot(3,3,7);
plot(output_struct.marginal_likelihoods.s.values,output_struct.marginal_likelihoods.s.likelihoods,'b-o'),
xlabel("s"),
ylabel("probability of s");
box off;

subplot(3,3,1);
plot(output_struct.other_parameters.independent_variable_value, output_struct.other_parameters.mean_responses,'b-o');
set(gca,'xscale','linear')
xlabel("Contrast"),
ylabel("Response");
hold on
c_high_res = 0:0.01:1;
p=output_struct.maximum_likelihood_parameters.parameters;
nr_high_res = p.r100*(1+0.5)*vis.contrast.naka_rushton_func(c_high_res,p.c50,p.n,p.s);
plot(c_high_res,nr_high_res,'r--'); 
box off;

subplot(3,3,3);
plot(output_struct.descriptors.RelativeMaximumGain.values,output_struct.descriptors.RelativeMaximumGain.likelihoods,'b-o'),
xlabel("Relative maximum gain"),
ylabel("probability of RMG");
box off;

subplot(3,3,8);
plot(output_struct.descriptors.ContrastThreshold.values,output_struct.descriptors.ContrastThreshold.likelihoods,'b-o'),
xlabel("Contrast Threshold"),
ylabel("probability of Contrast Threshold");
box off;


subplot(3,3,9);
plot(output_struct.descriptors.SaturationIndex.values,output_struct.descriptors.SaturationIndex.likelihoods,'b-o'),
xlabel("Saturation Index"),
ylabel("Prob of Saturation index");
box off;


return;

subplot(3,3,9);
plot(output_struct.descriptors.dir_cv.values,output_struct.descriptors.dir_cv.likelihoods,'b-o'),
xlabel("DCV"),
ylabel("probability of DCV");

box off;
