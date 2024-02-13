function output = test()
% test- test contrast bayes analysis
%
% OUTPUT = test();
% 

test_params.rmax = 30*(1+0.5);
test_params.c50 = 0.5;
test_params.n = 1;
test_params.s = 2;

resp_struct.contrast = [0 0.1 0.2 0.4 0.8 1];
resp_struct.mean_responses = test_params.rmax*vis.contrast.naka_rushton_func(...
	resp_struct.contrast,...
	test_params.c50,test_params.n,test_params.s);
resp_struct.num_trials = 5 * ones(size(resp_struct.contrast));

noise_mdl = [0.25 0.73];

param_grid.r100 = linspace(0.001,150,150);
param_grid.c50 =  linspace(0.1,1,10);
param_grid.n =    linspace(0.5,4,20);
param_grid.s =    linspace(1,2,11);

[output_struct,lik] = vis.bayes.naka_rushton.grid_proportional_noise(...
	param_grid, resp_struct, noise_mdl);

output.output_struct = output_struct;
output.lik = lik;
output.test_params = test_params;
output.resp_struct = resp_struct;
output.noise_mdl = noise_mdl;

figure;
vis.bayes.naka_rushton.plot(output.output_struct,output)

