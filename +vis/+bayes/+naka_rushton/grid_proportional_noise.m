function [output_struct,Lik] = grid_proportional_noise(param_grid, resp_struct, noise_mdl)

%% DESCRIPTION
% - Uses supersaturation-modified Naka-Rusthon from Pierce 2007 without parameter B as mean responses have baseline response subtracted
% - Limits on parameters set according to Cortes et al. 2022
% - Response variance is proportional to mean response for pooled data in log space
%
% INPUTS:
%   - PARAM_GRID: A structure with fields:
%       r100 - the values of r100 to examine
%       c50 - the values of c50 to examine
%       n - the values of n to examine
%       s - the values of s to examine
%   - RESP_STRUCT: a structure with fields:
%       contrast - the contrast that were used for stimulation
%       mean_responses - the mean responses to each stimulus
%       num_trials - the number of trials of each stimulus
%   - NOISE_MDL: a log linear regression model of noise for the given data [offset slope]
%
% OUTPUTS:
%   - OUTPUT_STRUCT: a structure with fields:
%       noise_model - details type of fit, offset, and slope
%       other_parameters - independent variable details
%       marginal_likelihods - marginal likelihoods of parameters
%       maximum_likelihoods_parameters - maximum likelihoods of parameters
%       descriptors - descriptors like relative maximum gain
%   - LIK
%
% supersaturation-modified NR fxn, baseline subtracted out:
%   R(c) = Rm*c^n/c50^(s*n)+c^(s*n)
%   Rm = r100 * (1+0.5)
%    % This assumes R(1) == r100 if c50 is 0.5 and n, s == 1,
%    %   and prevents c50 and Rm from playing off each other
%
% The prior probabilities of the values of the parameters are assumed
% to be uniform over the values called 
%
% See also:
%   vis.bayes.noise.proportional()

%% DEFINE VARIABLES
% from data_struct
resp = resp_struct.mean_responses(:);
c_value = resp_struct.contrast(:);
num_trials = resp_struct.num_trials(:);

% from grid_size
r100_values = param_grid.r100;
c50_values = param_grid.c50;
n_values = param_grid.n;
s_values = param_grid.s;

% from noise_mdl
offset = noise_mdl(1);
slope = noise_mdl(2);
prop_mdl = [offset slope];

%% BUILD BAYES GRID MATRIX
Lik = NaN(numel(param_grid.r100),numel(param_grid.c50),numel(param_grid.n),numel(param_grid.s));

%% descriptors, same dimensions as Lik
RelativeMaximumGain = Lik;
SaturationIndex = Lik;
ContrastThreshold = Lik;

%% POSTERIOR PROBABILITY FOR 4 PARAMETERS 
total_steps = (length(param_grid.r100)*length(param_grid.c50)*length(param_grid.n)*length(param_grid.s));
current_step = 0;
step_landmarks = round(linspace(1,total_steps,40)); 
next_landmark = 1;

c_high_res = [0:0.01:1]';
noise_trials = num_trials(1) * ones(size(c_high_res));


for r100 = 1:length(param_grid.r100)
	for c50 = 1:length(param_grid.c50)
		for n = 1:length(param_grid.n)
			for s = 1:length(param_grid.s)
				cresp = r100_values(r100)*(1+0.5) * ...
					vis.contrast.naka_rushton_func(c_value,...
					c50_values(c50), n_values(n), s_values(s));

				% probability density of contrast response
				presp = normpdf(cresp(:)-resp,zeros(size(cresp(:))),...
					vis.bayes.noise.proportional(prop_mdl,cresp,num_trials));
				multipresp = squeeze(prod(presp));
				Lik(r100,c50,n,s) = multipresp;

				nr_high_res = r100_values(r100)*(1+0.5) * ...
					vis.contrast.naka_rushton_func(c_high_res, ...
					c50_values(c50), n_values(n), s_values(s));
				RelativeMaximumGain(r100,c50,n,s) = vis.contrast.indexes.contrastfit2relativemaximumgain(...
					c_high_res, nr_high_res);
				SaturationIndex(r100,c50,n,s) = (max(nr_high_res)-nr_high_res(end))/nr_high_res(end);
				noise_high_res = vis.bayes.noise.proportional(prop_mdl,nr_high_res,noise_trials);
				index = find(nr_high_res>=noise_high_res*2,1,'first');
				cthresh = Inf;
				if ~isempty(index), 
					cthresh = c_high_res(index);
				end;
				ContrastThreshold(r100,c50,n,s) = cthresh;

				current_step = current_step + 1;
				if (current_step == step_landmarks(next_landmark))
					%disp(['Progress: ' int2str(current_step) ' step(s) of ' int2str(total_steps) ', or ' num2str(100*current_step/total_steps) '%']);
					next_landmark = next_landmark+1;
				end
			end
		end
	end
end

Lik = Lik./sum(Lik(:));

%% EXTRACT MARGINAL LIKELIHOOD OF EACH PARAMETER 
lik_r100 = squeeze(sum(sum(sum(Lik,4),3),2));
lik_r100 = lik_r100./sum(lik_r100,"all");

lik_c50 = squeeze(sum(sum(sum(Lik,4),3),1))';
lik_c50 = lik_c50./sum(lik_c50,"all");

lik_n = squeeze(sum(sum(sum(Lik,4),2),1));
lik_n = lik_n./sum(lik_n,"all");

lik_s = squeeze(sum(sum(sum(Lik,3),2),1));
lik_s = lik_s./sum(lik_s,"all");

%% MAXIMUM LIKELIHOOD CONDITIONS
[~,ind] = max(Lik,[],'all');
ind = squeeze(ind);
[r100,c50,n,s] = ind2sub(size(Lik),ind);
cresp = (param_grid.r100(r100))*(1+0.5) * ...
	vis.contrast.naka_rushton_func(c_value, c50_values(c50), n_values(n), s_values(s));
SAT = SaturationIndex(r100,c50,n,s);
CTHRESH = ContrastThreshold(r100,c50,n,s);
RMG = RelativeMaximumGain(r100,c50,n,s);

%% DESCRIPTOR estimation

sat_bins = 0:0.05:10;
cthresh_bins = 0:0.01:1;
rmg_bins = 0:0.1:10;

sat_lik = 0 * sat_bins;
cthresh_lik = 0 * cthresh_bins;
rmg_lik = 0 * rmg_bins;

[sat_index] = discretize(SaturationIndex(:),sat_bins);
for bin = 1:max(sat_index)
    F = sat_index==bin;
    sat_lik(bin) = sum(Lik(F));
end
sat_lik = sat_lik./sum(sat_lik);

[cthresh_index] = discretize(ContrastThreshold(:),cthresh_bins);
for bin = 1:max(cthresh_index)
    F = cthresh_index==bin;
    cthresh_lik(bin) = sum(Lik(F));
end
cthresh_lik = cthresh_lik./sum(cthresh_lik);

[rmg_index] = discretize(RelativeMaximumGain(:),rmg_bins);
for bin = 1:max(rmg_index)
    F = rmg_index==bin;
    rmg_lik(bin) = sum(Lik(F));
end
rmg_lik = rmg_lik./sum(rmg_lik);


%% CREATE OUTPUT_STRUCT
output_struct = struct( ...
	'noise_model',struct('type',{'proportional'},'offset',offset,'slope',slope), ...
	'other_parameters',struct('independent_variable',{'contrast'},'independent_variable_value',resp_struct.contrast,...
		'mean_responses',resp_struct.mean_responses,'num_trials',resp_struct.num_trials), ...
	'marginal_likelihoods',struct( ...
		'r100',struct('values',param_grid.r100,'likelihoods',lik_r100), ...
		'c50',struct('values',param_grid.c50,'likelihoods',lik_c50), ...
		'n',struct('values',param_grid.n,'likelihoods',lik_n), ...
		's',struct('values',param_grid.s,'likelihoods',lik_s)), ... 
	'maximum_likelihood_parameters',struct(...
			'parameters',struct('r100',param_grid.r100(r100),'c50',param_grid.c50(c50),...
				'n',param_grid.n(n),'s',param_grid.s(s)),...
			'curve',cresp,...
			'descriptors', struct('RelativeMaximumGain',RMG,'ContrastThreshold',CTHRESH, 'SaturationIndex',SAT)),...
	'descriptors', struct(...
		'RelativeMaximumGain', struct('values',rmg_bins,'likelihoods',rmg_lik), ...
		'ContrastThreshold', struct('values',cthresh_bins,'likelihoods',cthresh_lik), ...
		'SaturationIndex', struct('values',sat_bins,'likelihoods',sat_lik) ) );


