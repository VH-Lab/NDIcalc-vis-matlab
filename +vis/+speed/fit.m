function [P,sse] = fit(SF,TF,R, min_xi, max_xi)
% FIT Fit two-dimensional Gaussian function to data set
%
%  [P,SSE] = FIT(SF,TF,R, min_xi)
%
%  Fit a Gaussian to a set of responses generated by a speed tuning curve.
%  SSE is the sum of squared errors between the fit and the data.
%
%  Inputs:
%    SF is an array of spatial frequency values
%    TF is an array of temporal frequency values
%    R is an array of measured responses driven by the spatial and
%    temporal frequency values
%    MIN_XI: if provided, provides the lower limit on XI. If not provided,
%    assumed to be 0. (Might provide -1, for example.)
%
%  Outputs:
%    P is an array with parameters:
%    ---------------------------------------------
%    | A         | Peak response of the neuron   |
%    |           |                               |
%    | zeta      | Skew of the temporal          |
%    |           | frequency tuning curve        |
%    |           |                               |
%    | xi        | Speed parameter               |
%    |           |                               |
%    | sigma_sf  | Tuning width of the neuron    |
%    |           | for spatial frequency         |
%    |           |                               |
%    | sigma_tf  | Tuning width of the neuron    |
%    |           | for temporal frequency        |
%    |           |                               |
%    | sf0       | Preferred spatial frequency   |
%    |           | averaged across temporal      |
%    |           | frequencies                   |
%    |           |                               |
%    | tf0       | Preferred temporal frequency  |
%    |           | averaged across spatial       |
%    |           | frequencies                   |
%    ---------------------------------------------
%
%  See: Priebe et al. 2006
%
%  By Noah Lasky-Nielson
%

if nargin<4,
    min_xi = 0;
end;

if nargin<5,
    max_xi = 1;
end;

SF = SF(:); TF = TF(:); R = R(:); % Make them column vectors

% Define constraints on the fit
Upper = [2*max(R); 3; max_xi; 4; 4; max(SF); max(TF)];
Lower = [0; -3; min_xi; 0.1; 0.1; min(SF); min(TF)];

zeta_rand = (Upper(2) - Lower(2)).*rand(10,1) + Lower(2); % Generate an array of 10 random values for fitting
xi_rand = (Upper(3) - Lower(3)).*rand(10,1) + Lower(3); % Generate an array of 10 random values for fitting
sigma_sf_rand = (Upper(4) - Lower(4)).*rand(10,1) + Lower(4); % Generate an array of 10 random values for fitting
sigma_tf_rand = (Upper(5) - Lower(5)).*rand(10,1) + Lower(5); % Generate an array of 10 random values for fitting

 % StartPoint = [2; 1; 0.6; 3.2; 1.4; 2.8; 3.1]; % Arbitrarily chosen  % replace with more random selections

StartPoint = [[0.2*max(R); zeta_rand(1,:); xi_rand(1,:); sigma_sf_rand(1,:); sigma_tf_rand(1,:); 0.1*max(SF); 0.1*max(TF)]...
	[0.4*max(R); zeta_rand(2,:); xi_rand(2,:); sigma_sf_rand(2,:); sigma_tf_rand(2,:); 0.2*max(SF); 0.2*max(TF)]...
	[0.6*max(R); zeta_rand(3,:); xi_rand(3,:); sigma_sf_rand(3,:); sigma_tf_rand(3,:); 0.3*max(SF); 0.3*max(TF)]...
	[0.8*max(R); zeta_rand(4,:); xi_rand(4,:); sigma_sf_rand(4,:); sigma_tf_rand(4,:); 0.4*max(SF); 0.4*max(TF)]...
	[max(R); zeta_rand(5,:); xi_rand(5,:); sigma_sf_rand(5,:); sigma_tf_rand(5,:); 0.5*max(SF); 0.5*max(TF)]...
	[1.2*max(R); zeta_rand(6,:); xi_rand(6,:); sigma_sf_rand(6,:); sigma_tf_rand(6,:); 0.6*max(SF); 0.6*max(TF)]...
	[1.4*max(R); zeta_rand(7,:); xi_rand(7,:); sigma_sf_rand(7,:); sigma_tf_rand(7,:); 0.7*max(SF); 0.7*max(TF)]...
	[1.6*max(R); zeta_rand(8,:); xi_rand(8,:); sigma_sf_rand(8,:); sigma_tf_rand(8,:); 0.8*max(SF); 0.8*max(TF)]...
	[1.8*max(R); zeta_rand(9,:); xi_rand(9,:); sigma_sf_rand(9,:); sigma_tf_rand(9,:); 0.9*max(SF); 0.9*max(TF)]...
	[2*max(R); zeta_rand(10,:); xi_rand(10,:); sigma_sf_rand(10,:); sigma_tf_rand(10,:); max(SF); max(TF)]]; % Some are abritrarily chosen and others randomly generated


% Use lsqcurvefit
optimal_error = Inf; % The initial optimal error value, will be reduced later
P_optimal = []; % Optimal P based on the various supplied starting points
opts = optimset('Display','off');
for i=1:size(StartPoint,2), % Loop through the start point values
	[P,error] = lsqcurvefit(@(P,SF) vis.speed.tuningfunc(SF,TF,P),StartPoint(:,i),SF,R,Lower,Upper,opts);
	if (error < optimal_error) % If the error from the fit is lower than the current optimal error
		optimal_error = error; % Set the optimal error to the new error
		P_optimal = P; % Set the optimal P to be the fitted P with the lower error
	end
end

P = P_optimal; % Once the loop has finished, set P to be the most optimal fit
sse = optimal_error;

