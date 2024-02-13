function [params,err] = dog_fit(x,y,s,varargin)
% DOG_FIT - perform a difference of gaussians fit
%
% [PARAMS, ERR] = DOG_FIT(X,Y,S, ...)
%
% Perform a difference of gaussians fit Y_ for the equation
%
% Y_ = a1*exp(-X.^2/(2*b1^2)) - a2*exp(-X.^2/(2*b2^2))
%
% X and Y are column vectors with values of the fit. If S is provided,
% it is expected to be the standard deviation of the measurement and
% the fit will be weighted by 1/(1+s) for each entry (so more variable
% points get less weight).
%
% PARAMS are the parameters [A1 B1 A2 B2] of the DOG fit (see help DOG).
% ERR is the averaged squared error over all entries of Y, weighted
% by the weights calculated with S if provided.
%
% This function takes name/value pairs that influence its behavior
% |-----------------------------------------------------------------|
% | Parameter (default)  | Description                              |
% |----------------------|------------------------------------------|
% | start_positions (10) | Number of random start positions to use  |
% |_---------------------|------------------------------------------|
%
%

if nargin<3,
	s = [];
end;

if isempty(s),
	w = 1+0*y(:);,
else,
	w = 1./(1+s(:));
end;

start_positions = 10;

did.datastructures.assign(varargin{:});

mydog = fittype('a*exp(-(x.^2)/(2*b^2)) - c*exp(-(x.^2)/(2*d^2))');

my = max(y(:));
mx = max(x(:));

fo = fitoptions(mydog);
fo.Lower = [0; eps; 0; eps];
fo.Upper = [10*my; mx; 10*my; mx;];

best_error = Inf;
dog_param_best = [];

warn_state = warning('off');

for jj=1:start_positions,
	mydog = setoptions(mydog,fo);

	[mydog_fit,mydog_gof] = fit(x(:),y(:),mydog,'weight',w);

	if mean(mydog_gof.sse)<best_error,
		best_error = mean(mydog_gof.sse);
		dog_param_best(1) = mydog_fit.a;
		dog_param_best(2) = mydog_fit.b;
		dog_param_best(3) = mydog_fit.c;
		dog_param_best(4) = mydog_fit.d;
	end;
end;

params = dog_param_best;
err = best_error; 

warning(warn_state);
