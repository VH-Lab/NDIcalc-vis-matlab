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
% |-----------------------------------------------------------------------|
% | Parameter (default)  | Description                                    |
% |----------------------|------------------------------------------------|
% | start_positions (40) | Number of random start positions to use        |
% | RandomSeed (0)       | Seed for the random start positions. The draws |
% |                      |   use a private stream, so the fit is          |
% |                      |   reproducible and the caller's global random  |
% |                      |   stream is left untouched. Pass 'shuffle' to  |
% |                      |   vary the search between runs, or [] to draw  |
% |                      |   from the global stream. See vis.randomstream.|
% |_---------------------|------------------------------------------------|
%
% See also: vis.randomstream, vis.frequency.dog
%

if nargin<3,
	s = [];
end;

if isempty(s),
	w = 1+0*y(:);,
else,
	w = 1./(1+s(:));
end;

start_positions = 40;
RandomSeed = 0;

vlt.data.assign(varargin{:});

mydog = fittype('a*exp(-(x.^2)/(2*b^2)) - c*exp(-(x.^2)/(2*d^2))');

my = max(y(:));
mx = max(x(:)) - min(x(:));

fo = fitoptions(mydog);
fo.Lower = [0; eps; 0; eps];
fo.Upper = [10*max(0,my); 10*mx; 10*max(0,my); 10*mx;];

best_error = Inf;
dog_param_best = [];

 % The start positions are drawn here, from a private seeded stream, rather
 % than left to the Curve Fitting Toolbox. With no StartPoint set, FIT picks
 % one at random on every call, so this loop returned a different answer each
 % time it was run on the same data. Drawing them explicitly makes the fit
 % reproducible; the stream is private, so the caller's global stream is
 % untouched. Bounds are the same ones FIT was already searching within.
lb = fo.Lower(:).';
ub = fo.Upper(:).';
startStream = vis.randomstream(RandomSeed);

warn_state = warning('off');

for jj=1:start_positions,
	 % clamped because ub can fall below lb when the data are degenerate
	 % (all responses <= 0, or a single distinct frequency)
	fo.StartPoint = min(ub,max(lb,(ub-lb).*rand(startStream,1,4)+lb));
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
