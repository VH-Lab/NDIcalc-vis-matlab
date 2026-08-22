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
% | start_positions (0)  | Number of EXTRA random start positions to try  |
% |                      |   on top of the deterministic grid below. The   |
% |                      |   grid alone reaches the best fit on every mock |
% |                      |   in this repository, so the default is none.   |
% | RandomSeed (0)       | Seed for those extra start positions. The draws |
% |                      |   use a private stream, so the fit is           |
% |                      |   reproducible and the caller's global random   |
% |                      |   stream is left untouched. Pass 'shuffle' to   |
% |                      |   vary the search between runs, or [] to draw   |
% |                      |   from the global stream. See vis.randomstream. |
% |_---------------------|------------------------------------------------|
%
% The start positions are a deterministic grid derived from the data, not
% random draws over the parameter bounds. See the note in the code for why.
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

start_positions = 0;
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

lb = fo.Lower(:).';
ub = fo.Upper(:).';

 % Where the start positions come from, and why they are not random.
 %
 % This loop used to set no StartPoint at all, which left the choice to the
 % Curve Fitting Toolbox: FIT then picks one at random on every call, so the
 % function returned a different answer each time it was run on the same
 % data. That is the defect this addresses, but simply seeding those draws is
 % not enough, because the box they are drawn from is the wrong box. The
 % width bounds are ten times the full frequency range -- for the mocks here,
 % (0, 599.9] against data spanning 0.01 to 60 cyc/deg -- and a Gaussian that
 % wide is flat to within half a percent over the data. Most of the volume
 % being sampled is a plateau where the gradient is nearly zero, so a large
 % minority of draws never leave it: fitting the mock spatial frequency
 % curves from uniform starts over these bounds reaches the good optimum on
 % roughly two calls in three and collapses to a flat fit on the rest.
 % Seeding that would only pin down which outcome you get.
 %
 % So the starts are a grid built from the data instead, over the region a
 % band-pass difference of gaussians must live in: the centre width spans the
 % measured frequencies, the surround is narrower than the centre in
 % frequency (d < b, which is what makes the response fall off at low
 % frequency rather than rise), and the two amplitudes are comparable so that
 % the difference is small near zero. 30 starts, and on every mock in this
 % repository the best of them reaches the same optimum as the best of the
 % old 40 random ones -- every time rather than most of the time.
 %
 % vis.frequency.movshon2005_fit reached the same conclusion for its own fit,
 % where the comment reads 'go full random -- did not work well'.
fpos = x(x>0);
if isempty(fpos),
	 % nothing to build a grid from; fall back to the bounds
	startGrid = repmat(0.5*(lb+ub),30,1);
else,
	b_cand = logspace(log10(min(fpos(:))),log10(max(fpos(:))),5);
	bd_ratio = [2 4 8];      % how much narrower the surround is, in frequency
	c_frac = [0.5 0.9];      % surround amplitude as a fraction of the centre
	startGrid = [];
	for ib=1:numel(b_cand),
		for ir=1:numel(bd_ratio),
			for ic=1:numel(c_frac),
				startGrid(end+1,:) = [my b_cand(ib) c_frac(ic)*my b_cand(ib)/bd_ratio(ir)];
			end;
		end;
	end;
end;

 % anyone who wants the old wide random search can ask for it by name
if start_positions>0,
	startStream = vis.randomstream(RandomSeed);
	startGrid = [startGrid; (ub-lb).*rand(startStream,start_positions,4)+lb];
end;

 % clamped because ub can fall below lb when the data are degenerate
 % (all responses <= 0, or a single distinct frequency)
startGrid = min(repmat(ub,size(startGrid,1),1),max(repmat(lb,size(startGrid,1),1),startGrid));

warn_state = warning('off');

for jj=1:size(startGrid,1),
	fo.StartPoint = startGrid(jj,:);
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
