function s = randomstream(seed)
% RANDOMSTREAM - a reproducible random number stream
%
%  S = vis.randomstream(SEED)
%
%  Returns a random number stream to be passed to RAND, RANDN and friends
%  (for example rand(S) or randn(S,3,4)), so that a calculation which depends
%  on random numbers can be repeated exactly.
%
%  SEED may be:
%    a non-negative integer | a private stream with that seed. Repeated calls
%                           |   with the same seed draw the same numbers.
%    'shuffle'              | a private stream seeded from the clock. Use this
%                           |   to deliberately vary a calculation run to run.
%    [] (empty)             | the global stream, that is, no seeding at all.
%                           |   Numbers drawn will not be reproducible.
%
%  A private stream leaves the caller's global stream untouched, so seeding a
%  calculation does not disturb any other random numbers the caller draws.
%
%  Example:
%     s = vis.randomstream(1);
%     x = rand(s,1,10);   % the same ten numbers on every run
%
%  See also: RandStream, vis.speed.fit
%

if nargin<1,
	seed = 0;
end;

if isempty(seed),
	s = RandStream.getGlobalStream();
else,
	s = RandStream('twister','Seed',seed);
end;
