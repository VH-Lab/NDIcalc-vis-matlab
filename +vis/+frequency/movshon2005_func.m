function r = movshon2005_func(f, P)
% MOVSHON2005_FUNC - Movshen et al. 2005 frequency response function
%
% R = MOVSHON2005_FUNC(F, P)
%
% Computes responses to frequencies according to Movshon et al. 2005,
% J Neurosci 25:2712-2722.
%
% F is a set of frequencies over which to make the calculation.
% P is a vector with the parameters where
%   k - P(1) - scaling factor
%   fc - P(2) - characteristic temporal frequency
%   fh - P(3) - corner frequency of the low-frequency limb
%   B - P(4) - slope of the low-frequency limb
%
% If P has 5 parameters, then there is a constant term C. Otherwise C is 0.
%   C - P(5)
% 
% The function has the form:
%
%    R(f) = k * exp(-(f./fc).^2) ./ (1+(fh./f).^B) + C
%
%
% Example:
%   f = logspace(log10(0.01),log10(32),50);
%   k = 10;
%   fc = 5;
%   fh = 12;
%   B = 1;
%   R = vis.frequency.movshon2005_func(f,[k fc fh B]);
%
%   figure;
%   plot(f,R,'ko-');
%   set(gca,'xscale','log');
%   box off;
%   ylabel('Response');
%   xlabel('Frequency');
%
%  

k = P(1);
fc = P(2);
fh = P(3);
B = P(4);
if numel(P)>4,
	C = P(5);
else,
	C = 0;
end;

r = k * exp(-(f./fc).^2) ./ (1+(fh./f).^B) + C;

r = r(:);

