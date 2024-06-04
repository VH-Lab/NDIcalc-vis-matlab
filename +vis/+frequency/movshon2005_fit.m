function [p,rfit,mse,rsquared] = movshon2005_fit(f, r)
% MOVSHON2005_FIT - fit a Movshen et al. 2005 frequency response function
%
% [P,RFIT,MSE,RSQUARED] = MOVSHON2005_FUNC(F, R)
%
% Computes a least squares fit of responses (R) to frequencies (F) according to
% Movshon et al. 2005, J Neurosci 25:2712-2722.
%
% F is a vector of frequencies over which to make the calculation.
% R is a corresponding vector of responses for each value of F
% Outpus:
% P is a vector with the parameters where
%   k - P(1) - scaling factor
%   fc - P(2) - characteristic temporal frequency
%   fh - P(3) - corner frequency of the low-frequency limb
%   B - P(4) - slope of the low-frequency limb
% MSE - mean squared error (per point)
% RSQUARED - the RSQUARED value of the fit.
% 
% The function has the form:
%
%    R(f) = k * exp(-(f./fc).^2) ./ (1+(fh./f).^B)
%
% See MOVSHON2005_FUNC (vis.frequency.movshon2005_func)
%
% Example:
%   [cell1,cell2] = vis.frequency.movshon2005_cells();
%   [P1,rfit1,mse1,rsq1] = vis.frequency.movshon2005_fit(cell1(:,1),cell1(:,2));
%   [P2,rfit2,mse2,rsq2] = vis.frequency.movshon2005_fit(cell2(:,1),cell2(:,2));
%   freq_points = logspace(log10(0.02),log10(64),50);
%
%   figure;
%      % the data points
%   plot(cell1(:,1),cell1(:,2),'ro'); 
%   hold on;
%      % the fit
%   plot(freq_points,vis.frequency.movshon2005_func(freq_points,P1),'r-'); 
%      % the second data points
%   plot(cell2(:,1),cell2(:,2),'bo');
%      % the second fit 
%   plot(freq_points,vis.frequency.movshon2005_func(freq_points,P2),'b-'); 
%   box off;
%   ylabel('Response (ips)');
%   xlabel('Temporal Frequency (Hz)');
%   set(gca,'xscale','log');
%

myfun = @(P,x) vis.frequency.movshon2005_func(x,P);

startPoint = [max(r);median(f);median(f);1];
startPoint(:,2) = [max(r); 0; max(f); 1];
startPoint(:,3) = [max(r); max(f); 0; 1];
startPoint(:,4) = [max(r); median(f); median(f)/2; 1];
startPoint(:,5) = [max(r); median(f); 2*median(f); 1];

lowerB = [0; min(f)/10; min(f)/10; 1/10000];
upperB = [2 * max(r); max(f)*10; max(f)*10; 10000];

best_err = Inf;
P_best = startPoint(:,1);

opts = optimset('Display','off');

for i=1:size(startPoint,1),
	[P_here,err_here] = lsqcurvefit(myfun, startPoint(:,i), f, r, lowerB, upperB,opts);
	if err_here < best_err,
		best_err = err_here;
		P_best = P_here;
	end;
end;

rfit = vis.frequency.movshon2005_func(f,P_best);

p = P_best;

mse = mean( (rfit-r).^2);

rsquared = 1 - mse/(mean( (r-mean(r)).^2  ));
