function r2 = rsquared(y, yfit)
% RSQUARED - coefficient of determination for a frequency-tuning fit
%
%  R2 = vis.frequency.rsquared(Y, YFIT)
%
%  Returns the coefficient of determination
%
%     R2 = 1 - mean((YFIT-Y).^2) / mean((Y-mean(Y)).^2)
%
%  where Y is the vector of measured responses and YFIT is the model
%  evaluated at the same frequencies at which Y was measured. This is the
%  same expression used by vis.frequency.movshon2005_fit, so the R2 values
%  reported by every fit in vis.frequency.spatial_frequency_analysis and
%  vis.frequency.temporal_frequency_analysis are directly comparable.
%
%  R2 is NaN if the fit failed (YFIT empty or not the same length as Y) or
%  if the responses have no variance to explain (mean((Y-mean(Y)).^2)==0),
%  in which case the ratio is undefined rather than perfect.
%
%  Note that the cubic spline of STEP 5 interpolates the measured responses
%  rather than fitting them, so it passes exactly through every data point
%  and its R2 is identically 1. That value is reported for completeness and
%  for a uniform document structure; it is not evidence that the spline
%  describes the data better than the fitted models, and it should not be
%  compared against them.
%
%  See also: vis.frequency.movshon2005_fit
%

y = y(:);
yfit = yfit(:);

if isempty(yfit) | numel(yfit)~=numel(y),
	r2 = NaN;
	return;
end;

sst = mean( (y-mean(y)).^2 );

if sst==0,
	r2 = NaN;
	return;
end;

r2 = 1 - mean( (yfit-y).^2 ) / sst;
