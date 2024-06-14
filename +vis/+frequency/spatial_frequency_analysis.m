function [sf_props]=spatial_frequency_analysis(resp)

%  SFANALYSIS_COMPUTE  Analyze reponses to spatial frequencies
% %  [ASSOC]=SFANALYSIS_COMPUTE(RESP)
%
%  Analyzes spatial frequency responses.
%
%  RESP is a structure list of response properties with fields:
%  curve    |    4xnumber of spatial frequencies tested,
%           |      curve(1,:) is spatial frequencies tested
%           |      curve(2,:) is mean responses
%           |      curve(3,:) is standard deviation
%           |      curve(4,:) is standard error
%  ind      |    cell list of individual trial responses for each SF
%  spont    |    spontaneous responses [mean stddev stderr]
%  spontind |    individual spontaneous responses
%  Optionally:
%  blank    |    response to a blank trial: [mean stddev stderr]
%  blankind |    individual responses to blank
%   
%  
%  If the function is called with no arguments, then the of associate
%  names that are computed by the function is returned.
%
%  Returns data in the form of 'associates' that can be added
%  to a measured data object:
%  'SF Response curve'        |   Response curve (tfs;mean;stddev;stderr)
%  'SF Pref'                  |   SF w/ max response
%  'SF Low'                   |   low SF with half of max response 
%  'SF High'                  |   high SF with half of max response
%
%        Same as above with 'blank' or 'spont' rate subtracted
%  'SF Low SF'                |   low SF with half of max response 
%  'SF High SF'               |   high SF with half of max response
%
%  Difference of gaussians fit:
%  'SF DOG params'            |   'r0 re se ri si'
%  'SF DOG Fit'               |   1st row has SF values, 2nd has responses
%  'SF DOG R2'                |   R^2 error
%  'SF DOG Low'               |   Low cut-off, as measured with DOG
%  'SF DOG High'              |   High cut-off, as measured with DOG
%  'SF DOG Pref'              |   SF Pref, as measured with DOG
%
%  Cubic spline "Fit":
%  'SF spline Fit'            |   1st row has SF values, 2nd has responses
%  'SF spline Pref'           |   SF Pref, as measured with spline
%  'SF spline Low'            |   Low cut-off, as measured with spline
%  'SF spline High'           |   High cut-off, as measured with spline
%
%  Movshon et al 2005 fit: 
%  'SF MV params'             |   [k fc fh B]
%  'SF MV Fit'                |   First row has SF values, 2nd has responses
%  'SF MV R2'                 |   R2 error.
%  'SF MV L50'                |   Low cut-off
%  'SF MV Pref'               |   SF Preference
%  'SF MV H50'                |   High cut-off
%


  % STEP 1: empirical parameters (fitless)

[mf,preftf]=max(resp.curve(2,:));
preftf = resp.curve(1,preftf(1));
[lowv, maxv, highv] = vis.compute_halfwidth_interp(resp.curve(1,:),resp.curve(2,:));

[sf_props.fitless.L50,sf_props.fitless.Pref,sf_props.fitless.H50] = ...
	vis.compute_halfwidth_interp(resp.curve(1,:),resp.curve(2,:));

sf_props.fitless.bandwidth = vis.frequency.bandwidth(sf_props.fitless.L50, sf_props.fitless.H50);
[minValue,minLoc] = min(resp.curve(2,:));
[maxValue,maxLoc] = max(resp.curve(2,:));
[minFValue,minFLoc] = min(resp.curve(1,:));
[maxFValue,maxFLoc] = max(resp.curve(1,:));

sf_props.fitless.low_pass_index = rectify(resp.curve(2,minFLoc))/rectify(maxValue);
sf_props.fitless.high_pass_index = rectify(resp.curve(2,maxFLoc))/rectify(maxValue);

 % STEP 2: DOG fit
rcurve = resp.curve;

[dog_par,norm_error] = vis.dog_fit(rcurve(1,:),rcurve(2,:),rcurve(3,:));

norm_error_overall = Inf;
dog_par_overall = [];

for jj=1:10,
	r0 = 0; re = mf; ri = mf; se = maxv; si = 5+5*randn;
    %forcing a better fit by manually adding (100,0), and (120,0)
    %dog_par=fminsearch('dog_error',[r0 re se ri si],search_options,...
	%	[rcurve(1,:) 100 120],[rcurve(2,:) 0 0], ...
	%	[rcurve(4,:) mean(rcurve(4,:)) mean(rcurve(4,:))] 	    )';
	 dog_par=fminsearch('dog_error',[r0 re se ri si],search_options,...
	 	[rcurve(1,:)],[rcurve(2,:)], ...
	 	[rcurve(4,:) mean(rcurve(4,:)) mean(rcurve(4,:))] 	    )';
    %forcing a better fit by manually adding (100,0), and (120,0)
	%norm_error=dog_error(dog_par, [rcurve(1,:) 100 120],[rcurve(2,:) 0 0 ]);
    norm_error=dog_error(dog_par, [rcurve(1,:)],[rcurve(2,:)]);

	if norm_error<norm_error_overall,
		dog_par_overall = dog_par;
	end;
end;

norm_error = norm_error_overall;
dog_par = dog_par_overall;
dog_par = [0 dog_par]; % add in the 0 to comply with old dog

sfrange_interp=logspace( log10(min( min(rcurve(1,:)),0.01)),log10(120),100);
if isempty(dog_par),
	norm_error = Inf;
	r2 = -Inf;
	response=NaN*sfrange_interp;
else,
	%forcing a better fit by manually adding (0,0), (100,0), and (120,0)
  norm_error=dog_error(dog_par, [rcurve(1,:)],[rcurve(2,:)]);
	r2 = norm_error - ((rcurve(2,:)-mean(rcurve(2,:)))*(rcurve(2,:)'-mean(rcurve(2,:))));
	response=dog(dog_par',sfrange_interp);
end;

[lowdog, prefdog, highdog] = vis.compute_halfwidth(sfrange_interp,response);

fit_dog.parameters = dog_par;
fit_dog.values = sfrange_interp;
fit_dog.fit = response;
fit_dog.L50 = lowdog;
fit_dog.Pref = prefdog;
fit_dog.H50 = highdog;
fit_dog.bandwidth = vis.frequency.bandwidth(fit_dog.L50, fit_dog.H50);

sf_props.fit_dog = fit_dog;

 % STEP 3: Movshon fit without constant term

[MV_P,mvfit,mv_mse,MV_R2] = vis.frequency.movshon2005_fit(rcurve(1,:)', rcurve(2,:)');

fit_movshon.parameters = MV_P;
fit_movshon.values = sfrange_interp(:);
fit_movshon.fit = vis.frequency.movshon2005_func(sfrange_interp(:),MV_P);
[lowmov, prefmov, highmov] = vis.compute_halfwidth(sfrange_interp,fit_movshon.fit);
fit_movshon.L50 = lowmov;
fit_movshon.Pref = prefmov;
fit_movshon.H50 = highmov;
fit_movshon.R2 = MV_R2;
fit_movshon.bandwidth = vis.frequency.bandwidth(fit_movshon.L50, fit_movshon.H50);

sf_props.fit_movshon = fit_movshon;

 % STEP 4: Movshon fit without constant term

[MV_P_c,mvfit_c,mv_mse_c,MV_R2_c] = vis.frequency.movshon2005_fit(rcurve(1,:)', rcurve(2,:)',1);

fit_movshon_c.parameters = MV_P_c;
fit_movshon_c.values = sfrange_interp(:);
fit_movshon_c.fit = vis.frequency.movshon2005_func(sfrange_interp(:),MV_P_c);
[lowmov_c, prefmov_c, highmov_c] = vis.compute_halfwidth(sfrange_interp,fit_movshon_c.fit);
fit_movshon_c.L50 = lowmov_c;
fit_movshon_c.Pref = prefmov_c;
fit_movshon_c.H50 = highmov_c;
fit_movshon_c.R2 = MV_R2_c;
fit_movshon_c.bandwidth = vis.frequency.bandwidth(fit_movshon_c.L50, fit_movshon_c.H50);

sf_props.fit_movshon_c = fit_movshon_c;

 % STEP 5: spline fitting

fitx = sfrange_interp;
%forcing a better fit by manually adding (0,0), (100,0), and (120,0)
%fity = interp1([0 rcurve(1,:) 100 120],[0 rcurve(2,:) 0 0], fitx,'spline');
fity = interp1([rcurve(1,:)],[rcurve(2,:) ], fitx,'spline');
[lowspline, prefspline, highspline] = vis.compute_halfwidth(fitx,fity);

fit_spline.values = fitx;
fit_spline.fit = fity;
fit_spline.L50 = lowspline;
fit_spline.Pref = prefspline;
fit_spline.H50 = highspline;
fit_spline.bandwidth = vis.frequency.bandwidth(fit_spline.L50, fit_spline.H50);

sf_props.fit_spline = fit_spline;

% STEP 6: gausslog fitting

a = 0;
a_range = [0 0];
b = mf;
b_range = [0 2*max(0,mf)];
c = maxv;
d = rand*(highv - lowv);
if isnan(d)
    d = 1;
    %d = randn;
end
e = 0;
e_range = [ 0 0 ];


[gausslog_par,gof,gausslog_fitcurve] = vlt.fit.gausslogfit([rcurve(1,:) ]',[rcurve(2,:) ]',...
	'a_hint',a,'a_range',a_range,'b_hint',b,'b_range',b_range,'c_hint',c,'d_hint',d,'e_hint',e,'e_range',e_range);

[low_gausslog,pref_gausslog,high_gausslog] = vis.compute_halfwidth(gausslog_fitcurve(1,:),gausslog_fitcurve(2,:));

fit_gausslog.parameters = gausslog_par;
fit_gausslog.values = sfrange_interp;
fit_gausslog.fit = vlt.math.gausslog(fit_gausslog.values,fit_gausslog.parameters);
fit_gausslog.L50 = low_gausslog;
fit_gausslog.Pref = pref_gausslog;
fit_gausslog.H50 = high_gausslog;
fit_gausslog.bandwidth = vis.frequency.bandwidth(fit_gausslog.L50,fit_gausslog.H50);

sf_props.fit_gausslog = fit_gausslog;


