function [tf_props]=temporal_frequency_analysis(resp)

%  TFANALYSIS_COMPUTE  Analyze reponses to temporal frequencies
% %  [ASSOC]=TFANALYSIS_COMPUTE(RESP)
%
%  Analyzes temporal frequency responses.
%
%  RESP is a structure list of response properties with fields:
%  curve    |    4xnumber of temporal frequencies tested,
%           |      curve(1,:) is temporal frequencies tested
%           |      curve(2,:) is mean responses
%           |      curve(3,:) is standard deviation
%           |      curve(4,:) is standard error
%  ind      |    cell list of individual trial responses for each TF
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
%  'TF Response curve'        |   Response curve (tfs;mean;stddev;stderr)
%  'TF Pref'                  |   TF w/ max response
%  'TF Low'                   |   low TF with half of max response 
%  'TF High'                  |   high TF with half of max response
%
%        Same as above with 'blank' or 'spont' rate subtracted
%  'TF Low TF'                |   low TF with half of max response 
%  'TF High TF'               |   high TF with half of max response
%
%  Difference of gaussians fit:
%  'TF DOG params'            |   'r0 re se ri si'
%  'TF DOG Fit'               |   1st row has TF values, 2nd has responses
%  'TF DOG R2'                |   R^2 error
%  'TF DOG Low'               |   Low cut-off, as measured with DOG
%  'TF DOG High'              |   High cut-off, as measured with DOG
%  'TF DOG Pref'              |   TF Pref, as measured with DOG
%
%  Cubic spline "Fit":
%  'TF spline Fit'            |   1st row has TF values, 2nd has responses
%  'TF spline Pref'           |   TF Pref, as measured with spline
%  'TF spline Low'            |   Low cut-off, as measured with spline
%  'TF spline High'           |   High cut-off, as measured with spline
%


  % STEP 1: empirical parameters (fitless)

[mf,preftf]=max(resp.curve(2,:));
preftf = resp.curve(1,preftf(1));
[lowv, maxv, highv] = compute_halfwidth(resp.curve(1,:),resp.curve(2,:));

[tf_props.fitless.L50,tf_props.fitless.Pref,tf_props.fitless.H50] = ...
	ndi.fun.vis.compute_halfwidth(resp.curve(1,:),resp.curve(2,:));

 % STEP 2: DOG fit
rcurve = resp.curve;

search_options=optimset('fminsearch');
search_options.TolFun=1e-3;
search_options.TolX=1e-3;
%search_options.MaxFunEvals='300*numberOfVariables';
search_options.Display='off';

norm_error_overall = Inf;
dog_par_overall = [];

for jj=1:10,
	r0 = 0; re = mf; ri = mf; se = maxv; si = 5+5*randn;
	dog_par=fminsearch('dog_error',[r0 re se ri si],search_options,...
		[rcurve(1,:) 100 120],[rcurve(2,:) 0 0], ...
		[rcurve(4,:) mean(rcurve(4,:)) mean(rcurve(4,:))] 	    )';

	norm_error=dog_error(dog_par, [rcurve(1,:) 100 120],[rcurve(2,:) 0 0 ]);

	if norm_error<norm_error_overall,
		dog_par_overall = dog_par;
	end;
end;

norm_error = norm_error_overall;
dog_par = dog_par_overall;

tfrange_interp=logspace( log10(min( min(rcurve(1,:)),0.01)),log10(120),100);
if isempty(dog_par),
	norm_error = Inf;
	r2 = -Inf;
	response=NaN*tfrange_interp;
else,
	norm_error=dog_error(dog_par, [0 rcurve(1,:) 100 120 ],[0 rcurve(2,:) 0 0]);
	r2 = norm_error - ((rcurve(2,:)-mean(rcurve(2,:)))*(rcurve(2,:)'-mean(rcurve(2,:))));
	response=dog(dog_par',tfrange_interp);
end;

	
[lowdog, prefdog, highdog] = ndi.fun.vis.compute_halfwidth(tfrange_interp,response);

fit_dog.parameters = dog_par;
fit_dog.values = tfrange_interp;
fit_dog.fit = response;
fit_dog.L50 = lowdog;
fit_dog.Pref = prefdog;
fit_dog.H50 = highdog;

tf_props.fit_dog = fit_dog;

 % STEP 3: spline fitting

fitx = tfrange_interp;
fity = interp1([0 rcurve(1,:) 100 120],[0 rcurve(2,:) 0 0], fitx,'spline');
[lowspline, prefspline, highspline] = ndi.fun.vis.compute_halfwidth(fitx,fity);

fit_spline.values = fitx;
fit_spline.fit = fity;
fit_spline.L50 = lowspline;
fit_spline.Pref = prefspline;
fit_spline.H50 = highspline;

tf_props.fit_spline = fit_spline;

% STEP 4: gausslog fitting

a = 0;
a_range = [0 0];
b = mf;
b_range = [0 2*max(mf,0)];
c = maxv;
d = rand*(highv - lowv);
e = 0;
e_range = [ 0 0 ];

[gausslog_par,gof,gausslog_fitcurve] = vlt.fit.gausslogfit([0 rcurve(1,:) 100 120]',[0 rcurve(2,:) 0 0]',...
	'a_hint',a,'a_range',a_range,'b_hint',b,'b_range',b_range,'c_hint',c,'d_hint',d,'e_hint',e,'e_range',e_range);

[low_gausslog,pref_gausslog,high_gausslog] = ndi.fun.vis.compute_halfwidth(gausslog_fitcurve(1,:),gausslog_fitcurve(2,:));

fit_gausslog.parameters = gausslog_par;
fit_gausslog.values = tfrange_interp;
fit_gausslog.fit = vlt.math.gausslog(fit_gausslog.values,fit_gausslog.parameters);
fit_gausslog.L50 = low_gausslog;
fit_gausslog.Pref = pref_gausslog;
fit_gausslog.H50 = high_gausslog;

tf_props.fit_gausslog = fit_gausslog;

