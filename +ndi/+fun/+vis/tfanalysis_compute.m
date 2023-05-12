function [assoc]=tfanalysis_compute(resp)

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
%  'TF sig resp p'            |   p value of ANOVA across conditions
%  'TF sig resp'              |   (0/1) Is above < 0.05?
%  'TF visual response p'     |   p value of ANOVA across conditions + 
%                             |     blank, if available
%  'TF visual response'       |   (0/1) Is above < 0.05?
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

assoclist = { 'TF Low TF','TF High TF',...
	'TF Response curve','TF Pref','TF Low','TF High',...
	'TF sig resp p','TF sig resp','TF visual response p','TF visual response','TF Blank Response',...
	'TF DOG params','TF DOG Fit','TF DOG R2','TF DOG Low','TF DOG High','TF DOG Pref',...
	'TF spline Fit','TF spline Pref','TF spline Low','TF spline High'};

assoc=struct('type','t','owner','t','data',0,'desc',0); assoc=assoc([]);

if nargin==0, assoc = assoclist; return; end;

r0=resp.spont(1);
if isfield(resp,'blankresp'),
	r0 = resp.blankresp(1);
end;

% from the raw curve
[mf,preftf]=max(resp.curve(2,:));
preftf = resp.curve(1,preftf(1));
[lowv, maxv, highv] = compute_halfwidth(resp.curve(1,:),resp.curve(2,:));
[lowvsp, maxvsp, highvsp] = compute_halfwidth(resp.curve(1,:),resp.curve(2,:)-r0);

 % significance
 
[tfsigp,tfvp] = neural_response_significance(resp);

 % DOG fit
re = mf; ri = mf; se = maxv; si = 5+5*randn;
rcurve = resp.curve;

search_options=optimset('fminsearch');
search_options.TolFun=1e-3;
search_options.TolX=1e-3;
%search_options.MaxFunEvals='300*numberOfVariables';
search_options.Display='off';

norm_error_overall = Inf;
dog_par_overall = [];

for jj=1:10,
	re = mf; ri = mf; se = maxv; si = 5+5*randn;
	dog_par=fminsearch('dog_error',[r0 re se ri si],search_options,...
		[rcurve(1,:) 60 70],[rcurve(2,:) r0 r0], ...
		[rcurve(4,:) mean(rcurve(4,:)) mean(rcurve(4,:))] 	    )';

	norm_error=dog_error(dog_par, [rcurve(1,:) 60 70],[rcurve(2,:) r0 r0]);

	if norm_error<norm_error_overall,
		dog_par_overall = dog_par;
	end;
end;

norm_error = norm_error_overall;
dog_par = dog_par_overall;

tfrange_interp=logspace( log10(min( min(rcurve(1,:)),0.01)),log10(50),50);
if isempty(dog_par),
	norm_error = Inf;
	r2 = -Inf;
	response=NaN*tfrange_interp;
else,
	norm_error=dog_error(dog_par, [rcurve(1,:) ],[rcurve(2,:)]);
	r2 = norm_error - ((rcurve(2,:)-mean(rcurve(2,:)))*(rcurve(2,:)'-mean(rcurve(2,:))));
	response=dog(dog_par',tfrange_interp);
end;

	
[lowdog, prefdog, highdog] =   compute_halfwidth(tfrange_interp,response);

fitx = min(rcurve(1,:)):1:max(rcurve(1,:));
if fitx(end)~=max(rcurve(1,:)), fitx(end+1) = max(rcurve(1,:)); end;
fity = interp1([rcurve(1,:)],[rcurve(2,:)], fitx,'spline');
[lowspline, prefspline, highspline] = compute_halfwidth(fitx,fity);

assoc(end+1) = tfassoc('TF Response curve',rcurve,'TF Response curve');
if isfield(resp,'blankresp'),
        assoc(end+1)=tfassoc('TF Blank Response',resp.blankresp,'TF blank resp');
end;
assoc(end+1) = tfassoc('TF Pref',preftf,'TF Pref');
assoc(end+1) = tfassoc('TF Low',lowv,'TF Low cut-off (half-max)');
assoc(end+1) = tfassoc('TF High',highv,'TF High cut-off (half-max)');

assoc(end+1) = tfassoc('TF Low TF',lowvsp,'TF Low cut-off (half-max) TF');
assoc(end+1) = tfassoc('TF High TF',highvsp,'TF High cut-off (half-max) TF');

assoc(end+1) = tfassoc('TF sig resp p',tfsigp,'TF response p value');
assoc(end+1) = tfassoc('TF sig resp',tfsigp<=0.05,'Is TF response significant?');
assoc(end+1) = tfassoc('TF visual response p',tfvp,'TF visual response (including blank) p value');
assoc(end+1) = tfassoc('TF visual response',tfvp<=0.05,'Is TF response across stims and blank significant?');

assoc(end+1) = tfassoc('TF DOG Pref',prefdog,'TF Pref difference of gaussians');
assoc(end+1) = tfassoc('TF DOG Low',lowdog,'TF Low cut-off (half-max from DOG)');
assoc(end+1) = tfassoc('TF DOG High',highdog,'TF High cut-off (half-max from DOG)');
assoc(end+1) = tfassoc('TF DOG R2',r2,'TF DOG r^2 of fit');
assoc(end+1) = tfassoc('TF DOG Fit',[tfrange_interp; response],'TF DOG Fit, 1st row is TF, 2nd row is response');

assoc(end+1) = tfassoc('TF spline Pref',prefspline,'TF Pref spline');
assoc(end+1) = tfassoc('TF spline Low',lowspline,'TF Low cut-off (half-max from spline)');
assoc(end+1) = tfassoc('TF spline High',highspline,'TF High cut-off (half-max from spline)');
assoc(end+1) = tfassoc('TF spline Fit',[fitx; fity],'TF spline Fit, 1st row is TF, 2nd row is response');

function myassoc = tfassoc(type,data,desc)
myassoc  = struct('type',type,'owner','tf','data',data,'desc',desc);
