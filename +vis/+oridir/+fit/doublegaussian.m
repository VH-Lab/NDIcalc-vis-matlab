function [P,fitcurve,msse,R2] = doublegaussian(angles,responses,options)
% DOUBLEGAUSSIAN - perform a double gaussian fit for stimulus direction
%
%  [P, FITCURVE, MSSE, R2] = DOUBLEGAUSSIAN(ANGLES, RESPONSES)
%
% Given ANGLES and RESPONSES, produce a least-squares fit to the function
%
% R=Rsp+Rp*EXP(-(ANGLES-Op)^2)/(2*sigm^2))+Rn*EXP(-(ANGLES-Op+180)^2/(2*sigm^2))
%
% ANGLES should be in degrees.
% 
% where 
%   Rsp = P(1)
%   Rp = P(2)
%   Op = P(3)
%   sigm = P(4)
%   Rn = P(5)
%
% Rsp is constrained to be within -span and + span, where span is max(responses)-min(responses).
% Rp is constrained to be between 0 and 3 times the peak response.
% Rn is constrained to be between 0 and Rp
% Op is between 0 and 2*pi, excluding 2*pi
% sigm must be larger than the median difference between angles divided by 2 and 90.
%
% Example:
%    angles = 0:30:360-30;
%    P = [ -0.5 20 10 55 39];
%    responses = vlt.neuro.vision.oridir.doublegaussianfunc(angles,P);
%    [P_fit,fitcurve,msse,R2] = vis.oridir.fit.doublegaussian(angles,responses);
%    figure;
%    plot(angles,responses,'o');
%    hold on
%    plot(fitcurve(1,:),fitcurve(2,:),'b-');
%    box off;
%
 

arguments
	angles {mustBeVector}
	responses {mustBeVector}
	options.widthHints = [];
end

 % go columnar!
angles = angles(:);
responses = responses(:); 

if isempty(options.widthHints),
	options.widthHints = linspace(median(diff(angles))/2,90,6);
end;

dg = fittype(@(a,b,c,d,e,x) vlt.data.colvec(vlt.neuro.vision.oridir.doublegaussianfunc(x,[a;b;c*b;d;e])));
fo = fitoptions(dg);

[peak,loc] = max(responses);
rmin = min(responses);

span = peak - rmin;

fo.Lower = [  -span;       0;   0;    -4*180;    median(diff(angles))/2];
fo.Upper = [   span;  max(3*peak,0);   1;     4*180;    90 ];

bestErr = Inf;
bestP = [];
msse = Inf;
R2 = Inf;

for w = 1:numel(options.widthHints),
	fo.StartPoint = [0; peak; 0.5; angles(loc); options.widthHints(w)];
	dg = setoptions(dg,fo);
	[orifit,gof] = fit(angles,responses,dg);
	if gof.sse < bestErr,
		bestErr = gof.sse;
		bestP = [ orifit.a; orifit.b; orifit.b*orifit.c; mod(orifit.d,360); orifit.e];
		msse = bestErr / numel(responses);
		R2 = gof.rsquare;
	end;
end;

angles_fitcurve = 0:359;

fitcurve = [ angles_fitcurve; vlt.neuro.vision.oridir.doublegaussianfunc(angles_fitcurve,bestP)];

P = bestP;
