function [a, mu, C] = gauss3fit(X, Y, T, response, on)
% ERROR: this function is not working
% GAUSSSPOTFIT - Fit a 2d gaussian to data
%
%  [MU,C,AMP,FIT_RESPONSES] = vlt.fit.gaussspotfit(X, Y, T, RESPONSE, ON)
%
%  Fits a 3d gaussian PDF to responses to circle stimulation at different positions.
%
%  INPUTS:
%    X, Y, and T specify the mesh grid of the responses.
%    RESPONSE is an size(X,1) x size(Y,2) x size(T,3) matrix of responses
%    ON is a boolean that indicates response on or response off is going to
%    be calculated.
%  OUTPUTS:
%    MU - The mean of the best-fit gaussian PDF, in X and Y
%    C  - The covariance matrix of the best-fit gaussian PDF
%    AMP - The amplitude of the response for each circle size.
%    FIT_RESPONSES - the fit responses

error("ERROR: this function is not working");

if (on)
    response = vlt.math.rectify(response);
else
    response = vlt.math.rectify(-response);
end

[max_res, idx_res] = max(response(:));
x_col = X(:); y_col = Y(:); t_col = T(:);
mu = [x_col(idx_res), y_col(idx_res), t_col(idx_res)];
t_range = max(t_col) - min(t_col);
t_resolution = t_range / (size(T, 3) - 1);
% c_unscaled = [size(x, 1) 0 0; 0 size(y, 2) 0; 0 0 t_range];
C = [100*100 0 0; 0 100*100 0; 0 0 0.5*0.5]
C = 0.5*(C+C');
C = C + 3*eye(3);
%C = 0.05*(t_resolution / t_range) * [max(x_col) 0 0; 0 max(y_col) 0; 0 0 t_range];
amp_initial = t_resolution / t_range;
% Upper = [ max(x_col); max(y_col); max(t_col); max(x_col); max(y_col); ...
%     max(t_col); Inf ];
% Lower = [ min(x_col); min(y_col); min(t_col); 0; 0];
StartPoint = [ mu(1); mu(2); mu(3); C(1,1); C(2,2); C(3,3); amp_initial];
mgrid = reshape([x_col, y_col, t_col], 240000, 1, 3);
Upper = [ Inf; Inf; Inf; Inf; Inf; Inf; Inf; ];
Lower = [ -Inf; -Inf; -Inf; -Inf; -Inf; -Inf; 0];

revcorr.mvnpdf_3(StartPoint,mgrid),

x = lsqcurvefit(@(x,xdata) revcorr.mvnpdf_3(x,xdata),StartPoint,mgrid,response(:), Lower, Upper);

mu = [x(1) x(2) x(3)];
C = [x(4) x(5) x(6); x(5) x(7) x(8); x(6) x(8) x(9)];
a = x(10);

% 
% radii = unique(radius);
% response = response(:);  % make sure we have a vector
% x = lsqcurvefit(@(x,xdata) vlt.math.ellipse_on_mvnpdf_x0(x,xdata,X,Y),StartPoint,ellipse_params,response,Lower,Upper);
% 
%  % initial guesses
% amp_initial = [];
% for r=1:length(radii),
% 	inds = find(radius==radii(r));
% 	[amp_initial(r),themaxind] = max(response(inds));
% 	mu = [x_ctr(themaxind) y_ctr(themaxind)];
% 	C = radii(r).^2 * [1 0; 0 1];
% end;
% 
% Upper = [ max(xrange); max(yrange); 10*max(radii)^2; 10*(max(radii)^2); 10*max(radii)^2; Inf ];
% Lower = [ min(xrange); min(yrange); 0; -10*(max(radii)^2); 0;0];
% StartPoint = [ mu(1); mu(2); C(1);C(2);C(4); max(amp_initial)];
% 
% x = lsqcurvefit(@(x) mvnpdf([X(:) Y(:) Z(:)],[x(1),x(2),x(3)],vlt.math.rot3d(x(4),1)*vlt.math.rot3d(x(5),2)*diag([x(6) x(7) x(8)])),StartPoint,ellipse_params,response,Lower,Upper);
% 
% fit_responses = vlt.math.ellipse_on_mvnpdf_x0(x,ellipse_params,X,Y);
% 
% mu = [x(1) x(2)];
% C = [x(3) x(4);x(4) x(5)];
% a = x(6);
% 
% 
