function r = mvnpdf_3(x, xdata)
% MVNPDF_3 - fit the data to a mvnpdf function
%
% R = vis.revcorr.mvnpdf_3(X, XDATA)
%
% Inputs:
%  X - parameters for the gaussian [mu_x mu_y mu_t sigma_1 ... sigma_3 scale]
%  XDATA - the data points [X Y T]
%
% Output:
%  R - the evaluated function values
%
% ERROR: this function is not working

error("ERROR: this function is not working");

X = xdata(:,1);
Y = xdata(:,2);
T = xdata(:,3);

mu = [x(1) x(2) x(3)];
%C = vlt.math.rot3d(x(7),1)*vlt.math.rot3d(x(8),2)*[ x(4) 0 0 ; 0 x(5) 0 ; 0 0 x(6)];
C = [x(4) ; x(5); x(6)] * [x(4); x(5); x(6)]';
%x(4)=max(1e-3,x(4));
%x(7)=max(1e-3,x(7));
%x(9)=max(1e-3,x(9));
%C = [x(4) x(5) x(6); x(5) x(7) x(8); x(6) x(8) x(9)];
a = x(7);
C
r = a * mvnpdf([X Y T], mu, C);
end

