function y = dog(x,p_dog)
% DOG - difference of gaussians
%
% Y = DOG(X, P_DOG)
%
% Given values of X and difference-of-gaussian parameters
% P_DOG = [ A1 B1 A2 B2]
%
% returns 
%
% y = A1*exp(-x.^2/(2*B1^2)) - A2*exp(-x.^2/(2*B2^2))
%
% (By Alexander Heimel)
%

a1 = p_dog(1);
b1 = p_dog(2);
a2 = p_dog(3);
b2 = p_dog(4);

y = a1*exp(-x.^2/(2*b1^2)) - a2*exp(-x.^2/(2*b2^2));
