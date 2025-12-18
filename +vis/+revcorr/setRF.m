function rf = setRF()
%SETRF Summary of this function goes here
%   Detailed explanation goes here
x = 0:199;
y = 0:199;
[X,~] = meshgrid(x,y);
rf_ = 100 * sin(2 * pi * X * (3 / 200));
rf = zeros(200, 200, 41);
rf(:, :, 20) = rf_;
end

