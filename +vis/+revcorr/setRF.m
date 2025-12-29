function [rf, timelags] = setRF(M, num_timesteps, dt)
%SETRF - set the receptive field for the simulation
%
% [RF, TIMELAGS] = vis.revcorr.setRF(M, NUM_TIMESTEPS, DT)
%
% OUTPUTS:
%  RF - the receptive field (MxMxNUM_TIMESTEPS), where the first two dimensions
%       are space and the 3rd dimension is time.
%  TIMELAGS - a vector of kernel times (going from 0 to positive), corresponding
%             to the 3rd dimension of RF.
%
% INPUTS:
%  M - the spatial dimension of the receptive field (result is MxM)
%  NUM_TIMESTEPS - the number of time steps in the receptive field
%  DT - the duration of each time step

x = 0:M-1;
y = 0:M-1;
[X,~] = meshgrid(x,y);
rf_ = 100 * sin(2 * pi * X * (3 / M));
rf = zeros(M, M, num_timesteps);

% Place the RF impulse roughly in the middle or at a fixed index
% The original code used index 20 out of 42.
% We will use round(num_timesteps/2) as a default position for the impulse.

impulse_idx = max(1, round(num_timesteps * (1/4)));

rf(:, :, impulse_idx) = rf_;

timelags = 0:dt:(num_timesteps-1)*dt;

end
