function [t0,t1] = get_t0(t1,delta_t)
% GET_T0 - returns the time to begin and end the reconstruction
%
% [T0,T1] = GET_T0(T1,DELTA_T)
%
% Inputs:
%  T1 - the time to end the reconstruction
%  DELTA_T - the time duration of the reconstruction
% 
% OUTPUTS:
%  T0 - the time to begin the reconstruction
%  T1 - the time to end the reconstruction

t0 = t1 - delta_t;
end