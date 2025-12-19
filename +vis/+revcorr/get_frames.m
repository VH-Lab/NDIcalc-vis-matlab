function [hartley_stimulus_parameters, hartley_stimulus_times] = get_frames(s,kx_v, ky_v, frameTimes, t0, t1)
%GET_FRAMES - extract the hartley stimulus parameters and time points
% in a given interval [t0, t1]
%
% [HARTLEY_STIMULUS_PARAMETERS, HARTLEY_STIMULUS_TIMES] = vis.revcorr.get_frames(S, KX_V, KY_V, FRAMETIMES, T0, T1)
%
% Inputs:
%  S - the s parameters of the hartley stimulus
%  KX_V - the kx parameters of the hartley stimulus
%  KY_V - the ky parameters of the hartley stimulus
%  FRAMETIMES - the time points of the hartley stimulus
%  T0 - the start time of the interval
%  T1 - the end time of the interval
%
% OUTPUTS:
%  HARTLEY_STIMULUS_PARAMETERS - a structure containing the s, kx, and ky parameters
%  HARTLEY_STIMULUS_TIMES - the time points of the hartley stimulus in the interval
    
   idx = frameTimes >= t0 & frameTimes <= t1;
   hartley_stimulus_parameters = struct("s", s(idx), "kx_v", kx_v(idx), "ky_v", ky_v(idx));
   hartley_stimulus_times = frameTimes(idx);

   %% check the corner cases
   if t0 < frameTimes(1)
       hartley_stimulus_parameters.s = vertcat(0,hartley_stimulus_parameters.s);
       hartley_stimulus_parameters.kx_v = vertcat(0,hartley_stimulus_parameters.kx_v);
       hartley_stimulus_parameters.ky_v = vertcat(0,hartley_stimulus_parameters.ky_v);
       hartley_stimulus_times = vertcat(t0,hartley_stimulus_times);
   end

   if t1 > frameTimes(end)
       hartley_stimulus_parameters.s = vertcat(hartley_stimulus_parameters.s, 0);
       hartley_stimulus_parameters.kx_v = vertcat(hartley_stimulus_parameters.kx_v, 0);
       hartley_stimulus_parameters.ky_v = vertcat(hartley_stimulus_parameters.ky_v, 0);
       hartley_stimulus_times = vertcat(hartley_stimulus_times, t1);
   end

   if t0 > frameTimes(end)
       hartley_stimulus_parameters.s = vertcat(0,hartley_stimulus_parameters.s);
       hartley_stimulus_parameters.kx_v = vertcat(0,hartley_stimulus_parameters.kx_v);
       hartley_stimulus_parameters.ky_v = vertcat(0,hartley_stimulus_parameters.ky_v);
       hartley_stimulus_times = vertcat(t0,hartley_stimulus_times);
   end

   %% insert the frames that is currently displayed at time t0
   if t0 ~= hartley_stimulus_times(1)
       i = find(idx, 1, 'first');
       hartley_stimulus_parameters.s = vertcat(s(i - 1),hartley_stimulus_parameters.s);
       hartley_stimulus_parameters.kx_v = vertcat(kx_v(i - 1),hartley_stimulus_parameters.kx_v);
       hartley_stimulus_parameters.ky_v = vertcat(ky_v(i - 1),hartley_stimulus_parameters.ky_v);
       hartley_stimulus_times = vertcat(t0,hartley_stimulus_times);
   end

   %% insert the frame after the interval
   if t1 ~= hartley_stimulus_times(end)
       i = find(idx, 1, 'last');
       if isempty(i)
           i = 0;
       end
       i = min(i, size(frameTimes, 1) - 1);
       hartley_stimulus_parameters.s = vertcat(hartley_stimulus_parameters.s, s(i + 1));
       hartley_stimulus_parameters.kx_v = vertcat(hartley_stimulus_parameters.kx_v, kx_v(i + 1));
       hartley_stimulus_parameters.ky_v = vertcat(hartley_stimulus_parameters.ky_v, ky_v(i + 1));
       hartley_stimulus_times = vertcat(hartley_stimulus_times,frameTimes(i + 1));
   end
end

