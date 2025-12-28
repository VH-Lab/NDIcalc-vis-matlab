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

   if t0 < frameTimes(1) && t1 < frameTimes(1)
       hartley_stimulus_parameters.s = [0; 0];
       hartley_stimulus_parameters.kx_v = [0; 0];
       hartley_stimulus_parameters.ky_v = [0; 0];
       hartley_stimulus_times = [t0; t1];
       return;
   end

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
   if isempty(hartley_stimulus_times) || t0 ~= hartley_stimulus_times(1)
       i = find(idx, 1, 'first');
       if isempty(i)
           % if the interval is empty, we need to find the frame that is currently displayed at t0
           % this corresponds to the last frame before t0
           i = find(frameTimes < t0, 1, 'last');
           if isempty(i)
               % if there is no frame before t0, we use dummy values
               hartley_stimulus_parameters.s = vertcat(0,hartley_stimulus_parameters.s);
               hartley_stimulus_parameters.kx_v = vertcat(0,hartley_stimulus_parameters.kx_v);
               hartley_stimulus_parameters.ky_v = vertcat(0,hartley_stimulus_parameters.ky_v);
               hartley_stimulus_times = vertcat(t0,hartley_stimulus_times);
               return;
           else
               % if there is a frame before t0
               % however, we need to be careful if t0 > t1. In this case,
               % the loop in test.m provides (responseTimes(i), responseTimes(i)-rfTimeRange).
               % Ideally get_frames should be called with (t_start, t_end).
               % If t0 > t1, idx is empty.
               % Let's assume standard usage t0 <= t1.
               % If t0 > t1 is passed, the caller might be expecting backward time?
               % In test.m: get_frames(..., responseTimes(i), responseTimes(i)-rfTimeRange)
               % This implies t0 > t1.
               % If the intention is to get frames in the past of responseTime(i),
               % the call should probably be (responseTimes(i)-rfTimeRange, responseTimes(i)).

               % But assuming we must fix get_frames to not crash:
               hartley_stimulus_parameters.s = vertcat(s(i),hartley_stimulus_parameters.s);
               hartley_stimulus_parameters.kx_v = vertcat(kx_v(i),hartley_stimulus_parameters.kx_v);
               hartley_stimulus_parameters.ky_v = vertcat(ky_v(i),hartley_stimulus_parameters.ky_v);
               hartley_stimulus_times = vertcat(t0,hartley_stimulus_times);
               return;
           end
       end
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
