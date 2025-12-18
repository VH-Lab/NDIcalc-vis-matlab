function [hartley_stimulus_parameters, hartley_stimulus_times] = get_frames(s,kx_v, ky_v, frameTimes, t0, t1)
%GET_FRAMES Summary of this function goes here
%   Detailed explanation goes here
    
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

