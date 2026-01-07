function [b,t] = hartley_stimulus_resampled_time(M, hartley_stimulus_parameters, hartley_stimulus_times, t0, t1, numTimesteps)
% HARTLEY_STIMULUS_RESAMPLED_TIME - generate a space-time picture of the stimulus in a time interval
%
% [B,T] = vis.revcorr.hartley_stimulus_resampled_time(M, HARTLEY_STIMULUS_PARAMETERS, HARTLEY_STIMULUS_TIMES, T0, T1, NUMTIMESTEPS)
%
% Inputs:
%  M - the size of the Hartley stimulus 
%  HARTLEY_STIMULUS_PARAMETERS, a structure with fields
%    s - an array of the values of S for each stimulus
%    kx_v - an array of the values of kx for each stimulus
%    ky_v - an array of the values of ky for each stimulus
%  HARTLEY_STIMULUS_TIMES - an array of the stimulus times in seconds
%  T0 - the time to begin the reconstruction
%  T1 - the time to stop the reconstruction
%  NUMTIMESTEPS - the number of timesteps to be used in the reconstruction
% 
% OUTPUTS:
%  B - an NUMTIMESTEPS array
%    B(:,:,s) is the stimulus that was on the screen at time sample s where
%    the time of stimulus sample s is t0 + (s-1)*deltaT.
%  T = the time of each reconstructed stimulus sample in the same units as
%     T0, T1, HARTLEY_STIMULUS_TIMES

    
    t1 = t1-eps;
    t = linspace(t0,t1,numTimesteps);
    
    S = [];
    
    for i = 1:size(hartley_stimulus_times,1)
        IM = vlt.neuro.reverse_correlation.hartley.hartley_image(hartley_stimulus_parameters.s(i), ...
            hartley_stimulus_parameters.kx_v(i), hartley_stimulus_parameters.ky_v(i), M);
        IM_col = reshape(IM, M*M, 1);
        S = cat(2, S, IM_col);
    end
    
    Si = vlt.math.stepfunc(reshape(hartley_stimulus_times,1,[]),S,t); 
    
    b = [];
    for i = 1:size(t,2)
        IM = reshape(Si(:,i), M, M);
        b = cat(3, b, IM);
    end
end