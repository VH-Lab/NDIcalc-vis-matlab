function [t_profile,t_profile_pval] = peak_time_profile(T_coords,sta,pval_sta)
% PEAK_TIME_PROFILE - compute a space x time profile for a receptive field
%
% [T_PROFILE,T_PROFILE_PVAL] = vis.revcorr.peak_time_profile(T_COORDS,STA,PVAL_STA)
%
% Computes a space x time slice through an STA. The first dimension is
% exampled for a peak location. The signal is first blurred by sliding a
% 5-pixel boxcar filter over the data.  Then, the 5 slices at the 
% spatial peak location are averaged together to create a space x time
% profile.
%
% The same pixels are selected for PVAL_STA, which is assumed to be an 
% alternative projection of the STA onto probability space, where
% PVAL_STA(x,y,t) indicates the likelihood that the pixel is significant.
%

low_time = 0;
high_time = 0.1;
blur_window = 5;
avg_window = 4;

t_good = find(T_coords>=low_time & T_coords <= high_time);

sta_mean = mean(sta(:,:,t_good),3);

myvector = sum(abs(sta_mean),2);
myvector = conv(myvector(:),ones(blur_window,1)/blur_window,'same');

[myvalue,location] = max(myvector);

span = max(1,location-avg_window) : min(numel(myvector),location+avg_window) ;

t_profile = squeeze(mean(sta(span,:,:),1));
t_profile_pval = squeeze(mean(pval_sta(span,:,:),1));


