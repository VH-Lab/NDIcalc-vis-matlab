function [theta,sta,prob_sta] = rotate_sta(T_coords,sta,prob_sta)
% ROTATE_STA - rotate the spike-triggered average to the prevailing orientation
%
% [THETA,STA,PROB_STA] = ROTATE_STA(T_COORDS,STA,PROB_STA)
%
% Find the prevailing orientation THETA (in radians) for a spike-triggered
% average. T_COORDS are the time coordinates (in seconds) of the spike-triggered% average and STA is an XxYxT matrix with the STA. PROB_STA is an XxYxT matrix
% with the significance of each value of the STA.
%
% The prevailing orientation is computed from image frames with time lags
% between 0 and 0.1 seconds.
% 

low_time = 0;
high_time = 0.1;


t_good = T_coords>=low_time & T_coords <= high_time;

theta = vlt.image.prevailing_orientation(mean(sta(:,:,t_good),3));
theta_degrees = vlt.math.rad2deg(theta);

for i=1:size(sta,3),
	sta(:,:,i) = imrotate(sta(:,:,i),theta_degrees,'bilinear','crop');
end;
for i=1:size(prob_sta,3),
	prob_sta(:,:,i) = imrotate(prob_sta(:,:,i),theta_degrees,...
		'bilinear','crop');
end;



