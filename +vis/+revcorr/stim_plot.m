function stim_plot(B,cmap)
%STIM_PLOT Summary of this function goes here
%   Detailed explanation goes here
figure();
ax = [];
for i = 1:size(B,3)
    if nargin == 1
        ax(end+1) = subplot(7,8,i);
	imshow(B(:,:,i), []);
    else
        ax(end+1) = subplot(7,8,i);
	imshow(B(:,:,i), cmap);
    end
end
linkaxes(ax);
end

