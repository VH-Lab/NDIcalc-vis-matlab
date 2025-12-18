function stim_plot(B,cmap)
%STIM_PLOT - plot the stimulus reconstruction
%
% vis.revcorr.stim_plot(B,CMAP)
%
% Inputs:
%  B - the stimulus reconstruction (MxMxTmax)
%  CMAP - the colormap to use (optional)

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

