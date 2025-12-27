function stim_plot(B,cmap,timelags)
%STIM_PLOT - plot the stimulus reconstruction
%
% vis.revcorr.stim_plot(B,CMAP,TIMELAGS)
%
% Inputs:
%  B - the stimulus reconstruction (MxMxTmax)
%  CMAP - the colormap to use (optional, default [])
%  TIMELAGS - the time lag of each frame (optional, default [])

if nargin < 2
    cmap = [];
end
if nargin < 3
    timelags = [];
end

figure();
ax = [];
for i = 1:size(B,3)
    ax(end+1) = subplot(7,8,i);
    if isempty(cmap)
        imshow(B(:,:,i), []);
    else
        imshow(B(:,:,i), cmap);
    end
    if ~isempty(timelags)
        if i <= length(timelags)
            xlabel(sprintf('%.3f', timelags(i)));
        end
    end
end
linkaxes(ax);
end
