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

num_plots = size(B,3);

if num_plots <= 56
    num_rows = 7;
    num_cols = 8;
else
    num_cols = 8;
    num_rows = ceil(num_plots / num_cols);
end

figure();
ax = [];
for i = 1:num_plots
    ax(end+1) = subplot(num_rows,num_cols,i);
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
