function plottuning(SF,TF,R,varargin)
% PLOTTUNING - Plot speed tuning
%
%  PLOTTUNING(SF,TF,R,...)
%
%  Creates a plot like Priebe et al. 2006.
%
%  On the left side, plot the speed tuning for each spatial frequency.
%  On the right side, plot a "heat map" of response as a function of
%  spatial frequency.
%
%  This function also takes name/value pairs that modify its
%  default behavior:
%  |--------------------------------------------------------------------|
%  |Parameter (default)         | Description                           |
%  |----------------------------|---------------------------------------|
%  |'marker' ('o')              | Marker type to use in plot            |
%  |'linestyle' ('none')        | Line style to use                     |
%  |'do_surf' (1)               | 0/1 Should we do the surface plot?    |
%  |'surf_colormap' (gray(256)) | Colormap for the heat map. The default|
%  |                            |   ramps from black (low) to white     |
%  |                            |   (high).                             |
%  |'half_contour_params' ([])  | If non-empty, a 7-element speed tuning |
%  |                            |   parameter vector (see               |
%  |                            |   vis.speed.tuningfunc). The half-peak |
%  |                            |   (half-height) contour of that model  |
%  |                            |   is overlaid on the heat map.         |
%  |'contour_color' ([1 0 0])   | Color of the overlaid half-peak contour|
%  |'contour_linewidth' (1.5)   | Line width of the overlaid contour     |
%  |----------------------------|---------------------------------------|
%
%  Example:
%    out = vis.speed.test.responseplay();
%    figure;
%    vis.speed.plottuning(out.SF,out.TF,out.R,'marker','d','linestyle','-');
%

marker = 'o';
linestyle = 'none';
do_surf = 1;
surf_colormap = gray(256);
half_contour_params = [];
contour_color = [1 0 0];
contour_linewidth = 1.5;
vlt.data.assign(varargin{:});

% On the left side
subplot(1,2,1);
hold on;

all_sfs = unique(SF); % All the spatial frequencies present

colors = jet(numel(all_sfs));

for s = 1:numel(all_sfs)
    indexes = find(SF==all_sfs(s));
    % Plot the responses and calculate the speed for each spot
    speed_here = TF(indexes)./SF(indexes);
    h = plot(speed_here,R(indexes),'marker',marker,'linestyle',linestyle,'color',colors(s,:));
end
set(gca,'XScale','log');
set(gca,'FontAngle','italic');

xlabel('Speed (deg/s)','FontAngle','italic');
ylabel('Response','FontAngle','italic');

% On the right side
if ~do_surf, 
    return; % stop if we aren't plotting the surface
end;

subplot(1,2,2);

% Use surf to draw the heat map. The data live in linear SF/TF coordinates
% but are displayed on log-scaled axes. To keep the flat-shaded coloring
% aligned with the true coordinates (so that an overlaid contour lines up
% with the colors), we place the cell *edges* at the geometric midpoints
% between adjacent SF/TF values in log space. This centers each data point
% within its colored cell rather than leaving it at a corner.
sf_edges = local_log_edges(SF(1,:)); % 1 x (n+1) vertex coordinates
tf_edges = local_log_edges(TF(:,1)); % 1 x (m+1) vertex coordinates
[SF_surf, TF_surf] = meshgrid(sf_edges, tf_edges); % (m+1) x (n+1) vertices

% Pad the color/height matrix to the vertex grid size. With flat shading,
% face (i,j) takes its color from C(i,j) = R(i,j); the extra row/column are
% never used for face color and are present only to match the grid size.
R_surf = [ R R(:,end) ; R(end,:) R(end,end)];

surf(SF_surf,TF_surf,R_surf,R_surf);
set(gca,'View',[0 90]);
set(gca,'XScale','log','YScale','log');
set(gca,'XTick',[SF(1,:)]);
set(gca,'YTick',[TF(:,1)]);
set(gca,'FontAngle','italic');
axis([sf_edges(1) sf_edges(end) tf_edges(1) tf_edges(end)]);
shading faceted;
colormap(gca,surf_colormap);

% Optionally overlay the half-peak (half-height) contour of a fitted model.
if ~isempty(half_contour_params)
    [cX, cY, peak_sf, peak_tf] = vis.speed.halfpeakcontour(half_contour_params);
    cX(end+1) = cX(1); % close the contour loop
    cY(end+1) = cY(1);
    z_top = max(R_surf(:)) + 1; % draw above the surface so it is visible
    hold on;
    plot3(cX, cY, z_top*ones(size(cX)), '-', ...
        'Color', contour_color, 'LineWidth', contour_linewidth);
    plot3(peak_sf, peak_tf, z_top, '+', ...
        'Color', contour_color, 'LineWidth', contour_linewidth, 'MarkerSize', 8);
end

xlabel('Spatial frequency (c/deg)','FontAngle','italic');
ylabel('Temporal frequency (Hz)','FontAngle','italic');

end % plottuning

% -------------------------------------------------------------------------

function edges = local_log_edges(v)
% LOCAL_LOG_EDGES - compute cell-edge coordinates centered on the values V.
%
% Given a monotonic vector V of (positive) coordinate values, returns a
% vector of length numel(V)+1 whose elements bracket each value at the
% geometric midpoint (midpoint in log10 space) between neighbors, so that
% each value sits at the center of its cell.
    v = v(:).';
    lv = log10(v);
    if numel(lv) == 1
        e = [lv - 0.5, lv + 0.5];
    else
        mids = (lv(1:end-1) + lv(2:end)) / 2;
        e = [lv(1) - (lv(2)-lv(1))/2, mids, lv(end) + (lv(end)-lv(end-1))/2];
    end
    edges = 10.^e;
end % local_log_edges
