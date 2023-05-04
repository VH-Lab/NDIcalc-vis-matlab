function hartley_plot_demo(S, epoch_id, elements, varargin)
% HARTLEY_PLOT_DEMO - plot receptive fields for several Hartley reconstructions 
%
% HARTLEY_PLOT_DEMO(S, EPOCH_ID, ELEMENTS, ...)
%
% 
% Inputs: 
%   S - an ndi.session object
%   EPOCHID - an epochid to plot:w
%   ELEMENTS - a cell array of ndi.element objects (type: spikes) to plot
%
% This function also takes name/value pairs that modify its default behavior:
% ---------------------------------------------------------------------------
% | Parameter (default)       | Description                                 |
% |---------------------------|---------------------------------------------|
% | latency_xy (0.050)        | Latency (in seconds) to show in the XY      |
% |                           |   panel.                                    |
% | Y_axis (1/2 of height)    | Y axis plotting location                    |
% | X_axis (1/2 of width)     | X axis plotting location                    |
% | axis_color ([1 1 0])      | Color to use to plot axis                   |
% ---------------------------------------------------------------------------
%  
%  

latency_xy = 0.050;
Y_axis = [];
X_axis = [];
axis_color = [1 1 0];

vlt.data.assign(varargin{:});

if ~isa(S,'ndi.session'),
	error(['S must be an ndi.session object.']);
end;

hc = ndi.calc.vis.hartley_calc(S);
shc = ndi.calc.vis.spike_shape(S);

q1 = ndi.query('','isa','hartley_calc');
q2 = ndi.query('','isa','spike_shape_calc');

p_q = ndi.query('','isa','stimulus_presentation') & ndi.query('epochid.epochid','exact_string',epoch_id);

p_obj = S.database_search(p_q);
if numel(p_obj)~=1,
	error(['Did not find exactly one stimulus presentation document.']);
end;
p_obj = p_obj{1};

xy_axes = {};
rt_axes = {};
sp_axes = {};

for i=1:numel(elements),
	h_q1 = ndi.query('','depends_on','element_id',elements{i}.id());
	h_q2 = ndi.query('','depends_on','stimulus_presentation_id',p_obj.id());
	element_epoch = S.database_search(ndi.query('','isa','element_epoch') ...
		& h_q1 & ndi.query('epochid.epochid','exact_string',epoch_id));
	if numel(element_epoch)~=1, error('wrong number of element_epochs.'); end;
	s_q2 = ndi.query('','depends_on','element_epoch_id',element_epoch{1}.id());

	h = S.database_search(q1 & h_q1 & h_q2);
	sh = S.database_search(q2 & h_q1 & s_q2); 

	if numel(h)~=1, error(['Found wrong number of hartley_calc docs: ' int2str(numel(h)) '.']); end;
	if numel(sh)~=1, error(['Found wrong number of spike_shape_calc docs: ' int2str(numel(sh)) '.']); end;

	h = h{1};
	sh = sh{1};

	frame_index = vlt.data.findclosest(h.document_properties.hartley_reverse_correlation.reconstruction_properties.T_coords,latency_xy);

	% read Hartley reconstruction
	[sta,pval] = hc.read_sta(h);
	significance_plot = revcorr.rescale_p_image(pval);
	cmap = revcorr.get_cmap();

	rp = h.document_properties.hartley_reverse_correlation.reconstruction_properties;
	[theta,sta_r,pval_r] = revcorr.rotate_sta(rp.T_coords, sta, significance_plot);
	[t_profile,t_profile_pval] = revcorr.peak_time_profile(rp.T_coords,sta,significance_plot);

	if isempty(Y_axis),
		Y_axis_here = 0.5 * size(sta,2); % height is second axis in images
	else,
		Y_axis_here = Y_axis;
	end;

	if isempty(X_axis),
		X_axis_here = 0.5 * size(sta,1); % height is second axis in images
	else,
		X_axis_here = XY_axis;
	end;

	xy_axes{i} = subplot(numel(elements), 3, 1+(i-1)*3);
	
	image(significance_plot(:,:,frame_index));
	colormap(xy_axes{i},cmap);
	axis equal off;

	hold on;
	plot(X_axis_here*[1 1],[0 size(sta,2)],'-','color',axis_color);
	plot([0 size(sta,1)],Y_axis_here*[1 1],'-','color',axis_color);

	rt_axes{i} = subplot(numel(elements), 3, 2+(i-1)*3);
	image(t_profile_pval,'XData',rp.T_coords,'YData',rp.Y_coords);
	colormap(rt_axes{i},cmap);
	axis normal on
	grid on;

	sp_axes{i} = subplot(numel(elements), 3, 3+(i-1)*3);
	[mean_waves,std_waves,sample_times] = shc.load(sh);

	if 0,
		ndi.fun.plot.multichan(mean_waves(:,:,1),sample_times,30);
	else,
		plot(sample_times,mean_waves(:,:,1));
	end;
 
end;

linkaxes(xy_axes{:});

linkaxes(rt_axes{:});


