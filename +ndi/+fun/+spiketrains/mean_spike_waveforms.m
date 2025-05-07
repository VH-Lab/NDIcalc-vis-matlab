function [mean_waves,std_waves,parameters] = mean_spike_waveforms(spiking_element, epoch_id, varargin)
% MEAN_SPIKE_WAVEFORMS - find mean spike waveforms in a voltage recording
%
% [MEAN_WAVES,STD_WAVES,PARAMETERS] = MEAN_SPIKE_WAVEFORMS(S, SPIKING_ELEMENT, EPOCH_ID, ...)
%
% Given an ndi.element of type 'spikes' that is derived from a voltage (SPIKING_ELEMENT),
% recording, returns mean and standard deviation waveforms from different places
% in the recording. The mean is computed in different intervals in the recording EPOCH_ID 
% specified by a parameter 'averaging_window', which defaults to 60 s. This means that
% mean waveforms are computed in the first averaging_window, the second averaging_window,
% etc., until the end of the recording. 
%
% WAVES is an MxCxT matrix where M is the number of samples around the spike
%   that are examined, C is the number of channels, and T is the number of
%   times that the mean waveform is examined. Any bins with no spikes will be
%   coded as NaN values.
% STD_WAVES is an MxCxT matrix, the same size as WAVES, that has the standard
%   deviation of the spike waves. Any bins with no spikes will be coded as 
%   NaN values.
% PARAMETERS is a structure with the following fields:
%   interval_center_times - the center time of each averaging_window
%   number_of_spikes_per_interval - the number of spikes in each interval
%   sample_times - the relative times of each spike waveform
%   s0 - number of samples before spike sample that were read
%   s1 - number of samples after spike sample that were read
%   sample_rate - the sampling rate of the probe
%   
% 
% This function takes name/value parameters that modify its behavior:
% ---------------------------------------------------------------------------
% | Parameter (default)      | Description                                  |
% |--------------------------|----------------------------------------------|
% | spike_window_before_time | How long before each spike time should we    |
% |  (-0.001)                |   retrieve for our waveforms? (In seconds.)  |
% |                          |   Usually negative to indicate times before  |
% |                          |   the spike.                                 |
% | spike_window_after_time  | How long after each spike time should we     |
% |  ( 0.002)                |   retrieve for our waveforms? (In seconds.)  |
% |                          |   Usually positive to indicate times after   |
% |                          |   the spike.                                 |
% | averaging_window (60)    | How long should the averaging window be?     |
% |                          |   (Seconds.) A mean and standard deviation of|
% |                          |   the waveform is constructed every interval.|
% | averaging_window_step    | How long should the step be between averaging|
% |       (300)              |   windows? (Seconds.)
% | filter_padding (0.100)   | Extra time to read before and after each     |
% |                          |   spike for filtering.                       |
% | cheby_order (4)          | Chebyshev Type I filter order                |
% | cheby_R (0.5)            | Chebyshev Type I filter roll off R parameter |
% | cheby_cutoff (300)       | Chebyshev Type I high pass cut off frequency |
% |--------------------------|----------------------------------------------|
%

spike_window_before_time = -0.001;
spike_window_after_time = 0.002;
averaging_window= 60;
averaging_window_step = 300;
filter_padding = 0.050;

cheby_order = 4;
cheby_R = 0.5;
cheby_cutoff = 300;

vlt.data.assign(varargin{:});

ndi.globals();

et = spiking_element.epochtableentry(epoch_id);

underlying_probe = et.underlying_epochs(1).underlying;
underlying_epoch = et.underlying_epochs(1).epoch_id;

match = 0;

for i=1:numel(et.underlying_epochs.epoch_clock),
	if strcmp(et.underlying_epochs.epoch_clock{i}.type,'dev_local_time'),
		match = i;
		break;
	end;
end;

sample_rate = underlying_probe.samplerate(underlying_epoch);
t0t1 = et.underlying_epochs.t0_t1{match};
t0 = t0t1(1);
t1 = t0t1(2);

data = underlying_probe.readtimeseries(underlying_epoch,t0,t0+1/sample_rate);
num_channels = size(data,2);

[B,A] = cheby1(cheby_order,cheby_R,cheby_cutoff/(0.5*sample_rate),'high');

[dummy,spiketimes] = spiking_element.readtimeseries(epoch_id,-inf,inf);

if t1<averaging_window,
	interval_edges = [t0 t1];
else,
	interval_edges = t0:averaging_window_step:t1;
	if interval_edges(end)~=t1,
		interval_edges(end+1) = t1;
	end;
end;

start_times = [];
end_times = [];
interval_center_times = [];
number_of_spikes_per_interval = [];
s0 = round(spike_window_before_time*sample_rate);
s1 = round(spike_window_after_time*sample_rate);
s000 = round((spike_window_before_time-filter_padding)*sample_rate);
s001 = round((spike_window_after_time+filter_padding)*sample_rate);

relative_sample_read = s000:s001;
relative_time_reads = relative_sample_read/sample_rate;

first_sample = s0 - s000 + 1;

sample_times = [s0:s1]/sample_rate;

mean_waves = [];
std_waves = [];

for t=1:numel(interval_edges)-1,
	if t~=numel(interval_edges)-1, % not last one
		start_times(t) = interval_edges(t);
		end_times(t) = start_times(t)+averaging_window;
	else, % last one
		end_times(t) = interval_edges(end);
		start_times(t) = end_times(t) - averaging_window;
	end;
	spikes_here = find(spiketimes>=start_times(t) & spiketimes<=end_times(t));
	number_of_spikes_per_interval(t) = numel(spikes_here);
	ndi_globals.log.msg('debug',5,['mean_spike_waveform: Analyzing interval ' int2str(t) ' of ' int2str(numel(interval_edges)-1) ', ' int2str(numel(spikes_here)) ' spike waveforms.']);
	interval_center_times(t) = mean([start_times(t) end_times(t)]);
	waves_here = NaN(s1-s0+1,num_channels,0);
	for s = 1:numel(spikes_here),
		[data,t_here] = underlying_probe.readtimeseries(underlying_epoch,...
			spiketimes(spikes_here(s))+relative_time_reads(1),...
			spiketimes(spikes_here(s))+relative_time_reads(end));

		if numel(t_here)==numel(relative_sample_read), % we got a complete read
			data_here = filtfilt(B,A,data);
			% now have to find samples we actually want
			waves_here = cat(3,waves_here,data_here(first_sample:first_sample+(s1-s0),:));
		else, % we have an edge case, we can't filter it, don't count it
			number_of_spikes_per_interval(t) = number_of_spikes_per_interval(t) - 1;
		end;
	end;
	mean_waves = cat(3, mean_waves,nanmean(waves_here,3));
	std_waves = cat(3, std_waves, nanstd(waves_here,[],3));
end;

parameters = struct('interval_center_times',interval_center_times,...
	'number_of_spikes_per_interval',number_of_spikes_per_interval,...
	'sample_times',sample_times,'s0',s0,'s1',s1,'sample_rate',sample_rate);

