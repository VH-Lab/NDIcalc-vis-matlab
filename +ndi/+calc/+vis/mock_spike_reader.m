classdef mock_spike_reader < ndi.element
    methods
        function obj = mock_spike_reader(session, varargin)
            if nargin==0, return; end
            if nargin==2 && isa(varargin{1}, 'ndi.document')
                obj = obj@ndi.element(session, varargin{1});
            else
                obj = obj@ndi.element(session, varargin{:});
            end
        end
        function [data, t, timeref] = readtimeseries(obj, timeref_or_epoch, t0, t1)
             % Use the element name to determine the filename
             fname = fullfile(obj.session.path, [obj.name '.mat']);
             if isfile(fname)
                 d = load(fname);
                 st = d.spike_times;
             else
                 st = [];
             end
             if nargin < 3, t0 = -inf; end
             if nargin < 4, t1 = inf; end
             inds = st >= t0 & st <= t1;
             t = st(inds);
             data = ones(size(t));
             timeref = timeref_or_epoch;
        end
    end
end
