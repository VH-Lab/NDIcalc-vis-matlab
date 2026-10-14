function fh = plotResults(D, options)
% PLOTRESULTS - Plot the stored calculator results in the paper's example dataset
%
% FH = NDI.CALC.PAPERFIGURES.PLOTRESULTS()
% FH = NDI.CALC.PAPERFIGURES.PLOTRESULTS(D)
%
% Plots the calculator output documents that are stored in the example
% dataset that accompanies the calculators paper. For each session in the
% dataset and each calculator, the documents of that calculator's output type
% are found with a database search, and each one is drawn in its own figure
% with the calculator's own PLOT method.
%
% D is the ndi.dataset to plot from. If it is not given, the dataset is opened
% (and downloaded if necessary) with ndi.calc.paperFigures.openDataset.
%
% FH is a vector of the figure handles that were created.
%
% This function also takes name/value pairs that modify its behavior:
% |----------------------------|----------------------------------------------|
% | Parameter (default)        | Description                                  |
% |----------------------------|----------------------------------------------|
% | calculators                | The calculators whose results to plot, by    |
% |   (all seven)              | class name, e.g. "ndi.calc.vis.hartley".     |
% | sessions ([])              | Which of the dataset's sessions to plot, by  |
% |                            | index into D.session_list(). Empty means all.|
% |----------------------------|----------------------------------------------|
%
% Example:
%   D = ndi.calc.paperFigures.openDataset();
%   ndi.calc.paperFigures.plotResults(D, 'calculators', "ndi.calc.vis.hartley");
%
% See also: ndi.calc.paperFigures.openDataset, ndi.calculator/plot

arguments
    D = []
    options.calculators (1,:) string = [ ...
        "ndi.calc.vis.oridir_tuning", ...
        "ndi.calc.vis.contrast_tuning", ...
        "ndi.calc.vis.contrast_sensitivity", ...
        "ndi.calc.vis.spatial_frequency_tuning", ...
        "ndi.calc.vis.temporal_frequency_tuning", ...
        "ndi.calc.vis.speed_tuning", ...
        "ndi.calc.vis.hartley"]
    options.sessions (1,:) double = []
end

if isempty(D)
    D = ndi.calc.paperFigures.openDataset();
end

% The output document type of each calculator
docTypes = dictionary( ...
    "ndi.calc.vis.oridir_tuning",             "oridirtuning_calc", ...
    "ndi.calc.vis.contrast_tuning",           "contrasttuning_calc", ...
    "ndi.calc.vis.contrast_sensitivity",      "contrastsensitivity_calc", ...
    "ndi.calc.vis.spatial_frequency_tuning",  "spatial_frequency_tuning_calc", ...
    "ndi.calc.vis.temporal_frequency_tuning", "temporal_frequency_tuning_calc", ...
    "ndi.calc.vis.speed_tuning",              "speedtuning_calc", ...
    "ndi.calc.vis.hartley",                   "hartley_calc");

unknown = options.calculators(~isKey(docTypes, options.calculators));
if ~isempty(unknown)
    error('NDICALC:paperFigures:plotResults:UnknownCalculator', ...
        'Unknown calculator(s): %s.', strjoin(unknown, ', '));
end

[~, sessionIds] = D.session_list();
if isempty(options.sessions)
    options.sessions = 1:numel(sessionIds);
end

fh = gobjects(0);

for s = options.sessions
    S = D.open_session(sessionIds{s});
    for calcName = options.calculators
        c = feval(char(calcName), S);
        q = ndi.query('', 'isa', char(docTypes(calcName)));
        docs = S.database_search(q);
        for i = 1:numel(docs)
            fh(end+1) = figure; %#ok<AGROW>
            c.plot(docs{i});
        end
    end
end

end
