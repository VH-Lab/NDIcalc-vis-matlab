function fh = plotResults(D, options)
% PLOTRESULTS - Plot the calculator results shown in the paper from the example dataset
%
% FH = NDI.CALC.PAPERFIGURES.PLOTRESULTS(D, ...)
%
% Plots the calculator output documents that the calculators paper shows,
% drawn from the example dataset that accompanies the paper. Each calculator
% has its own section, switched by its own name/value pair; all are on by
% default, since plotting is quick, so name the ones to leave out. Each section finds the output
% documents of the example cells (ndi.calc.paperFigures.exampleCells) and
% draws them with the calculator's own PLOT method.
%
% The orientation/direction, contrast, contrast sensitivity, spatial
% frequency and temporal frequency examples come from the dataset's session
% from Griswold and Van Hooser (2025) and are drawn in windows of 4 x 4 plots.
% The speed tuning and Hartley examples come from the session from the Suarez
% Casanova et al. study; their plots have several panels, so each gets its own
% window. A cell has one Hartley document per epoch of Hartley stimuli; only
% the last is shown (ndi.calc.paperFigures.lastHartley).
%
% D is the ndi.dataset to plot from. If it is empty, the dataset is opened
% (and downloaded if necessary) with ndi.calc.paperFigures.openDataset.
%
% FH is a vector of the figure handles that were created.
%
% This function also takes name/value pairs that modify its behavior:
% |----------------------------|----------------------------------------------|
% | Parameter (default)        | Description                                  |
% |----------------------------|----------------------------------------------|
% | doOridir (true)            | Plot the ndi.calc.vis.oridir_tuning examples.|
% | doContrast (true)          | Plot the ndi.calc.vis.contrast_tuning        |
% |                            |   examples.                                  |
% | doContrastSensitivity      | Plot the ndi.calc.vis.contrast_sensitivity   |
% |   (true)                   |   examples.                                  |
% | doSpatialFrequency (true)  | Plot the                                     |
% |                            |   ndi.calc.vis.spatial_frequency_tuning      |
% |                            |   examples.                                  |
% | doTemporalFrequency (true) | Plot the                                     |
% |                            |   ndi.calc.vis.temporal_frequency_tuning     |
% |                            |   examples.                                  |
% | doSpeed (true)             | Plot the ndi.calc.vis.speed_tuning examples. |
% | doHartley (true)           | Plot the ndi.calc.vis.hartley examples.      |
% |----------------------------|----------------------------------------------|
%
% Example:
%   D = ndi.calc.paperFigures.openDataset();
%   fh = ndi.calc.paperFigures.plotResults(D);                       % everything
%   fh = ndi.calc.paperFigures.plotResults(D, 'doContrast', false);  % all but contrast
%
% See also: ndi.calc.paperFigures.openDataset, ndi.calc.paperFigures.exampleCells,
%   ndi.calc.paperFigures.recomputeResults, ndi.calculator/plot

arguments
    D = []
    options.doOridir (1,1) logical = true
    options.doContrast (1,1) logical = true
    options.doContrastSensitivity (1,1) logical = true
    options.doSpatialFrequency (1,1) logical = true
    options.doTemporalFrequency (1,1) logical = true
    options.doSpeed (1,1) logical = true
    options.doHartley (1,1) logical = true
end

if isempty(D)
    D = ndi.calc.paperFigures.openDataset();
end

fh = gobjects(1,0);

% The Griswold and Van Hooser session: windows of 4 x 4 plots per calculator

griswold = { ...
    options.doOridir,              'ndi.calc.vis.oridir_tuning',             'oridirtuning_calc'; ...
    options.doContrast,            'ndi.calc.vis.contrast_tuning',           'contrasttuning_calc'; ...
    options.doContrastSensitivity, 'ndi.calc.vis.contrast_sensitivity',      'contrastsensitivity_calc'; ...
    options.doSpatialFrequency,    'ndi.calc.vis.spatial_frequency_tuning',  'spatial_frequency_tuning_calc'; ...
    options.doTemporalFrequency,   'ndi.calc.vis.temporal_frequency_tuning', 'temporal_frequency_tuning_calc'};

if any([griswold{:,1}])
    [S, ~, q_cells] = ndi.calc.paperFigures.findElements(D, ...
        ndi.calc.paperFigures.exampleCells("griswold"));
    for k = find([griswold{:,1}])
        c = feval(griswold{k,2}, S);
        docs = S.database_search(ndi.query('','isa',griswold{k,3}) & q_cells);
        fh = cat(2, fh, ndi.calc.paperFigures.plotGrid(c, docs, 4, 4));
    end
end

% The Suarez Casanova session: one window per plot

suarez = { ...
    options.doSpeed,   'ndi.calc.vis.speed_tuning', 'speedtuning_calc'; ...
    options.doHartley, 'ndi.calc.vis.hartley',      'hartley_calc'};

if any([suarez{:,1}])
    [S, ~, q_cells] = ndi.calc.paperFigures.findElements(D, ...
        ndi.calc.paperFigures.exampleCells("suarezCasanova"));
    for k = find([suarez{:,1}])
        c = feval(suarez{k,2}, S);
        docs = S.database_search(ndi.query('','isa',suarez{k,3}) & q_cells);
        if strcmp(suarez{k,2}, 'ndi.calc.vis.hartley')
            % a cell has one Hartley document per Hartley epoch; show the last
            docs = ndi.calc.paperFigures.lastHartley(S, docs);
        end
        for i = 1:numel(docs)
            fh(end+1) = figure; %#ok<AGROW>
            c.plot(docs{i});
        end
    end
end

end
