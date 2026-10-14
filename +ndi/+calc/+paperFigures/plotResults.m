function fh = plotResults(D, options)
% PLOTRESULTS - Plot the calculator results shown in the paper from the example dataset
%
% FH = NDI.CALC.PAPERFIGURES.PLOTRESULTS()
% FH = NDI.CALC.PAPERFIGURES.PLOTRESULTS(D)
%
% Plots the calculator output documents that the calculators paper shows,
% drawn from the example dataset that accompanies the paper. Each calculator
% has its own section, which picks out the example cells for that calculator
% and draws each of their output documents with the calculator's own PLOT
% method.
%
% The orientation/direction, contrast, contrast sensitivity, spatial
% frequency and temporal frequency examples come from the dataset's session
% from Griswold and Van Hooser (2025) and are drawn in windows of 4 x 4 plots.
% Until example cells are chosen for that session, every document of each
% type is drawn. The speed tuning and Hartley examples come from the session
% from the Suarez Casanova et al. study; their plots have several panels, so
% each gets its own window.
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
%   fh = ndi.calc.paperFigures.plotResults(D);
%
% See also: ndi.calc.paperFigures.openDataset, ndi.calculator/plot

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

fh = gobjects(0);

% The Griswold and Van Hooser session: one 4 x 4 grid of windows per calculator

griswold = { ...
    options.doOridir,              'ndi.calc.vis.oridir_tuning',             'oridirtuning_calc'; ...
    options.doContrast,            'ndi.calc.vis.contrast_tuning',           'contrasttuning_calc'; ...
    options.doContrastSensitivity, 'ndi.calc.vis.contrast_sensitivity',      'contrastsensitivity_calc'; ...
    options.doSpatialFrequency,    'ndi.calc.vis.spatial_frequency_tuning',  'spatial_frequency_tuning_calc'; ...
    options.doTemporalFrequency,   'ndi.calc.vis.temporal_frequency_tuning', 'temporal_frequency_tuning_calc'};

if any([griswold{:,1}])
    S = sessionByReference(D, griswoldSessionReference());
    q_cells = ndi.query.all();
    cells = griswoldCells();
    if ~isempty(cells)
        [~, elementIds] = findElements(D, cells);
        q_cells = dependsOnAny(elementIds);
    end
    for k = find([griswold{:,1}])
        c = feval(griswold{k,2}, S);
        docs = S.database_search(ndi.query('','isa',griswold{k,3}) & q_cells);
        fh = cat(2, fh, plotGrid(c, docs, 4, 4));
    end
end

% Speed tuning

if options.doSpeed
    % The same cells as the Hartley examples, from the Suarez Casanova session
    [S, elementIds] = findElements(D, suarezCasanovaCells());
    c = ndi.calc.vis.speed_tuning(S);
    docs = S.database_search(ndi.query('','isa','speedtuning_calc') & dependsOnAny(elementIds));
    for i = 1:numel(docs)
        fh(end+1) = figure; %#ok<AGROW>
        c.plot(docs{i});
    end
end

% Hartley reverse correlation

if options.doHartley
    [S, elementIds] = findElements(D, suarezCasanovaCells());
    c = ndi.calc.vis.hartley(S);
    docs = S.database_search(ndi.query('','isa','hartley_calc') & dependsOnAny(elementIds));
    for i = 1:numel(docs)
        fh(end+1) = figure; %#ok<AGROW>
        c.plot(docs{i});
    end
end

end

function ref = griswoldSessionReference()
% GRISWOLDSESSIONREFERENCE - Reference of the Griswold and Van Hooser session
ref = '2021-07-06';
end

function cells = griswoldCells()
% GRISWOLDCELLS - Example cells from the Griswold and Van Hooser session
%
% The cells, by element name and reference, whose orientation/direction,
% contrast, contrast sensitivity, spatial frequency and temporal frequency
% results are shown. Empty means every cell is shown.
cells = cell(0,2);
end

function cells = suarezCasanovaCells()
% SUAREZCASANOVACELLS - Example cells from the Suarez Casanova session
%
% The cells, by element name and reference, whose speed tuning and Hartley
% receptive fields are shown. They are recorded in the example dataset's
% session from the Suarez Casanova et al. study.
cells = { ...
    'leftcortex_11', 1; ...
    'leftcortex_22', 1; ...
    'leftcortex_26', 1; ...
    'leftcortex_31', 1; ...
    'leftcortex_35', 1};
end

function [S, elementIds] = findElements(D, cells)
% FINDELEMENTS - Look up the element documents of the example cells
%
% [S, ELEMENTIDS] = FINDELEMENTS(D, CELLS)
%
% CELLS is an N x 2 cell array of element names and references. Each must
% match exactly one element document in the dataset D, and all must belong to
% the same session. Returns that session S and the N element ids.

elementIds = cell(size(cells,1),1);
sessionIds = cell(size(cells,1),1);
for k = 1:size(cells,1)
    q = ndi.query('','isa','element') & ...
        ndi.query('element.name','exact_string',cells{k,1}) & ...
        ndi.query('element.reference','exact_number',cells{k,2});
    edoc = D.database_search(q);
    if numel(edoc) ~= 1
        error('NDICALC:paperFigures:plotResults:ElementNotFound', ...
            'Expected one element named %s | %d, found %d.', ...
            cells{k,1}, cells{k,2}, numel(edoc));
    end
    elementIds{k} = edoc{1}.id();
    sessionIds{k} = edoc{1}.document_properties.base.session_id;
end
sessionIds = unique(sessionIds);
if numel(sessionIds) ~= 1
    error('NDICALC:paperFigures:plotResults:MixedSessions', ...
        'The example cells come from %d sessions; expected one.', numel(sessionIds));
end
S = D.open_session(sessionIds{1});
end

function S = sessionByReference(D, ref)
% SESSIONBYREFERENCE - Open the session of the dataset D with reference REF
[refs, ids] = D.session_list();
match = find(strcmp(refs, ref));
if numel(match) ~= 1
    error('NDICALC:paperFigures:plotResults:SessionNotFound', ...
        'Expected one session with reference %s, found %d.', ref, numel(match));
end
S = D.open_session(ids{match});
end

function fh = plotGrid(c, docs, rows, cols)
% PLOTGRID - Plot documents in windows of ROWS x COLS axes
%
% FH = PLOTGRID(C, DOCS, ROWS, COLS)
%
% Draws each document in the cell array DOCS with the PLOT method of
% calculator C, one per axes, filling windows of ROWS x COLS axes with
% vlt.plot.supersubplot. Returns the handles of the windows.
fh = gobjects(1,0);
if isempty(docs)
    return;
end
fig = figure;
for i = 1:numel(docs)
    vlt.plot.supersubplot(fig, rows, cols, i);
    c.plot(docs{i});
end
fh = reshape(get(fig, 'userdata'), 1, []);
end

function q = dependsOnAny(elementIds)
% DEPENDSONANY - A query matching documents that depend on any of ELEMENTIDS
q = ndi.query.none();
for k = 1:numel(elementIds)
    q = q | ndi.query('','depends_on','element_id',elementIds{k});
end
end
