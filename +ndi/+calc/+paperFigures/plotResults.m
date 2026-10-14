function fh = plotResults(D, options)
% PLOTRESULTS - Plot the calculator results shown in the paper from the example dataset
%
% FH = NDI.CALC.PAPERFIGURES.PLOTRESULTS()
% FH = NDI.CALC.PAPERFIGURES.PLOTRESULTS(D)
%
% Plots the calculator output documents that the calculators paper shows,
% drawn from the example dataset that accompanies the paper. Each calculator
% has its own section, which picks out the example cells for that calculator
% and draws each of their output documents in its own figure with the
% calculator's own PLOT method.
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
    options.doSpeed (1,1) logical = true
    options.doHartley (1,1) logical = true
end

if isempty(D)
    D = ndi.calc.paperFigures.openDataset();
end

fh = gobjects(0);

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

function cells = suarezCasanovaCells()
% SUAREZCASANOVACELLS - Example cells from the Suarez Casanova session
%
% The cells, by element name and reference, whose speed tuning and Hartley
% receptive fields are shown. They are recorded in the example dataset's
% session from the Suarez Casanova et al. study.
cells = { ...
    'leftcortex_8',  1; ...
    'leftcortex_11', 1; ...
    'leftcortex_22', 1; ...
    'leftcortex_26', 1; ...
    'leftcortex_31', 1; ...
    'leftcortex_35', 1; ...
    'leftcortex_37', 1};
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

function q = dependsOnAny(elementIds)
% DEPENDSONANY - A query matching documents that depend on any of ELEMENTIDS
q = ndi.query.none();
for k = 1:numel(elementIds)
    q = q | ndi.query('','depends_on','element_id',elementIds{k});
end
end
