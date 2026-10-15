function [fh, docs] = recomputeResults(D, options)
% RECOMPUTERESULTS - Recompute the paper's example results and plot them
%
% [FH, DOCS] = NDI.CALC.PAPERFIGURES.RECOMPUTERESULTS(D, 'doX', true, ...)
%
% Runs the calculators again on the example cells of the calculators paper's
% example dataset and plots what they compute, to show that the calculators
% run end to end on real recordings. Each calculator has its own section,
% switched on by its own name/value pair; all are off by default, unlike in
% ndi.calc.paperFigures.plotResults, because recomputing is slow (Hartley
% especially), so name the ones to recompute.
%
% For each calculator, the calculator's own input search
% (SEARCH_FOR_INPUT_PARAMETERS), narrowed to the example cells, lists the sets
% of inputs it can run on; the sets that produced the stored results of the
% example cells
% (ndi.calc.paperFigures.exampleCells) are calculated again with CALCULATE.
% A cell has one Hartley result per epoch of Hartley stimuli; only the last
% (ndi.calc.paperFigures.lastHartley) is recomputed, with the same input
% parameters it was computed with, which saves most of the time.
%
% The new documents are plotted the way ndi.calc.paperFigures.plotResults
% plots the stored ones, in windows named "Recomputed: <calculator>". They
% are not saved: the stored results are left as they are. (A Hartley
% document's plot reads its spike-triggered average through the database, so
% each recomputed Hartley document is added to the session only while it is
% plotted, and removed straight after.)
%
% The recomputed results are not compared with the stored ones. The fits
% start from random points, so a fit on real data can settle somewhere
% slightly different on each run; that the calculators recover known answers
% is what their self-tests show, against stored expectations and tolerances.
%
% D is the ndi.dataset to use. If it is empty, the dataset is opened (and
% downloaded if necessary) with ndi.calc.paperFigures.openDataset. The
% calculators read the recorded data, so with file syncing off those files
% are downloaded on demand the first time they are needed.
%
% FH is a vector of the handles of the windows plotted. DOCS is a cell array
% of the recomputed ndi.documents.
%
% This function also takes name/value pairs that modify its behavior:
% |----------------------------|----------------------------------------------|
% | Parameter (default)        | Description                                  |
% |----------------------------|----------------------------------------------|
% | doOridir (false)           | Recompute ndi.calc.vis.oridir_tuning.        |
% | doContrast (false)         | Recompute ndi.calc.vis.contrast_tuning.      |
% | doContrastSensitivity      | Recompute ndi.calc.vis.contrast_sensitivity. |
% |   (false)                  |                                              |
% | doSpatialFrequency (false) | Recompute                                    |
% |                            |   ndi.calc.vis.spatial_frequency_tuning.     |
% | doTemporalFrequency (false)| Recompute                                    |
% |                            |   ndi.calc.vis.temporal_frequency_tuning.    |
% | doSpeed (false)            | Recompute ndi.calc.vis.speed_tuning.         |
% | doHartley (false)          | Recompute ndi.calc.vis.hartley. This one     |
% |                            |   takes minutes per cell.                    |
% |----------------------------|----------------------------------------------|
%
% Example:
%   D = ndi.calc.paperFigures.openDataset();
%   fh = ndi.calc.paperFigures.recomputeResults(D, 'doOridir', true);
%
% See also: ndi.calc.paperFigures.plotResults, ndi.calc.paperFigures.exampleCells,
%   ndi.calculator/calculate

arguments
    D = []
    options.doOridir (1,1) logical = false
    options.doContrast (1,1) logical = false
    options.doContrastSensitivity (1,1) logical = false
    options.doSpatialFrequency (1,1) logical = false
    options.doTemporalFrequency (1,1) logical = false
    options.doSpeed (1,1) logical = false
    options.doHartley (1,1) logical = false
end

if isempty(D)
    D = ndi.calc.paperFigures.openDataset();
end

sections = { ...
    options.doOridir,              'ndi.calc.vis.oridir_tuning',             "griswold"; ...
    options.doContrast,            'ndi.calc.vis.contrast_tuning',           "griswold"; ...
    options.doContrastSensitivity, 'ndi.calc.vis.contrast_sensitivity',      "griswold"; ...
    options.doSpatialFrequency,    'ndi.calc.vis.spatial_frequency_tuning',  "griswold"; ...
    options.doTemporalFrequency,   'ndi.calc.vis.temporal_frequency_tuning', "griswold"; ...
    options.doSpeed,               'ndi.calc.vis.speed_tuning',              "suarezCasanova"; ...
    options.doHartley,             'ndi.calc.vis.hartley',                   "suarezCasanova"};

fh = gobjects(1,0);
docs = {};

for k = find([sections{:,1}])
    [S, elementIds] = ndi.calc.paperFigures.findElements(D, ...
        ndi.calc.paperFigures.exampleCells(sections{k,3}));
    c = feval(sections{k,2}, S);

    % recompute the input sets whose stored results belong to the example cells
    docs_here = {};
    if isa(c, 'ndi.calc.vis.hartley')
        all_parameters = lastHartleyParameters(S, elementIds);
    else
        all_parameters = c.search_for_input_parameters(exampleCellSearch(c, elementIds));
    end
    for i = 1:numel(all_parameters)
        stored = c.search_for_calculator_docs(all_parameters{i});
        if ~any(cellfun(@(d) ismember(d.dependency_value('element_id'), elementIds), stored))
            continue; % not one of the example cells
        end
        recomputed = c.calculate(all_parameters{i});
        if ~iscell(recomputed)
            recomputed = {recomputed};
        end
        docs_here = cat(2, docs_here, reshape(recomputed, 1, []));
    end
    docs = cat(2, docs, docs_here);

    name = "Recomputed: " + string(sections{k,2});
    if sections{k,3} == "griswold"
        fh = cat(2, fh, ndi.calc.paperFigures.plotGrid(c, docs_here, 4, 4, name));
    else
        for i = 1:numel(docs_here)
            fh(end+1) = figure('Name', name); %#ok<AGROW>
            plotUnsaved(c, docs_here{i});
        end
    end
end

end

function parameters = exampleCellSearch(c, elementIds)
% EXAMPLECELLSEARCH - The calculator's input search, narrowed to the example cells
%
% A calculator's default search can match far more than cells:
% contrast_sensitivity's offers every document in the session as a candidate
% element, so listing its inputs for the whole session takes a very long
% time. Each query in the default search is narrowed here before the search
% runs: a query for the element itself to the example cells' ids, and any
% other (such as a tuning curve, which depends on its element) to documents
% that depend on one of the example cells.
parameters = c.default_search_for_input_parameters();
q_ids = ndi.query.none();
q_depends = ndi.query.none();
for k = 1:numel(elementIds)
    q_ids = q_ids | ndi.query('base.id','exact_string',elementIds{k},'');
    q_depends = q_depends | ndi.query('','depends_on','element_id',elementIds{k});
end
for i = 1:numel(parameters.query)
    if strcmp(parameters.query(i).name, 'element_id')
        parameters.query(i).query = parameters.query(i).query & q_ids;
    else
        parameters.query(i).query = parameters.query(i).query & q_depends;
    end
end
end

function plotUnsaved(c, doc)
% PLOTUNSAVED - Plot a document that is not in the database
%
% A Hartley document's PLOT reads its spike-triggered average through the
% database, so the document is added for the moment it takes to plot.
if isa(c, 'ndi.calc.vis.hartley')
    doc = doc.set_session_id(c.session.id());
    c.session.database_add(doc);
    cleanup = onCleanup(@() c.session.database_rm(doc.id())); %#ok<NASGU>
end
c.plot(doc);
end

function all_parameters = lastHartleyParameters(S, elementIds)
% LASTHARTLEYPARAMETERS - Inputs that recompute each cell's last Hartley result
%
% For each cell's last stored hartley_calc document, the input parameters it
% was computed with and the element and stimulus presentation it depends on.
q = ndi.query.none();
for k = 1:numel(elementIds)
    q = q | ndi.query('','depends_on','element_id',elementIds{k});
end
stored = ndi.calc.paperFigures.lastHartley(S, ...
    S.database_search(ndi.query('','isa','hartley_calc') & q));
all_parameters = cell(1, numel(stored));
for i = 1:numel(stored)
    all_parameters{i} = struct( ...
        'input_parameters', stored{i}.document_properties.hartley_calc.input_parameters, ...
        'depends_on', struct( ...
            'name',  {'element_id', 'stimulus_presentation_id'}, ...
            'value', {stored{i}.dependency_value('element_id'), ...
                      stored{i}.dependency_value('stimulus_presentation_id')}));
end
end
