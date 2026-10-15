function [results, docs] = recomputeResults(D, options)
% RECOMPUTERESULTS - Recompute the paper's example results and compare them with the stored ones
%
% [RESULTS, DOCS] = NDI.CALC.PAPERFIGURES.RECOMPUTERESULTS(D, 'doX', true, ...)
%
% Runs the calculators again on the example cells of the calculators paper's
% example dataset, and compares each newly computed output document with the
% one stored in the dataset. Each calculator has its own section, switched on
% by its own name/value pair; all are off by default, unlike in
% ndi.calc.paperFigures.plotResults, because recomputing is slow (Hartley
% especially), so name the ones to recompute.
%
% For each calculator, the calculator's own input search
% (SEARCH_FOR_INPUT_PARAMETERS) lists every set of inputs it can run on; the
% sets whose stored output belongs to one of the example cells
% (ndi.calc.paperFigures.exampleCells) are calculated again with CALCULATE.
% The new documents are not added to the dataset, so the stored results are
% left as they are.
%
% Each new document is matched with the stored document that has the same
% dependencies and compared with the calculator's own comparison: the
% field-by-field tolerances of its first self-test
% (mock/<calculator>/mock.1.compare.json). Those tolerances allow for the
% fits' unseeded random start points, so a recomputed result is expected to
% agree within them rather than digit for digit. For ndi.calc.vis.hartley,
% the spike-triggered averages themselves are also compared.
%
% D is the ndi.dataset to use. If it is empty, the dataset is opened (and
% downloaded if necessary) with ndi.calc.paperFigures.openDataset. The
% calculators read the recorded data, so with file syncing off those files
% are downloaded on demand the first time they are needed.
%
% RESULTS is a table with one row per stored document compared:
%   calculator - the calculator class
%   element    - the cell's element name and reference
%   stored_id  - the id of the stored document
%   match      - true if the recomputed document agrees within tolerance
%   report     - where they disagree, if they do
%
% DOCS is a cell array of the recomputed ndi.documents, one per row of
% RESULTS (empty where no recomputed document matched the stored one). They
% can be drawn with the calculator's PLOT method, for example
%   c = ndi.calc.vis.oridir_tuning(S); figure; c.plot(docs{1});
% where S is the session the cell belongs to. A recomputed Hartley document
% keeps its spike-triggered average in a binary file that PLOT reads through
% the database, so it has to be added to the session before it can be
% plotted.
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
%   [results, docs] = ndi.calc.paperFigures.recomputeResults(D, 'doOridir', true)
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

results = table('Size',[0 5], ...
    'VariableTypes',{'string','string','string','logical','string'}, ...
    'VariableNames',{'calculator','element','stored_id','match','report'});
docs = {};

for k = find([sections{:,1}])
    cells = ndi.calc.paperFigures.exampleCells(sections{k,3});
    [S, elementIds] = ndi.calc.paperFigures.findElements(D, cells);
    elementNames = containers.Map(elementIds, ...
        arrayfun(@(i) sprintf('%s | %d', cells{i,1}, cells{i,2}), 1:size(cells,1), 'UniformOutput', false));

    c = feval(sections{k,2}, S);
    docCompare = c.load_mock_comparison(1);

    all_parameters = c.search_for_input_parameters(c.default_search_for_input_parameters());
    for i = 1:numel(all_parameters)
        stored = c.search_for_calculator_docs(all_parameters{i});
        stored = stored(cellfun(@(d) isKey(elementNames, d.dependency_value('element_id')), stored));
        if isempty(stored)
            continue; % not one of the example cells
        end
        recomputed = c.calculate(all_parameters{i});
        if ~iscell(recomputed)
            recomputed = {recomputed};
        end
        for j = 1:numel(stored)
            n = matchByDependencies(stored{j}, recomputed);
            if isempty(n)
                b = false;
                report = "No recomputed document has the same dependencies.";
                docs{end+1} = []; %#ok<AGROW>
            else
                [b, report] = compareDocs(c, stored{j}, recomputed{n}, docCompare);
                docs{end+1} = recomputed{n}; %#ok<AGROW>
            end
            results(end+1,:) = {string(sections{k,2}), ...
                string(elementNames(stored{j}.dependency_value('element_id'))), ...
                string(stored{j}.id()), logical(b), string(report)}; %#ok<AGROW>
        end
    end
end

end

function n = matchByDependencies(stored, recomputed)
% MATCHBYDEPENDENCIES - Find the recomputed document with the same dependencies as STORED
n = [];
key = dependencyKey(stored);
for i = 1:numel(recomputed)
    if strcmp(dependencyKey(recomputed{i}), key)
        n = i;
        return;
    end
end
end

function key = dependencyKey(doc)
% DEPENDENCYKEY - A string listing a document's dependencies in name order
deps = doc.document_properties.depends_on;
pairs = arrayfun(@(d) [char(d.name) '=' char(d.value)], deps(:), 'UniformOutput', false);
key = strjoin(sort(pairs), ';');
end

function [b, report] = compareDocs(c, stored, recomputed, docCompare)
% COMPAREDOCS - Compare a recomputed document with the stored one
%
% Uses the calculator's self-test tolerances (DOCCOMPARE); for Hartley
% documents, also compares the spike-triggered averages.

b = true;
report = "";

if ~isempty(docCompare)
    recomputed_struct = vlt.data.columnize_struct(recomputed.document_properties, 'columnizeNumericVectors', true);
    [b, r] = docCompare.compare(ndi.document(recomputed_struct), stored);
    if ~b
        report = string(formattedDisplayText(r));
        return;
    end
end

if isa(c, 'ndi.calc.vis.hartley')
    % read_sta finds the binary file through the database, so add the
    % recomputed document for the moment it takes to read it
    recomputed = recomputed.set_session_id(c.session.id());
    c.session.database_add(recomputed);
    cleanup = onCleanup(@() c.session.database_rm(recomputed.id()));
    sta_stored = c.read_sta(stored);
    sta_new = c.read_sta(recomputed);
    clear cleanup;
    max_diff = max(abs(sta_stored(:) - sta_new(:)));
    scale = max(abs(sta_stored(:)));
    if ~isequal(size(sta_stored), size(sta_new)) || max_diff > 1e-6 * max(scale, eps)
        b = false;
        report = sprintf("Spike-triggered averages differ: max difference %g (max |STA| %g).", max_diff, scale);
    end
end

end
