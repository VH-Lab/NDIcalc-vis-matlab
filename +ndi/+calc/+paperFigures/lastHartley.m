function docs = lastHartley(S, docs)
% LASTHARTLEY - Keep only each cell's last Hartley document
%
% DOCS = NDI.CALC.PAPERFIGURES.LASTHARTLEY(S, DOCS)
%
% A cell can have more than one hartley_calc document, one for each
% stimulus presentation (epoch) of Hartley stimuli it was recorded during.
% Given a cell array DOCS of hartley_calc documents from the ndi.session S,
% returns one document per cell (element): the one whose stimulus
% presentation has the last epoch id. The returned documents are in the
% order their cells first appear in DOCS.
%
% See also: ndi.calc.paperFigures.plotResults, ndi.calc.paperFigures.recomputeResults

elementIds = cellfun(@(d) d.dependency_value('element_id'), docs, 'UniformOutput', false);
epochIds = cell(size(docs));
for i = 1:numel(docs)
    sp = S.database_search(ndi.query('base.id','exact_string', ...
        docs{i}.dependency_value('stimulus_presentation_id')));
    if isempty(sp)
        epochIds{i} = '';
    else
        epochIds{i} = sp{1}.document_properties.epochid.epochid;
    end
end

[~, first] = unique(elementIds, 'stable');
keep = zeros(1, numel(first));
for k = 1:numel(first)
    here = find(strcmp(elementIds, elementIds{first(k)}));
    [~, order] = sort(epochIds(here));
    keep(k) = here(order(end));
end
docs = docs(keep);
end
