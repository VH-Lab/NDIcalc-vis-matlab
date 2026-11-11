function [S, elementIds, q] = findElements(D, cells)
% FINDELEMENTS - Look up the element documents of example cells
%
% [S, ELEMENTIDS, Q] = NDI.CALC.PAPERFIGURES.FINDELEMENTS(D, CELLS)
%
% CELLS is an N x 2 cell array of element names and references, such as
% the output of ndi.calc.paperFigures.exampleCells. Each must match exactly
% one element document in the ndi.dataset D, and all must belong to the same
% session.
%
% Returns that session S, the N element ids ELEMENTIDS, and an ndi.query Q
% that matches any document that depends on one of those elements through an
% 'element_id' dependency.
%
% See also: ndi.calc.paperFigures.exampleCells

elementIds = cell(size(cells,1),1);
sessionIds = cell(size(cells,1),1);
for k = 1:size(cells,1)
    qe = ndi.query('','isa','element') & ...
        ndi.query('element.name','exact_string',cells{k,1}) & ...
        ndi.query('element.reference','exact_number',cells{k,2});
    edoc = D.database_search(qe);
    if numel(edoc) ~= 1
        error('NDICALC:paperFigures:findElements:ElementNotFound', ...
            'Expected one element named %s | %d, found %d.', ...
            cells{k,1}, cells{k,2}, numel(edoc));
    end
    elementIds{k} = edoc{1}.id();
    sessionIds{k} = edoc{1}.document_properties.base.session_id;
end
sessionIds = unique(sessionIds);
if numel(sessionIds) ~= 1
    error('NDICALC:paperFigures:findElements:MixedSessions', ...
        'The cells come from %d sessions; expected one.', numel(sessionIds));
end
S = D.open_session(sessionIds{1});

q = ndi.query.none();
for k = 1:numel(elementIds)
    q = q | ndi.query('','depends_on','element_id',elementIds{k});
end

end
