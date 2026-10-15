function fh = plotGrid(c, docs, rows, cols, name)
% PLOTGRID - Plot calculator documents in windows of ROWS x COLS axes
%
% FH = NDI.CALC.PAPERFIGURES.PLOTGRID(C, DOCS, ROWS, COLS, [NAME])
%
% Draws each document in the cell array DOCS with the PLOT method of
% calculator C, one per axes, filling windows of ROWS x COLS axes with
% vlt.plot.supersubplot. NAME, if given, is used as the windows' name.
% Returns the handles of the windows.
%
% See also: ndi.calc.paperFigures.plotResults, ndi.calc.paperFigures.recomputeResults

arguments
    c (1,1) ndi.calculator
    docs cell
    rows (1,1) double
    cols (1,1) double
    name (1,1) string = ""
end

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
if strlength(name) > 0
    set(fh, 'Name', name);
end
end
