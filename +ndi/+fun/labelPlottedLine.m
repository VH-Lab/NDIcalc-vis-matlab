function h = labelPlottedLine(h, name)
% LABELPLOTTEDLINE - name a plotted line so a reader can identify it
%
%  H = ndi.fun.labelPlottedLine(H, NAME)
%
%  Gives the line H the display name NAME. The name appears in a legend, if
%  the caller draws one, and is added as a row of the line's data tip, so that
%  clicking or hovering over the curve names it.
%
%  A plot that draws several fits of the same data over one another is
%  otherwise unlabelled, and a reader cannot tell which curve came from which
%  fit. Naming the lines lets the figure answer that question itself.
%
%  H is returned with only these properties changed. This function never
%  raises: labelling is a convenience, and a caller drawing into an unusual
%  graphics context should not fail because a data tip could not be attached.
%
%  Example:
%     h = plot(x,y,'b-');
%     ndi.fun.labelPlottedLine(h,'Movshon et al. 2005 fit');
%     legend show;   % or hover over the curve
%
%  See also: dataTipTextRow, legend
%

arguments
	h
	name (1,:) char
end

try
	set(h,'DisplayName',name);

	 % A constant data-tip row needs one entry per plotted point.
	if isprop(h,'DataTipTemplate'),
		x = get(h,'XData');
		h.DataTipTemplate.DataTipRows(end+1) = ...
			dataTipTextRow('fit',repmat({name},size(x)));
	end;
catch
	 % Older graphics, a context without data tips, or a deleted handle.
end;
