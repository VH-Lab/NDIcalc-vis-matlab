function report = progressReporter(title, label)
% PROGRESSREPORTER - a progress bar that is always safe to call
%
%  REPORT = ndi.fun.progressReporter(TITLE, LABEL)
%
%  Returns a function handle REPORT(FRACTION) that advances a bar labeled
%  LABEL in an ndi.gui.component.ProgressBarWindow titled TITLE. FRACTION
%  runs from 0 to 1; REPORT(1) completes the bar, which then closes itself.
%  The window also estimates the time remaining from the rate of progress.
%
%  A long calculation that prints nothing is hard to tell apart from one that
%  has hung. This makes the difference visible.
%
%  No window is created when MATLAB cannot show figures -- started with
%  -batch or -nodisplay, as in continuous integration -- and REPORT is then a
%  handle that does nothing. Reporting also never raises: if the window cannot
%  be created the calculation carries on without it, after a warning saying so,
%  and if it cannot be updated -- because the viewer closed it -- the rest of
%  the calculation proceeds quietly. Callers never need to test anything before
%  calling REPORT.
%
%  Example:
%     report = ndi.fun.progressReporter('Hartley','Epoch t00001');
%     report(0.5);
%     report(1);
%
%  See also: ndi.gui.component.ProgressBarWindow
%

arguments
	title (1,:) char = ''
	label (1,:) char = ''
end

app = [];
tag = '';

try
	if ~batchStartupOptionUsed() && feature('ShowFigureWindows'),
		app = ndi.gui.component.ProgressBarWindow(title);
		tag = did.ido.unique_id();

		 % A bar that goes un-updated for longer than the window's timeout is
		 % flagged as timed out, and an 'Auto' bar is then closed. The stages
		 % reported here each run for minutes, so the bar must not be automatic
		 % and the timeout has to be longer than the gaps between updates. The
		 % bar is removed explicitly when the work reaches 1.
		app.setTimeout(minutes(60));
		app.addBar('Label',label,'Tag',tag,'Auto',false);
	end;
catch ME
	 % A missing display, no GUI toolkit, an older NDI without the component:
	 % none of these are reasons to stop the calculation. Say so rather than
	 % failing silently, because a reporter that quietly does nothing is
	 % indistinguishable from a calculation that is not reporting progress.
	app = [];
	warning('NDIcalc:progressReporter:unavailable',...
		'Could not open a progress window (%s). Continuing without one.',ME.message);
end;

report = @(fraction) ndi_fun_progressReporter_update(app,tag,fraction);

end % progressReporter()

function ndi_fun_progressReporter_update(app,tag,fraction)

if isempty(app),
	return;
end;

fraction = max(0,min(1,fraction));

try
	app.updateBar(tag,fraction);
	if fraction>=1,
		 % The bar is not automatic, so it is closed here rather than on a
		 % timeout that would otherwise fire during a long stage.
		app.removeBar(tag);
	end;
catch
	 % The user may have closed the window mid-calculation.
end;

end
