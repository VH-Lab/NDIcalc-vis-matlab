function report = progressReporter(title, label)
% PROGRESSREPORTER - a progress bar that is always safe to call
%
%  REPORT = ndi.fun.progressReporter(TITLE, LABEL)
%
%  Returns a function handle REPORT(FRACTION) that advances a bar labelled
%  LABEL in an ndi.gui.component.ProgressBarWindow titled TITLE. FRACTION
%  runs from 0 to 1; REPORT(1) completes the bar, which then closes itself.
%  The window also estimates the time remaining from the rate of progress.
%
%  A long calculation that prints nothing is hard to tell apart from one that
%  has hung. This makes the difference visible.
%
%  No window is created when MATLAB cannot show figures -- started with
%  -batch or -nodisplay, as in continuous integration -- and REPORT is then a
%  handle that does nothing. Reporting also never raises: if the window
%  cannot be created or updated, the calculation carries on without it.
%  Callers therefore never need to test anything before calling REPORT.
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
		tag = did.ido.unique();
		app.addBar('Label',label,'Tag',tag,'Auto',true);
	end;
catch
	 % A missing display, no GUI toolkit, an older NDI without the component:
	 % none of these are reasons to stop the calculation.
	app = [];
end;

report = @(fraction) ndi_fun_progressReporter_update(app,tag,fraction);

end % progressReporter()

function ndi_fun_progressReporter_update(app,tag,fraction)

if isempty(app),
	return;
end;

try
	app.updateBar(tag,max(0,min(1,fraction)));
catch
	 % The user may have closed the window mid-calculation.
end;

end
