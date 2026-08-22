function [tick, everyN] = parforProgress(total, progressFcn, numberOfWorkers, options)
% PARFORPROGRESS - report how much of a parfor loop has finished
%
%  [TICK, EVERYN] = ndi.fun.parforProgress(TOTAL, PROGRESSFCN, NUMBEROFWORKERS, ...)
%
%  Returns a handle TICK, and the interval EVERYN at which the caller should
%  call it. Each call reports one step of progress; PROGRESSFCN is called in
%  the client as PROGRESSFCN(FRACTION), with FRACTION the proportion of the
%  loop finished. Call it from the end of a parfor body like this:
%
%     if mod(i,everyN)==0 || i==total
%        tick();
%     end
%
%  The interval is the point. Reporting on every iteration is wasteful, and
%  when the loop has workers it is expensive: TICK then sends a message from
%  the worker to the client, far dearer than a counter increment, and there may
%  be tens of thousands of iterations. Gating in the loop body means most
%  iterations cost one mod and no call at all.
%
%  The gate belongs at the call site rather than inside TICK because parfor
%  cannot tell a function handle called with the loop variable from an array
%  sliced by it: TICK(I) inside a parfor raises
%  MATLAB:parfor:sliced_function_handle. Calling TICK with no arguments is
%  unambiguous. (MATLAB suggests feval(tick,i) instead, which also works, but
%  still pays a call on every iteration.)
%
%  A parfor hands its iterations out in an arbitrary order, so the loop index
%  says nothing about how much of the work is done: index 9000 of 10000 may be
%  the third slice to finish. Counting completions instead gives a fraction
%  that only ever moves forward. Reporting from inside the loop also has to
%  reach the client to be seen at all -- a worker's own output is buffered
%  until the loop ends, which is why an uninstrumented parfor appears to do
%  nothing and then finish all at once.
%
%  NUMBEROFWORKERS is the worker limit the loop was given, as returned by
%  ndi.fun.parallelWorkers. At 0 the loop runs in the client and TICK reports
%  directly; above 0 it reports through a parallel.pool.DataQueue, which the
%  client drains while the loop runs.
%
%  Optional name-value pair:
%     UpdateSteps - how many times PROGRESSFCN should be called over the whole
%        loop (default 100, i.e. about once per percent). This also sets how
%        often TICK does any work at all. The last iteration always reports, so
%        PROGRESSFCN(1) is always reached.
%
%  Because only every Nth index reports, the fraction is the count of those
%  that have finished scaled by N rather than an exact tally. It still only
%  moves forward and still ends at 1; it is a progress report, not a census.
%
%  PROGRESSFCN may be empty, in which case TICK does nothing. Reporting never
%  raises: a reporter that fails is not a reason to abandon the calculation.
%
%  Example:
%     n = ndi.fun.parallelWorkers();
%     [tick,everyN] = ndi.fun.parforProgress(N, ...
%        @(f) disp([num2str(round(100*f)) '%']), n);
%     parfor (i=1:N, n)
%        % ... work ...
%        if mod(i,everyN)==0 || i==N, tick(); end;
%     end
%
%  See also: ndi.fun.parallelWorkers, ndi.fun.progressReporter
%

arguments
	total (1,1) double {mustBeNonnegative}
	progressFcn = []
	numberOfWorkers (1,1) double = 0
	options.UpdateSteps (1,1) double {mustBePositive} = 100
end

if isempty(progressFcn) || total==0,
	 % Nothing to report. everyN is still returned so the caller's gate is
	 % valid; at this value it lets through only the last iteration, which
	 % then calls a handle that does nothing.
	tick = @() [];
	everyN = max(1,total);
	return;
end;

everyN = max(1,floor(total/options.UpdateSteps));

 % A containers.Map is a handle, so the count survives in the closure below
 % rather than being captured by value the way an ordinary variable would be.
counter = containers.Map({'done'},{0});

if numberOfWorkers>0,
	try
		queue = parallel.pool.DataQueue;
		afterEach(queue,@(~) ndi_fun_parforProgress_count(counter,total,everyN,progressFcn));
		tick = @() send(queue,1);
		return;
	catch ME
		 % No Parallel Computing Toolbox, or no pool to attach to. The loop
		 % still runs; only the reporting is lost, and silently losing it is
		 % what this function exists to prevent.
		warning('NDIcalc:parforProgress:unavailable',...
			'Could not open a progress queue (%s). Continuing without progress reports.',ME.message);
		tick = @() [];
		everyN = max(1,total);
		return;
	end;
end;

tick = @() ndi_fun_parforProgress_count(counter,total,everyN,progressFcn);

end % parforProgress()

function ndi_fun_parforProgress_count(counter,total,everyN,progressFcn)
 % Runs in the client, once per reporting iteration. Each call stands for
 % everyN iterations of the loop, which is why the fraction is scaled rather
 % than counted; min keeps a remainder at the end from overshooting 1.

n = counter('done') + 1;
counter('done') = n;

try
	progressFcn(min(1,n*everyN/total));
catch
	 % The viewer closed the window, or the reporter is otherwise gone.
end;

end
