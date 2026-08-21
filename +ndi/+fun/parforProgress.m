function tick = parforProgress(total, progressFcn, numberOfWorkers, options)
% PARFORPROGRESS - report how much of a parfor loop has finished
%
%  TICK = ndi.fun.parforProgress(TOTAL, PROGRESSFCN, NUMBEROFWORKERS, ...)
%
%  Returns a handle TICK to be called once per iteration, at the end of a
%  parfor body running TOTAL iterations. PROGRESSFCN is then called in the
%  client as PROGRESSFCN(FRACTION), with FRACTION the proportion of iterations
%  that have completed.
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
%        loop (default 100, i.e. about once per percent). The last iteration
%        always reports, so PROGRESSFCN(1) is always reached.
%
%  PROGRESSFCN may be empty, in which case TICK does nothing. Reporting never
%  raises: a reporter that fails is not a reason to abandon the calculation.
%
%  Example:
%     n = ndi.fun.parallelWorkers();
%     tick = ndi.fun.parforProgress(N, @(f) disp([num2str(round(100*f)) '%']), n);
%     parfor (i=1:N, n)
%        % ... work ...
%        tick();
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
	tick = @() [];
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
		return;
	end;
end;

tick = @() ndi_fun_parforProgress_count(counter,total,everyN,progressFcn);

end % parforProgress()

function ndi_fun_parforProgress_count(counter,total,everyN,progressFcn)

n = counter('done') + 1;
counter('done') = n;

if n<total && mod(n,everyN)~=0,
	return;
end;

try
	progressFcn(min(1,n/total));
catch
	 % The viewer closed the window, or the reporter is otherwise gone.
end;

end
