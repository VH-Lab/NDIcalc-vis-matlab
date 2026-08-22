function tick = parforProgress(total, progressFcn, numberOfWorkers, options)
% PARFORPROGRESS - report how much of a parfor loop has finished
%
%  TICK = ndi.fun.parforProgress(TOTAL, PROGRESSFCN, NUMBEROFWORKERS, ...)
%
%  Returns a handle TICK(I) to be called once per iteration, at the end of a
%  parfor body running TOTAL iterations, with I the loop index. PROGRESSFCN is
%  then called in the client as PROGRESSFCN(FRACTION), with FRACTION the
%  proportion of iterations that have completed.
%
%  TICK does nothing on most iterations. It acts only on every (TOTAL /
%  UpdateSteps)th index and on the last, so the work it costs the loop is
%  bounded by UpdateSteps however many iterations there are. This matters most
%  when the loop has workers: there TICK sends a message from the worker to the
%  client, which is far more expensive than a counter increment and would
%  otherwise happen on every one of what may be tens of thousands of iterations.
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
%     tick = ndi.fun.parforProgress(N, @(f) disp([num2str(round(100*f)) '%']), n);
%     parfor (i=1:N, n)
%        % ... work ...
%        tick(i);
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
	tick = @(i) [];
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
		tick = @(i) ndi_fun_parforProgress_send(queue,i,everyN,total);
		return;
	catch ME
		 % No Parallel Computing Toolbox, or no pool to attach to. The loop
		 % still runs; only the reporting is lost, and silently losing it is
		 % what this function exists to prevent.
		warning('NDIcalc:parforProgress:unavailable',...
			'Could not open a progress queue (%s). Continuing without progress reports.',ME.message);
		tick = @(i) [];
		return;
	end;
end;

tick = @(i) ndi_fun_parforProgress_serial(counter,i,everyN,total,progressFcn);

end % parforProgress()

function ndi_fun_parforProgress_send(queue,i,everyN,total)
 % Runs on the worker, so it must be as close to free as possible on the
 % iterations that do not report: one mod, then return.

if mod(i,everyN)~=0 && i~=total,
	return;
end;

send(queue,1);

end % ndi_fun_parforProgress_send()

function ndi_fun_parforProgress_serial(counter,i,everyN,total,progressFcn)
 % The no-worker path. Same gate, then report directly rather than through a
 % queue, since there is no worker to send from.

if mod(i,everyN)~=0 && i~=total,
	return;
end;

ndi_fun_parforProgress_count(counter,total,everyN,progressFcn);

end % ndi_fun_parforProgress_serial()

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
