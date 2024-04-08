function T_out = stimtunefitT2features(stimtunefit_table, tablefunc)


for i=1:size(stimtunefit_table,1),
	mystruct = feval(tablefunc,...
			stimtunefit_table.fit_doc(i),...
			stimtunefit_table.stim_resp_doc(i),...
			stimtunefit_table.stim_resp_type(i),...
			[]);
	if i==1,
		mys = mystruct;
	else,
		mys(end+1) = mystruct;
	end;
end;

t_inc = struct2table(mys);

T_out = [stimtunefit_table t_inc];


