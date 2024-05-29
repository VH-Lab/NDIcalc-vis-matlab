function [t_stimtunefit] = stimuluspurpose2tuning(S, purpose, purpose_short, independent_var, tuning_class, acceptable_response_types)
% STIMULUSPURPOSE2TUNING - get stimulus tuning and fit docs from stimulus purpose
%
% [TUNING_DOCS, FIT_DOCS, STIM_RESP_DOCS, RESP_TYPE] = STIMULUSPURPOSE2TUNING(S, PURPOSE, PURPOSE_SHORT, INDEPENDENT_VAR, TUNING_CLASS, ACCEPTABLE_RESPONSE_TYPES)
%
% Returns cell arrays of tuning documents, fit documents, stimulus response docs, and the text representation of response type.
% The ith entry of each cell array corresponds to the ith array of each other cell array.
% 

   % step 1: find the right tuning curve documents

q = ndi.viz.fun.stimuluspurpose2stimulusresponsequery(S, purpose);

i = 1; 
q_rt = ndi.query('stimulus_response_scalar.response_type','exact_string',acceptable_response_types{i}); % reducing cut and paste errors with the i here
for i=2:numel(acceptable_response_types),
	q_rt = q_rt |  ndi.query('stimulus_response_scalar.response_types','exact_string',acceptable_response_types{i});
end;
stim_response_docs = S.database_search(q & q_rt);

q2 = ndi.query('stimulus_tuningcurve.independent_variable_label','contains_string',independent_var);

if numel(stim_response_docs)>0,
	i = 1;
	q3 = ndi.query('','depends_on','stimulus_response_scalar_id',stim_response_docs{i}.document_properties.base.id);
	
	for i=2:numel(stim_response_docs),
		q3 = q3 | ndi.query('','depends_on','stimulus_response_scalar_id',stim_response_docs{i}.document_properties.base.id);
	end;
else,
	q3 = ndi.query('','isa','nothing');
end;

tuning_docs = S.database_search(q2&q3);

if numel(tuning_docs)>0,
	i = 1;
	q4 = ndi.query('','depends_on','stimulus_tuningcurve_id',tuning_docs{i}.document_properties.base.id);
	for i=2:numel(tuning_docs),
		q4 = q4 | ndi.query('','depends_on','stimulus_tuningcurve_id',tuning_docs{i}.document_properties.base.id);
	end;
else,
	q4 = ndi.query('','isa','nothing');
end;

q5 = ndi.query('','isa',tuning_class);

fit_docs = S.database_search(q4&q5);

 % okay, now sort into acceptable_response_type values

stim_r_struct = did.datastructures.emptystruct('stim_resp_id','stim_resp_element_id','stim_resp_type','stim_epoch','stim_resp_doc_index','stim_resp_doc');
for i=1:numel(stim_response_docs),
	stim_r_struct(end+1) = struct('stim_resp_id',stim_response_docs{i}.document_properties.base.id,...
		'stim_resp_element_id',stim_response_docs{i}.dependency_value('element_id'),...
		'stim_resp_type',stim_response_docs{i}.document_properties.stimulus_response_scalar.response_type,...
		'stim_epoch',stim_response_docs{i}.document_properties.stimulus_response.stimulator_epochid,...
		'stim_resp_doc_index',i,...
		'stim_resp_doc',stim_response_docs{i});
end;

stim_r_table = struct2table(stim_r_struct,'AsArray',true);

tuning_struct = did.datastructures.emptystruct('tune_doc_id','tune_stim_resp_dependency','tune_doc_index','tune_doc');
for i=1:numel(tuning_docs),
	tuning_struct(end+1) = struct('tune_doc_id',tuning_docs{i}.document_properties.base.id,...
		'tune_stim_resp_dependency',tuning_docs{i}.dependency_value('stimulus_response_scalar_id'),...
		'tune_doc_index',i,...
		'tune_doc', tuning_docs{i});
end;

tuning_table = struct2table(tuning_struct,'AsArray',true);

fit_struct = did.datastructures.emptystruct('fit_doc_id','fit_tuning_doc_dependency','fit_doc_index','fit_doc');
for i=1:numel(fit_docs),
	fit_struct(end+1) = struct('fit_doc_id',fit_docs{i}.document_properties.base.id,...
		'fit_tuning_doc_dependency',fit_docs{i}.dependency_value('stimulus_tuningcurve_id'),...
		'fit_doc_index',i,...
		'fit_doc',fit_docs{i});
end;

fit_table = struct2table(fit_struct,'AsArray',true);

t_stimtune = join(stim_r_table,tuning_table,'LeftKeys','stim_resp_id','RightKeys','tune_stim_resp_dependency');

t_stimtunefit = join(t_stimtune,fit_table,'LeftKeys','tune_doc_id','RightKeys','fit_tuning_doc_dependency');

P = ones(size(t_stimtunefit,1),numel(purpose_short));
t_here = array2table(P,'VariableNames',purpose_short);

t_stimtunefit = [t_stimtunefit t_here];

