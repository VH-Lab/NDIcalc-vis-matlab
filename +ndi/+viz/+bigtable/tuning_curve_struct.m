function [doc_struct] = tuning_curve_struct(S, element_name, element_ref, varargin)


% This function can be modified by the following name/value pairs:
% |--------------------------------|------------------------------------------|
% | Parameters (default)           | Description                              |
% |--------------------------------|------------------------------------------|
% | tuning_curve_search_string     | `independent_variable_label`  search     |
% |  ('Temporal Frequency')        |   string of `stimulus_tuningcurve`       |
% | fit_doc_type ('')              | Doc type to search for (e.g.,            |
% |                                |  `temporal_frequency_tuning_calc`)       |
% | element_type ('spikes')        | Element type?                            |
% | search_query ([])              | An additional ndi.query for the search   |
% | acceptable_response_types      | Response types to include in table       |
% |    ({'mean','F1'})             |  (e.g., 'mean', 'F1', 'F2')              |
% |--------------------------------|------------------------------------------|
%
%

tuning_curve_search_string = 'Temporal Frequency';
fit_doc_type = 'temporal_frequency_tuning_calc';
fit_type_property = 'temporal_frequency_tuning';
element_type = 'spikes';
search_query = [];
acceptable_response_types = {'mean','F1'};

did.datastructures.assign(varargin{:});

q1 = ndi.query('element.type','exact_string',element_type);
q2 = ndi.query('element.name','contains_string',element_name,'');
q3 = ndi.query('element.reference','exact_number',element_ref);
element_docs = S.database_search(q1&q2&q3);

if isempty(element_docs),
	S,element_name,element_ref
	error(['Empty element documents, no elements found with those parameters.']);
end;

q_o = ndi.query('','isa',fit_doc_type);
if ~isempty(search_query),
	q_o = q_o & search_query;
end;
all_fit_type = S.database_search(q_o);

q_otc = ndi.query('stimulus_tuningcurve.independent_variable_label','contains_string',tuning_curve_search_string);
all_tuning_curves = S.database_search(q_otc); % these are the stimulus tuning curves that form basis of responses

matching_calc = [];
matching_stim_response_doc = [];

for i=1:numel(all_fit_type),
	has_match = 0;
	for j=1:numel(all_tuning_curves),
		if strcmp(all_fit_type{i}.dependency_value('stimulus_tuningcurve_id'),...
			all_tuning_curves{j}.document_properties.base.id),
			% if we match, mark it
			[response_type,stim_response_scalar_doc] = ndi.fun.stimulus.tuning_curve_to_response_type(S,...
				all_tuning_curves{j});
			index = find(strcmp(response_type,acceptable_response_types));
			if ~isempty(index),
				has_match = 1;
				break;
			end;
		end;
	end;
	if has_match,
		matching_calc(end+1) = i; % i has a match and here is its stim_response_scalar
		matching_stim_response_doc{end+1} = stim_response_scalar_doc;
	end;
end;

 % now, we can limit ourselves to the matching ones

all_fit_type = all_fit_type(matching_calc);
 % and matching_stim_response_doc matches each one on a 1-1 basis

fit_docs = {};
best_stimulus_response_scalar_docs = {};
individual_fit_docs = {};
individual_stimulus_response_scalar_docs = {};

for i=1:numel(element_docs),
	bestValue = -Inf;
	bestValueIndex = 0;
	% find matches, and one with best response
	for j=1:numel(all_fit_type),
		if strcmp(all_fit_type{j}.dependency_value('element_id'),element_docs{i}.document_properties.base.id),
			info_here = getfield(all_fit_type{j}.document_properties,fit_type_property);
			[response_type,stim_response_scalar_doc_here] = ndi.fun.stimulus.tuning_curve_to_response_type(S,...
				all_fit_type{j});
			index = find(strcmp(response_type,acceptable_response_types));
			if ~isempty(index),
				individual_fit_docs{i,index} = all_fit_type{j};
				individual_stimulus_response_scalar_docs{i,index} = stim_response_scalar_doc_here;
			end;
			maxValue = nanmax(info_here.tuning_curve.mean);
			if maxValue > bestValue,
				bestValue = maxValue;
				bestValueIndex = j;
			end;
		end;
	end;
	if bestValueIndex > 0,
		fit_docs{i} = all_fit_type{bestValueIndex};
		best_stimulus_response_scalar_docs{i} = matching_stim_response_doc{bestValueIndex};
	end;
end;

include = [];
for i=1:numel(fit_docs),
        if ~isempty(fit_docs{i}),
                include(end+1) = i;
        end;
end;

fit_docs = fit_docs(include);
best_stimulus_response_scalar_docs = best_stimulus_response_scalar_docs(include);
individual_fit_docs = individual_fit_docs(include,:);
individual_stimulus_response_scalar_docs = individual_stimulus_response_scalar_docs(include,:);
element_docs = element_docs(include);

info = {};

f1f0 = vlt.data.emptystruct('f0','f1','f1f0_ratio','f1f0_2f1overf1f0');

for i=1:numel(fit_docs)
	clear f1f0_here;
	info{i} = getfield(fit_docs{i}.document_properties,fit_type_property);
	[f1f0_here.f0,f1f0_here.f1] = ndi.fun.stimulus.f0_f1_responses(S,fit_docs{i});
	f1f0_here.f1f0_ratio = f1f0_here.f1/f1f0_here.f0;
	f1f0_here.f1f0_2f1overf1f0 = 2*(f1f0_here.f1)/(f1f0_here.f0+f1f0_here.f1);
	f1f0(i) = f1f0_here;
end;

  % needs to have 

doc_struct.best_fit_docs = fit_docs;
doc_struct.best_stimulus_response_scalar_docs = best_stimulus_response_scalar_docs;

doc_struct.individual_response_scalar_docs = individual_stimulus_response_scalar_docs;
doc_struct.individual_fit_docs = individual_fit_docs;
doc_struct.best_stimulus_response_scalar_docs = best_stimulus_response_scalar_docs;
doc_struct.f1f0 = f1f0;
doc_struct.element_docs = element_docs;
doc_struct.options.fit_doc_type = fit_doc_type;
doc_struct.options.fit_type_property = fit_type_property;
doc_struct.options.tuning_curve_search_string = tuning_curve_search_string;
doc_struct.options.element_type = element_type ;
doc_struct.options.search_query = search_query;
doc_struct.options.acceptable_response_types = acceptable_response_types;


