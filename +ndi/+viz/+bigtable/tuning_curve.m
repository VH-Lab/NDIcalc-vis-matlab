function [bigtable] = tuning_curve(S, element_name, element_ref, table_func, varargin)


did.datastructures.assign(varargin{:});

doc_struct = ndi.viz.bigtable.tuning_curve_struct(S,element_name,element_ref,varargin{:});

N = numel(doc_struct.best_fit_docs);

if N==0,
    error(['No matching tuning curves.'])
end;

for i=1:N,
	element_struct(i) = ndi.viz.bigtable.element_struct(S,doc_struct.element_docs{i});
end;


for i=1:N,
	best_struct(i) = feval(table_func,doc_struct.best_fit_docs{i},...
			doc_struct.best_stimulus_response_scalar_docs{i},...
			doc_struct.best_stimulus_response_scalar_docs{i}.document_properties.stimulus_response_scalar.response_type,...
			doc_struct.f1f0(i),...
			'property_name',property_name,'x_axis_name',x_axis_name,...
			'prefix',['best_' test_prefix '_']);
end;

indiv_structs = {};
for k=1:size(doc_struct.individual_fit_docs,2),
	for i=1:N,
		rt = doc_struct.individual_response_scalar_docs{i,k}.document_properties.stimulus_response_scalar.response_type;
		indiv_structs{k}(i) = feval(table_func,doc_struct.individual_fit_docs{i,k},...
				doc_struct.individual_response_scalar_docs{i,k},...
				rt,...
				doc_struct.f1f0(i),...
				'property_name',property_name,'x_axis_name',x_axis_name,...
				'prefix',[rt '_' test_prefix '_']);
	end;
end;

bigtable = [struct2table(element_struct) struct2table(best_struct)];

for k=1:numel(indiv_structs),
	bigtable = [bigtable struct2table(indiv_structs{k})];
end;
