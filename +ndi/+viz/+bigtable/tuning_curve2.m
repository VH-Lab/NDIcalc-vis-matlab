function [bigtable] = tuning_curve2(S, element_name, element_ref, table_func, varargin)


did.datastructures.assign(varargin{:});

doc_struct = ndi.viz.bigtable.tuning_curve_struct(S,element_name,element_ref,varargin{:});

N = numel(doc_struct.best_fit_docs);

if N==0,
    error(['No matching tuning curves.'])
end;

for i=1:N,
	element_struct(i) = ndi.viz.bigtable.element_struct(S,doc_struct.element_docs{i});
end;

best_table = [];

for i=1:N,
	best_table_here = feval(table_func,doc_struct.best_fit_docs{i},...
			doc_struct.best_stimulus_response_scalar_docs{i},...
			doc_struct.best_stimulus_response_scalar_docs{i}.document_properties.stimulus_response_scalar.response_type,...
			doc_struct.f1f0(i),...
			'property_name',property_name,'x_axis_name',x_axis_name,...
			'prefix',['best_' test_prefix '_']);
	best_table = [best_table; best_table_here];
end;

indiv_table = [];

for i=1:N,
	indiv_table_here = [];
	for k=1:size(doc_struct.individual_fit_docs,2),
		rt = doc_struct.individual_response_scalar_docs{i,k}.document_properties.stimulus_response_scalar.response_type;
		table_here = feval(table_func,doc_struct.individual_fit_docs{i,k},...
				doc_struct.individual_response_scalar_docs{i,k},...
				rt,...
				doc_struct.f1f0(i),...
				'property_name',property_name,'x_axis_name',x_axis_name,...
				'prefix',[rt '_' test_prefix '_']);
		indiv_table_here = [indiv_table_here table_here]; % cat(2)
	end;
	indiv_table = [indiv_table; indiv_table_here]; % cat(1)
end;

bigtable = [struct2table(element_struct) best_table indiv_table ];


