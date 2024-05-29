function out = bayes_pref(S, ndi_neuron_id)
%
%
%
%

ele_q = ndi.query('','depends_on','element_id',ndi_neuron_id);

q1 = ndi.query('','isa','stimulus_tuningcurve','');
q2 = ndi.query('stimulus_tuningcurve.independent_variable_label','contains_string','irection');
q3 = ndi.query('stimulus_tuningcurve.independent_variable_label','contains_string','eal');
q_total = q1 & q2 & q3;

oribayes_q = ndi.query('','isa','oridir_bayes_calc');

bayes_docs = S.database_search(ele_q&oribayes_q);

tcs = S.database_search(ele_q&q_total);

sf_tf_list = [];
isf1 = [];

bayes_curves = {};

srs = {};

for i=1:numel(tcs),

	q_tc = ndi.query('','depends_on','stimulus_tuningcurve_id',tcs{i}.document_properties.base.id); 
	srs_q = ndi.query('base.id','exact_string',tcs{i}.dependency_value('stimulus_response_scalar_id'));

	srs{i} = S.database_search(srs_q);
	srs{i} = srs{i}{1};

	isf1(i) = strcmp(srs{i}.document_properties.stimulus_response_scalar.response_type,'F1');

	% kludge for grant deadline

	sf_tf_list(i,[1 2]) = tcs{i}.document_properties.tuningcurve_calc.stim_property_list.values(:)';

	bayes_curves{i} = S.database_search(oribayes_q&q_tc); %q_ele_error&oribayes_q);
	if ~isempty(bayes_curves{i}),
		bayes_curves{i} = bayes_curves{i}{1};
	end;

end;

sf_tf_unique = unique(sf_tf_list,'rows');
sf_list = unique(sf_tf_unique(:,1));
tf_list = unique(sf_tf_unique(:,2));

bayes_sftf{1} = cell(numel(sf_list),numel(tf_list));
bayes_sftf{2} = cell(numel(sf_list),numel(tf_list));

resp_sftf{1} = cell(numel(sf_list),numel(tf_list));
resp_sftf{2} = cell(numel(sf_list),numel(tf_list));

for i=1:size(sf_tf_unique,1),
	z = find( sf_tf_list(:,1)==sf_tf_unique(i,1) & sf_tf_list(:,2)==sf_tf_unique(i,2) );
	sf_list_index = find(sf_list==sf_tf_unique(i,1));
	tf_list_index = find(tf_list==sf_tf_unique(i,2));
	for j=1:numel(z),
		bayes_sftf{1+isf1(z(j))}{sf_list_index,tf_list_index} = bayes_curves{z(j)};
		 resp_sftf{1+isf1(z(j))}{sf_list_index,tf_list_index} = tcs{z(j)};
	end;
end;

out.sf_tf_unique = sf_tf_unique;
out.bayes_sftf = bayes_sftf;
out.resp_sftf = resp_sftf;

