function T = element_subject(S, element_type)
% ELEMENT_SUBJECT - return elements and subject information for a session/dataset
%
% T = ELEMENT_SUBJECT(S, ELEMENT_TYPE)
% 


q_e = ndi.query('element.type','exact_string',element_type);
el = S.database_search(q_e);

q_s = ndi.query('','isa','subject');
s = S.database_search(q_s);


e_struct = did.datastructures.emptystruct('element_full_name','element_name','element_reference',...
	'element_class','element_doc_id','element_type','element_doc_id','element_session_id','element_subject_id');
for i=1:numel(el),
	e = el{i}.document_properties; % save typing
	e_struct(end+1) = struct(...
		'element_full_name',  [e.element.name ' | ' int2str(e.element.reference) ] ,...
		'element_name',e.element.name,...
		'element_reference',e.element.reference,...
		'element_class',e.element.ndi_element_class,...
		'element_doc_id',e.base.id,...
		'element_type',e.element.type,...
		'element_session_id',e.base.session_id,...
		'element_subject_id',el{i}.dependency_value('subject_id'));
end;

s_struct = did.datastructures.emptystruct('subject_name','subject_id');
for i=1:numel(s),
	s_ = s{i}.document_properties;
	s_struct(end+1) = struct('subject_name',s_.subject.local_identifier,'subject_id',s_.base.id);
end;

e_table = struct2table(e_struct,'AsArray',true);
s_table = struct2table(s_struct,'AsArray',true);

T = join(e_table,s_table,'LeftKey',{'element_subject_id'},'RightKey',{'subject_id'});

