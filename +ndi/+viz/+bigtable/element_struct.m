function estruct = element_struct(S,element_doc)
% ELEMENT_STRUCT - make a structure based on an ndi element document
%
% ESTRUCT = ELEMENT_STRUCT(S, ELEMENT_DOC)
%
% Creates a structure with the element full name, name, reference,...
%  class, type, session ID, and document ID.
%

e = element_doc.document_properties; % save typing

estruct.element_full_name = [e.element.name ' | ' int2str(e.element.reference) ];
estruct.element_name = e.element.name;
estruct.element_reference = e.element.reference;
estruct.element_class = e.element.ndi_element_class;
estruct.element_type = e.element.type;
estruct.element_doc_id = e.base.id;
estruct.element_session_id = e.base.session_id;

q = ndi.query('base.id','exact_string',element_doc.dependency_value('subject_id'));

e_s = S.database_search(q);

if isempty(e_s),
	error(['Could not find subject for element ' estruct.element_full_name '.']);
end;

estruct.subject_name = e_s{1}.document_properties.subject.local_identifier;
