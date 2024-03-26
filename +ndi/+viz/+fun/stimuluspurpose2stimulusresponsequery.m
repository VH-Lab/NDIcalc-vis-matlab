function q = stimuluspurpose2stimulusresponsequery(S, stimulus_purpose)
% STIMULUSPURPOSE2STIMULUSRESPONSEQUERY - get query for stimulus response documents based on purpose
%
% Q = STIMULUSPURPOSE2STIMULUSRESPONSEQUERY(S, STIMULUS_PURPOSE)
%
% Return an ndi.query object that will find stimulus_response documents for the 
% stimulator epoch with purpose name STIMULUS_PURPOSE.
% If no stimulus purpose match is found, then q is ndi.query('','isa','nothing') .
%
% Example:
%    q = ndi.viz.fun.stimuluspurpose2stimulusresponsequery(S, 'Purpose: Assessing contrast tuning');
%

doc = S.database_search(ndi.query('openminds.fields.name','exact_string',purpose_string) & ...
	ndi.query('openminds.openminds_type','exact_string','https://openminds.ebrains.eu/controlledTerms/StimulationApproach'));

if ~isempty(doc),
	q = ndi.query('stimulus_response.stimulator_epochid','exact_string',doc{1}.document_properties.epochid.epochid);
	for i=2:numel(doc),
		q = q|ndi.query('stimulus_response.stimulator_epochid','exact_string',doc{i}.document_properties.epochid.epochid);
	end;
else,
	q = ndi.query('','isa','nothing');
end;

