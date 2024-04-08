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

q = ndi.query('','isa','nothing');

if isempty(stimulus_purpose),
	return;
end;

if ~iscell(stimulus_purpose),
	stimulus_purpose = {stimulus_purpose};
end;

q_approach = ndi.query('openminds.openminds_type','exact_string','https://openminds.ebrains.eu/controlledTerms/StimulationApproach') | ...
	ndi.query('openminds.openminds_type','exact_string','https://openminds.ebrains.eu/controlledTerms/StimulationPurpose');

i = 1;
q_purpose = ndi.query('openminds.fields.name','exact_string',stimulus_purpose{i});

for i=2:numel(stimulus_purpose),
	q_purpose = q_purpose | ndi.query('openminds.fields.name','exact_string',stimulus_purpose{i});
end;

doc = S.database_search(q_purpose & q_approach);

 % now we have all stimulus purposes individually, but have not ANDED them together

epochid = {};
for i=1:numel(stimulus_purpose),
	epochid{i} = {};
	for j=1:numel(doc),
		if strcmp(doc{j}.document_properties.openminds.fields.name,stimulus_purpose{i}),
			epochid{i}{end+1} = doc{j}.document_properties.epochid.epochid;
		end;
	end;
end;

unique_epochs = unique(epochid{1});
for i=2:numel(stimulus_purpose),
	unique_epochs = intersect(unique_epochs,epochid{i});
end;

if ~isempty(unique_epochs),
	i = 1;
	q = ndi.query('stimulus_response.stimulator_epochid','exact_string',unique_epochs{i});
	for i=2:numel(unique_epochs),
		q = q|ndi.query('stimulus_response.stimulator_epochid','exact_string',unique_epochs{i});
	end;
else,
	q = ndi.query('','isa','nothing');
end;


