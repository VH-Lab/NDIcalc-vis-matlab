function b=test_all_stims_present(d, dold)
% TEST_ALL_STIMS_PRESENT - test whether all expected stimuli are present in a newly made document
%
% B = TEST_ALL_STIMS_PRESENT(D, DOLD)
%
% Given a new document D and the old document DOLD, see if all stimuli are represented exactly once.
%
% This function assumes that the stimuli vary in fields sFrequency, angle, and contrast.
%
% It checks to make sure that each stimulus combination is represented exactly once, and with
% angles equal to the angles in dold and 180 added to all angles in dold.
%

 % only checking contrast, angle, sFrequency

SFs = [];
angles = [];
contrasts = [];

for i=1:numel(dold.document_properties.stimulus_presentation.stimuli),
	if ~isfield(dold.document_properties.stimulus_presentation.stimuli(i).parameters,'isblank'),
		SFs(end+1) = dold.document_properties.stimulus_presentation.stimuli(i).parameters.sFrequency;
		angles(end+1) = dold.document_properties.stimulus_presentation.stimuli(i).parameters.angle;
		contrasts(end+1) = dold.document_properties.stimulus_presentation.stimuli(i).parameters.contrast;
	end;
end;

SFs = unique(SFs);
angles = unique(angles);
contrasts = unique(contrasts);

N = numel(angles);
for i=1:N,
	angles(end+1) = angles(i) + 180;
end;

A = zeros(numel(SFs),numel(angles),numel(contrasts));

for i=1:numel(d.document_properties.stimulus_presentation.stimuli),
	if ~isfield(d.document_properties.stimulus_presentation.stimuli(i).parameters,'isblank'),
		sf_here = d.document_properties.stimulus_presentation.stimuli(i).parameters.sFrequency;
		angle_here = d.document_properties.stimulus_presentation.stimuli(i).parameters.angle;
		contrast_here = d.document_properties.stimulus_presentation.stimuli(i).parameters.contrast;
	
		sf_i = find(sf_here==SFs);
		angle_i = find(angle_here==angles);
		contrast_i = find(contrast_here==contrasts);

		if isempty(sf_i),
			error(['SF not found.']);
		end;

		if isempty(angle_i),
			error(['Angle not found.']);
		end;

		if isempty(contrast_i),
			error('contrast not found.');
		end;
		
		A(sf_i,angle_i,contrast_i) = A(sf_i,angle_i,contrast_i) + 1;
	
	end;
end;

b = vlt.data.eqlen(A,ones(numel(SFs),numel(angles),numel(contrasts)));

