function d = stimuluspurpose2stimulusresponse(S, stimulus_purpose)
% STIMULUSPURPOSE2STIMULUSRESPONSE - get stimulus response documents based on purpose
%
% D = STIMULUSPURPOSE2STIMULUSRESPONS(S, STIMULUS_PURPOSE)
%
% Return ndi.document objects of type stimulus_response for the
% stimulator epoch with purpose name STIMULUS_PURPOSE.
% If no stimulus purpose match is found, then docs is {}.
%
% See also: stimuluspurpose2stimulusresponsequery
%
% Example:
%    docs = ndi.viz.fun.stimuluspurpose2stimulusresponse(S, 'Purpose: Assessing contrast tuning');
%

q = ndi.viz.fun.stimuluspurpose2stimulusresponse(S, stimulus_purpose);

d = S.database_search(q);


