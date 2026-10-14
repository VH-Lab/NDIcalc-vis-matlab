function D = openDataset(targetPath, options)
% OPENDATASET - Open the calculators paper example dataset, downloading it if needed
%
% D = NDI.CALC.PAPERFIGURES.OPENDATASET()
% D = NDI.CALC.PAPERFIGURES.OPENDATASET(TARGETPATH)
%
% Returns an ndi.dataset object for the example dataset that accompanies the
% calculators paper ("A motif for domain-specific analysis applets that are
% easy to learn, reuse, test, and to compose into pipelines: application to
% vision science"), which contains two example sessions.
%
% The dataset lives in the folder [TARGETPATH filesep CLOUDDATASETID]. If that
% folder already holds an NDI dataset, it is opened with ndi.dataset.dir and
% nothing is downloaded. Otherwise the dataset is downloaded from NDI Cloud
% with ndi.cloud.downloadDataset, with file syncing off, so only the documents
% are downloaded and binary files are fetched on demand.
%
% TARGETPATH is the parent folder; it defaults to MATLAB's USERPATH.
%
% This function also takes name/value pairs that modify its behavior:
% |----------------------------|----------------------------------------------|
% | Parameter (default)        | Description                                  |
% |----------------------------|----------------------------------------------|
% | cloudDatasetId             | The NDI Cloud identifier of the dataset.     |
% |   (the paper's dataset)    |                                              |
% | Verbose (true)             | Print progress while downloading.            |
% |----------------------------|----------------------------------------------|
%
% Example:
%   D = ndi.calc.paperFigures.openDataset();
%
% See also: ndi.cloud.downloadDataset, ndi.dataset.dir

arguments
    targetPath (1,1) string = string(userpath)
    options.cloudDatasetId (1,1) string = paperCloudDatasetId()
    options.Verbose (1,1) logical = true
end

if strlength(targetPath) == 0
    error('NDICALC:paperFigures:openDataset:NoTargetPath', ...
        ['USERPATH is empty, so there is no default location for the dataset. ' ...
        'Pass a folder as TARGETPATH.']);
end

if ~isfolder(targetPath)
    mkdir(targetPath);
end

datasetFolder = fullfile(targetPath, options.cloudDatasetId);

if isfolder(fullfile(datasetFolder, '.ndi'))
    if options.Verbose
        disp(['Opening existing dataset at ' char(datasetFolder) '.']);
    end
    D = ndi.dataset.dir(char(datasetFolder));
else
    D = ndi.cloud.downloadDataset(options.cloudDatasetId, targetPath, ...
        'SyncFiles', false, 'Verbose', options.Verbose);
end

end

function id = paperCloudDatasetId()
% PAPERCLOUDDATASETID - The NDI Cloud identifier of the paper's example dataset
id = ""; % TODO: fill in the NDI Cloud dataset ID of the paper's example dataset
if strlength(id) == 0
    error('NDICALC:paperFigures:openDataset:NoDatasetId', ...
        ['The NDI Cloud ID of the paper''s example dataset has not been filled ' ...
        'in yet. Pass it as ''cloudDatasetId''.']);
end
end
