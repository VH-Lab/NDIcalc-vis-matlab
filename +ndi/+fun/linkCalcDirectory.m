function [linkPath, created] = linkCalcDirectory()
% LINKCALCDIRECTORY - make this package discoverable by NDI's calculator search
%
%  [LINKPATH, CREATED] = ndi.fun.linkCalcDirectory()
%
%  NDI-matlab discovers calculator schemas and database documents by scanning
%  the directory that contains NDI-matlab for sibling folders whose names match
%  NDIcalc*-matlab. A checkout of this repository that lives somewhere else, or
%  that was checked out under any other folder name -- the name a continuous
%  integration job chooses, or 'NDIcalc-vis-matlab-main' from a downloaded zip
%  archive -- is therefore not found, and the self-tests fail while reading the
%  JSON schema that defines each calculator's output type:
%
%     DID:Document:readjsonfilelocation: found no match for {calculator}
%
%  This function makes the current checkout discoverable, by creating a link to
%  it named NDIcalc-vis-matlab beside NDI-matlab (a symbolic link on Linux and
%  macOS, a directory junction on Windows; neither needs special privileges).
%  It is safe to call more than once: if the checkout can already be found,
%  nothing is created.
%
%  Returns LINKPATH, the location where NDI will look for this package, and
%  CREATED, 1 if a link was made on this call and 0 if none was needed.
%
%  Example:
%     ndi.fun.linkCalcDirectory();
%     ndi.calc.vis.oridir_tuning(S).test('highSNR',1,0);
%
%  See also: ndi.fun.ndiCalcVisPath
%

created = 0;

p = ndi.fun.ndiCalcVisPath();

if isempty(which('ndi.toolboxdir')),
	error('NDIcalc:linkCalcDirectory:noNDI', ...
		['NDI-matlab was not found on the MATLAB path. Add NDI-matlab to the ' ...
		 'path before calling ndi.fun.linkCalcDirectory.']);
end;

ndiParent = fileparts(fileparts(fileparts(ndi.toolboxdir)));
linkPath = fullfile(ndiParent,'NDIcalc-vis-matlab');

 % Case 1: this checkout already sits where NDI looks, under a matching name.

[thisParent,thisName] = fileparts(p);

if strcmp(thisParent,ndiParent) & ~isempty(regexp(thisName,'^NDIcalc.*-matlab$','once')),
	linkPath = p;
	return;
end;

 % Case 2: a folder of that name is already there. Leave it alone; it is either
 % this checkout reached by an existing link, or another copy that the user put
 % there deliberately. Either way NDI has something to find.

if isfolder(linkPath),
	return;
end;

 % Case 3: create the link.

if ispc,
	cmd = sprintf('mklink /J "%s" "%s"', linkPath, p);
else,
	cmd = sprintf('ln -s "%s" "%s"', p, linkPath);
end;

[status,msg] = system(cmd);

if status~=0,
	error('NDIcalc:linkCalcDirectory:linkFailed', ...
		['Could not create ''%s'' (%s).\n' ...
		 'NDI-matlab finds calculator packages by looking beside NDI-matlab ' ...
		 'for folders named NDIcalc*-matlab. Either move this checkout to ' ...
		 '''%s'', or create that link by hand, before running the self-tests.'], ...
		linkPath, strtrim(msg), linkPath);
end;

created = 1;
