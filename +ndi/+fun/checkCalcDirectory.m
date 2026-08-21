function [tf, expectedPath, msg] = checkCalcDirectory()
% CHECKCALCDIRECTORY - check that NDI can discover this calculator package
%
%  ndi.fun.checkCalcDirectory()
%  [TF, EXPECTEDPATH, MSG] = ndi.fun.checkCalcDirectory()
%
%  NDI-matlab discovers calculator schemas and database documents by scanning
%  the directory that contains NDI-matlab for sibling folders whose names match
%  NDIcalc*-matlab. A checkout of this repository that lives somewhere else, or
%  under any other folder name, is not found. The calculators still load, but
%  their self-tests fail while reading the JSON schema that defines each one's
%  output type, with an error that does not name the cause:
%
%     DID:Document:readjsonfilelocation: found no match for {calculator}
%
%  Folder names that do not match are easier to end up with than renaming
%  suggests: GitHub's "Download ZIP" expands to NDIcalc-vis-matlab-main, an
%  archived release tarball expands to NDIcalc-vis-matlab-1.0.0, and a
%  continuous integration job checks out under whatever name it chooses.
%
%  Called with no outputs, this function raises an error naming the checkout,
%  the directory NDI searches, and the exact command to link the two, or does
%  nothing at all if the checkout can already be found.
%
%  Called with outputs it never raises an error: TF is 1 if the package can be
%  discovered and 0 if it cannot, EXPECTEDPATH is where NDI looks for it, and
%  MSG is the explanation, empty when TF is 1.
%
%  This function does not create anything. To make a non-matching checkout
%  discoverable, follow the command it reports.
%
%  Every ndi.calc.vis calculator constructor calls this before its superclass
%  constructor, and the self-test harness calls it before the first test, so
%  an undiscoverable copy is reported at the point of use rather than as a
%  document-type error further down. It is cheap: two path splits and a
%  directory test.
%
%  Example:
%     ndi.fun.checkCalcDirectory();
%     ndi.calc.vis.oridir_tuning(S).test('highSNR',1,0);
%
%  See also: ndi.fun.ndiCalcVisPath
%

p = ndi.fun.ndiCalcVisPath();

if isempty(which('ndi.toolboxdir')),
	error('NDIcalc:checkCalcDirectory:noNDI', ...
		['NDI-matlab was not found on the MATLAB path. Add NDI-matlab to the ' ...
		 'path before running the calculators or their self-tests.']);
end;

ndiParent = fileparts(fileparts(fileparts(ndi.toolboxdir)));
expectedPath = fullfile(ndiParent,'NDIcalc-vis-matlab');

tf = 0;
msg = '';

 % Discoverable if this checkout sits beside NDI-matlab under a matching name,
 % or if some folder with a matching name is already there (this checkout
 % reached by a link, or another copy the user placed deliberately).

[thisParent,thisName] = fileparts(p);

if strcmp(thisParent,ndiParent) & ~isempty(regexp(thisName,'^NDIcalc.*-matlab$','once')),
	tf = 1;
	expectedPath = p;
elseif isfolder(expectedPath),
	tf = 1;
end;

if tf,
	return;
end;

if ispc,
	linkcmd = sprintf('mklink /J "%s" "%s"', expectedPath, p);
else,
	linkcmd = sprintf('ln -s "%s" "%s"', p, expectedPath);
end;

msg = sprintf([ ...
	'NDI cannot discover this copy of NDIcalc-vis-matlab.\n\n' ...
	'  This checkout:   %s\n' ...
	'  NDI searches:    %s\n\n' ...
	'NDI-matlab looks beside itself for folders named NDIcalc*-matlab, and ' ...
	'nothing there matches. The calculators will load, but their self-tests ' ...
	'will fail while reading each calculator''s output type, reporting\n' ...
	'  DID:Document:readjsonfilelocation: found no match for {calculator}\n\n' ...
	'To fix this, either move this checkout to\n  %s\nor link it there:\n  %s\n'], ...
	p, ndiParent, expectedPath, linkcmd);

if nargout==0,
	error('NDIcalc:checkCalcDirectory:notDiscoverable','%s',msg);
end;
