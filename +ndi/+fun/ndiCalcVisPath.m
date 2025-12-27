function p = ndiCalcVisPath()
% NDICALCVISPATH - return the path of the NDIcalc-vis-matlab root directory
%
% P = ndi.fun.ndiCalcVisPath()
%
% Returns the root directory of the NDIcalc-vis-matlab package.

p = fileparts(fileparts(fileparts(mfilename('fullpath'))));
