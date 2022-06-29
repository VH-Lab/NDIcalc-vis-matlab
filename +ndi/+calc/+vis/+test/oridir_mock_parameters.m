function P = oridir_mock_parameters()
% ndi.calc.vis.test.oridir_mock_parameters - generate parameters for a mock document test
%
% P = ndi.calc.vis.test.oridir_mock_parameters()
%
% Generate a set of parameters P = [Rsp Rp Rn sig op] that describe a
% direction tuning curve.
%
% %(Ray, describe how the parameters are chosen here
%
% See: Mazurek et al. (2014)

  % Step 1, generate random values

Rsp = randi([0 10]) % Ray, this will only choose integer values, use 0 + (10-0)*rand to generate a real-valued number between 0 and 10
Rp = randi([0 30]) % change to real-valued
Rn = randi([0 Rp]) % change to real-valued
sig = randi([10 90]) % change to real-valued
op = randi([0 360]) % change to real-valued

  % Step 2, assemble them as an array
  
P = [ Rsp Rp Rn sig op ];
