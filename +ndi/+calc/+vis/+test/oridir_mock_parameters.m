function P = oridir_mock_parameters()
% ndi.calc.vis.test.oridir_mock_parameters - generate parameters for a mock document test
%
% P = ndi.calc.vis.test.oridir_mock_parameters()
%
% Generate a set of parameters P = [Rsp Rp Rn sig op] that describe a direction tuning curve.
%
% The parameters are variables from the double Gaussian function provided in the work by Mazurek et al. (2014). The function was used to extract 
% parameters of orientation and direction tuning curves with fits. Below is the list of parameters and how they are chosen:
% - Rsp: the untuned responseof the neuron, chosen as a real-value between 0 and 10. 
% - Rp: the above-offset response to the preferred orientation, chosen as a real-value between 0 and 30.
% - Rn: the above-offset response to the null direction, chosen as a real-value between 0 and Rp since its response would not be higher than the preferred response. 
% - Sig: the tuning width parameter, chosen as a real-value between 10 and 90. 
% - op: the preferred angle, chosen as a real-value between between 0 and 360. 
%
% See: Mazurek et al. (2014)

  % Step 1, generate random values

Rsp = 0 + (10-0)*rand;
Rp = 0 + (30-0)*rand;
Rn = 0 + (Rp-0)*rand;
sig = 10 + (90-10)*rand;
op = 0 + (360-0)*rand;

  % Step 2, assemble them as an array
  
P = [ Rsp Rp Rn sig op ];
