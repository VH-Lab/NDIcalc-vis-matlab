classdef testSpikeShape < ndi.unittest.calc.sessionSetup
	methods (Test)
		function test_calculation(testCase)
			obj = ndi.calc.vis.spike_shape(testCase.S);
			obj.test('standard',1,0);
		end
	end
end
