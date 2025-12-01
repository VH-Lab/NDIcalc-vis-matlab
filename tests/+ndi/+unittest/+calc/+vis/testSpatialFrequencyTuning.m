classdef testSpatialFrequencyTuning < ndi.unittest.calc.sessionSetup
	methods (Test)
		function test_calculation(testCase)
			obj = ndi.calc.vis.spatial_frequency_tuning(testCase.S);
			obj.test('standard',1,0);
		end
	end
end
