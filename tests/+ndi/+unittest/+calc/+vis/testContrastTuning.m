classdef testContrastTuning < ndi.unittest.calc.sessionSetup
	methods (Test)
		function test_calculation(testCase)
			obj = ndi.calc.vis.contrast_tuning(testCase.S);
			obj.test();
		end
	end
end
