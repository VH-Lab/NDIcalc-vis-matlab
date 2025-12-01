classdef testContrastSensitivity < ndi.unittest.calc.sessionSetup
	methods (Test)
		function test_calculation(testCase)
			obj = ndi.calc.vis.contrast_sensitivity(testCase.S);
			obj.test();
		end
	end
end
