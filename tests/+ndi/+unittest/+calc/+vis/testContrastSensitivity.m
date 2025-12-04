classdef testContrastSensitivity < ndi.unittest.calc.sessionSetup
	methods (Test)
		function test_calculation(testCase)
			obj = ndi.calc.vis.contrast_sensitivity(testCase.S);
			obj.test('highSNR',1,0);
		end
	end
end
