classdef testSpeedTuning < ndi.unittest.calc.sessionSetup
	methods (Test)
		function test_calculation(testCase)
			obj = ndi.calc.vis.speed_tuning(testCase.S);
			obj.test();
		end
	end
end
