classdef testHartleyCalc < ndi.unittest.calc.sessionSetup
	methods (Test)
		function test_calculation(testCase)
			obj = ndi.calc.vis.hartley_calc(testCase.S);
			obj.test('standard',1,0);
		end
	end
end
