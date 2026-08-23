classdef testOridirTuning < ndi.unittest.calc.sessionSetup
	methods (Test)
		function test_calculation(testCase)
			obj = ndi.calc.vis.oridir_tuning(testCase.S);
			% Every self-test the calculator declares, not just the first.
			% Omitting testIndexes runs 1:numberOfSelfTests, so this stays
			% correct if that number changes.
			obj.verifySelfTests(testCase);
		end
	end
end
