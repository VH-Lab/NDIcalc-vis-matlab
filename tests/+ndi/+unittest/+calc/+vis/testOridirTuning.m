classdef testOridirTuning < ndi.unittest.calc.sessionSetup
	methods (Test)
		function test_calculation(testCase)
			obj = ndi.calc.vis.oridir_tuning(testCase.S);
			% verifySelfTests qualifies the comparison result. The previous
			% obj.test('highSNR',1,0) computed the same result and discarded it,
			% so the test passed whatever answer the calculator produced.
			% testIndexes is 1 to run exactly the self-test that was run before.
			obj.verifySelfTests(testCase,'testIndexes',1);
		end
	end
end
