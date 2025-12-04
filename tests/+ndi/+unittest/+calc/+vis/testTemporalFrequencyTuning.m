classdef testTemporalFrequencyTuning < ndi.unittest.calc.sessionSetup
	methods (Test)
		function test_calculation(testCase)
			obj = ndi.calc.vis.temporal_frequency_tuning(testCase.S);
			obj.test('highSNR',1,0);
		end
	end
end
