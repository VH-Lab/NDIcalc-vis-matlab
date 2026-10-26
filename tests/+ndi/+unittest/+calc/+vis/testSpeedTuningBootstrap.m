classdef testSpeedTuningBootstrap < ndi.unittest.calc.sessionSetup
	methods (Test)
		function test_construction(testCase)
			% Smoke test: constructing the calculator loads and validates the
			% new speedtuning_bootstrap_calc / speed_tuning_bootstrap document
			% types on the calculator document path.
			%
			% The field-by-field self-tests (generate_mock_parameters + stored
			% expected documents) are a follow-up: the expected mock documents
			% must be generated in a licensed-MATLAB session, as for every other
			% calculator here, then numberOfSelfTests raised above 0 and a
			% verifySelfTests call added here. verifySelfTests requires
			% numberOfSelfTests to be positive, so it is not called while the
			% count is 0.
			obj = ndi.calc.vis.speed_tuning_bootstrap(testCase.S);
			testCase.verifyClass(obj, 'ndi.calc.vis.speed_tuning_bootstrap');
			testCase.verifyEqual(obj.numberOfSelfTests, 0);
		end
	end
end
