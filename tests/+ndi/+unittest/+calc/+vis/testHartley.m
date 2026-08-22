classdef testHartley < ndi.unittest.calc.sessionSetup
	properties (Access = private)
		SavedPixelsPerCm
	end

	methods (TestClassSetup)
		function setupNewStimGlobals(testCase)
			% SETUPNEWSTIMGLOBALS - give the mock Hartley stimulus a display scale
			%
			% hartleystim converts spatial frequency to cycles per degree of
			% visual angle using NewStim's global pixels_per_cm. Unset, that
			% conversion factor is empty, and hartleyrange dies at
			%
			%    F_ = F*pixels_per_degree;
			%
			% with MATLAB:innerdim -- before the calculator runs, and naming
			% nothing the reader can act on. Continuous integration set this up
			% in the workflow, so the test passed there and failed for anyone
			% running it from MATLAB. It belongs with the test.
			%
			% The value is a display scale for synthetic data, so any plausible
			% one will do; it is fixed here so the mock stimulus is identical
			% on every machine. The caller's own value is put back afterwards.

			NewStimGlobals;
			global pixels_per_cm
			testCase.SavedPixelsPerCm = pixels_per_cm;
			pixels_per_cm = 21.0526;
		end
	end

	methods (TestClassTeardown)
		function restoreNewStimGlobals(testCase)
			global pixels_per_cm
			pixels_per_cm = testCase.SavedPixelsPerCm;
		end
	end

	methods (Test)
		function test_calculation(testCase)
			obj = ndi.calc.vis.hartley(testCase.S);
			% Every self-test the calculator declares, not just the first.
			% Omitting testIndexes runs 1:numberOfSelfTests, so this stays
			% correct if that number changes.
			obj.verifySelfTests(testCase);
		end
	end
end
