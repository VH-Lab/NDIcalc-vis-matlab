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
			% one serves. A round 20 is deliberate: it reads as a test fixture
			% rather than as a measurement from somebody's rig, which 21.0526
			% -- the value the workflow used -- did not.
			%
			% Fixing it also keeps the mock stimulus independent of whatever
			% calibration the caller has loaded. Without that, the stimulus,
			% and so the model responses and the spike-triggered average, would
			% differ from machine to machine, and the stored expectation would
			% only ever match on the machine that produced it.
			%
			% Nothing currently compared depends on this: the five fields
			% mock.1.compare.json checks -- M, L_max, K_max, sf_max, fps -- are
			% echoed from the input stimulus parameters, and everything derived
			% through pixels_per_degree sits under hartley_numbers, which is
			% 'none'. If those are ever compared, this constant becomes
			% load-bearing and the stored expectation must be regenerated at
			% whatever value is pinned here.
			%
			% The caller's own value is put back in TestClassTeardown, empty
			% included.

			NewStimGlobals;
			global pixels_per_cm
			testCase.SavedPixelsPerCm = pixels_per_cm;
			pixels_per_cm = 20;
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
