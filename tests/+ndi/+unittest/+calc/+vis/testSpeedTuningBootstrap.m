classdef testSpeedTuningBootstrap < ndi.unittest.calc.sessionSetup
	% Tests for ndi.calc.vis.speed_tuning_bootstrap
	%
	% The component functions (vis.speed.bootstrap.*) have their own unit
	% tests in tests/+vis/+speed/+bootstrap.

	methods (Test)

		function test_construction(testCase)
			% Constructing the calculator loads and validates the
			% speedtuning_bootstrap_calc / speed_tuning_bootstrap document types.
			obj = ndi.calc.vis.speed_tuning_bootstrap(testCase.S);
			testCase.verifyClass(obj, 'ndi.calc.vis.speed_tuning_bootstrap');
			testCase.verifyEqual(obj.numberOfSelfTests, 4);
		end

		function test_documentation_files(testCase)
			% The three documentation files exist and can be read.
			for t = {'general', 'searching for inputs', 'output'}
				txt = ndi.calculator.docfiletext('ndi.calc.vis.speed_tuning_bootstrap', t{1});
				testCase.verifyNotEmpty(txt, ['Empty documentation: ' t{1}]);
			end
		end

		function test_default_input_parameters(testCase)
			% Defaults are complete, and missing fields are filled in.
			ip = ndi.calc.vis.speed_tuning_bootstrap.default_input_parameters();
			testCase.verifyEqual(ip.min_xi, -4);
			testCase.verifyEqual(ip.max_xi, 4);
			testCase.verifyEqual(ip.method, 'percentile');
			ip2 = ndi.calc.vis.speed_tuning_bootstrap.complete_input_parameters( ...
				struct('numBootstrap', 50, 'deltaSpeedTuned', [0.25; 0.5]));
			testCase.verifyEqual(ip2.numBootstrap, 50);
			testCase.verifyEqual(ip2.deltaSpeedTuned, [0.25 0.5]);   % made a row
			testCase.verifyEqual(ip2.alpha, 0.05);                   % default filled in
			% a range that does not contain 0 and 1 in its interior warns
			testCase.verifyWarning(@() ndi.calc.vis.speed_tuning_bootstrap.complete_input_parameters( ...
				struct('min_xi', 0, 'max_xi', 1)), 'ndi:calc:vis:speed_tuning_bootstrap:xiRange');
		end

		function test_trials_matrix(testCase)
			ind = {[1 2 3], [4 5], [], 6};
			T = ndi.calc.vis.speed_tuning_bootstrap.trials_matrix(ind);
			testCase.verifySize(T, [4 3]);
			testCase.verifyEqual(T(1,:), [1 2 3]);
			testCase.verifyEqual(T(2,:), [4 5 NaN]);
			testCase.verifyTrue(all(isnan(T(3,:))));
			testCase.verifyEqual(T(4,1), 6);
		end

		function test_self_tests(testCase)
			% Field-by-field comparison with the stored expected documents.
			% Skipped until the expected documents mock.N.json have been
			% generated in a licensed-MATLAB session (see the README in
			% +ndi/+calc/+vis/mock/speed_tuning_bootstrap).
			obj = ndi.calc.vis.speed_tuning_bootstrap(testCase.S);
			for i = 1:obj.numberOfSelfTests
				testCase.assumeTrue(isfile(obj.mock_expected_filename(i)), ...
					['Expected document ' obj.mock_expected_filename(i) ...
					' has not been generated yet; see the mock README.']);
			end
			obj.verifySelfTests(testCase);
		end

		function test_run_and_read_fits(testCase)
			% Run the calculator on one mock cell and check the stored fits
			% and the recomputation of indexes from them.
			obj = ndi.calc.vis.speed_tuning_bootstrap(testCase.S);
			[param_struct, independent_variable, x, r] = obj.generate_mock_parameters('highSNR', 2);
			docs = ndi.mock.fun.stimulus_response(testCase.S, param_struct, ...
				independent_variable, x, r, 0.001, 5, 'RandomSeed', 'shuffle');
			parameters = obj.default_search_for_input_parameters();
			parameters.input_parameters.numBootstrap = 20;
			parameters.input_parameters.deltaSpeedTuned = [0.25 0.5];
			parameters.query.query = parameters.query.query & ...
				ndi.query('', 'depends_on', 'element_id', docs{3}.id());
			out = obj.run('Replace', parameters);
			testCase.verifyNumElements(out, 1);
			doc = out{1};
			sb = doc.document_properties.speed_tuning_bootstrap;

			B = obj.read_bootstrap_fits(doc);
			testCase.verifySize(B, [20 9]);
			testCase.verifyEqual(doc.document_properties.ngrid.data_dim(:).', [20 9]);
			% mock cell 2 is speed tuned (xi = 1)
			testCase.verifyEqual(sb.indexes.xi.estimate, 1, 'AbsTol', 0.1);
			testCase.verifyEqual(sb.indexes.inseparable.differentFrom0, 1);
			testCase.verifyNumElements(sb.indexes.speedTuned.equivalence, 2);   % one per margin
			% within 0.5 of 1 (a rule with a margin; "includes exactly 1" is
			% not checked, since at this noise level the interval is tiny)
			testCase.verifyEqual(sb.indexes.speedTuned.equivalence(2), 1);

			ix = obj.recompute_indexes(doc, 'deltaInseparable', [0.1 0.2 0.3]);
			testCase.verifyNumElements(ix.inseparable.minimumEffect, 3);
			testCase.verifyEqual(ix.xi.estimate, sb.indexes.xi.estimate, 'AbsTol', 1e-12);
		end

	end
end
