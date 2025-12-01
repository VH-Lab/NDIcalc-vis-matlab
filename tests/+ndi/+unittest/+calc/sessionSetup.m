classdef sessionSetup < matlab.unittest.TestCase
	properties
		S
		SessionPath
	end

	methods (TestMethodSetup)
		function setupSession(testCase)
			% SETUP_SESSION - create a temporary ndi.session.dir object for testing
			%
			% Creates an empty ndi.session.dir object at a temporary location.
			% The session object is stored in testCase.S and the path in testCase.SessionPath.

			ref = 'test_session';
			% Use tempname to generate a unique path.
			% tempname returns a path in the system temp directory.
			testCase.SessionPath = [tempname '_test_session'];

			if ~isfolder(testCase.SessionPath)
				mkdir(testCase.SessionPath);
			end

			% If we are creating an ndi.session.dir object for the first time,
			% we need to pass both the REFERENCE string (first input argument)
			% and the PATHNAME (second input argument).
			testCase.S = ndi.session.dir(ref, testCase.SessionPath);
		end
	end

	methods (TestMethodTeardown)
		function teardownSession(testCase)
			% TEARDOWN_SESSION - clean up the session directory
			%
			% Removes the temporary session directory created in setupSession.

			if isfolder(testCase.SessionPath)
				rmdir(testCase.SessionPath, 's');
			end
		end
	end
end
