classdef test_speed_nested_f < matlab.unittest.TestCase

    methods (Test)

        function testSpeedNestedF(testCase)
            % Test vis.speed.speed_nested_f F-test calculation

            % Case 1: SSE_withspeed << SSE_nospeed (Significant improvement)
            response_vector_length = 100;
            sse_with = 10;
            sse_no = 100; % Huge error without speed param

            p = vis.speed.speed_nested_f(response_vector_length, sse_with, sse_no);

            % Should be a very small p-value (significant)
            testCase.verifyLessThan(p, 0.05, 'Should be significant when speed model fits much better.');

            % Case 2: SSE_withspeed ~ SSE_nospeed (Not significant)
            sse_with_2 = 10;
            sse_no_2 = 10.1; % Tiny improvement

            p2 = vis.speed.speed_nested_f(response_vector_length, sse_with_2, sse_no_2);

            % Should be large p-value (not significant)
            testCase.verifyGreaterThan(p2, 0.05, 'Should not be significant when models fit similarly.');
        end

    end
end
