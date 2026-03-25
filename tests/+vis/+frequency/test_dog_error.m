classdef test_dog_error < matlab.unittest.TestCase

    methods (Test)

        function testZeroError(testCase)
            % When data matches the model exactly, error should be zero
            params = [0 10 1 5 2]; % R0=0, RE=10, SE=1, RI=5, SI=2
            x = linspace(-3, 3, 20);
            % Generate y from the model itself
            par = abs(params);
            y = par(1) + par(2).*exp(-x.^2/2./par(3)^2) ...
                - par(4).*exp(-x.^2/2./par(5)^2);
            err = vis.frequency.dog_error(params, x, y);
            testCase.verifyLessThan(abs(err), 1e-10, ...
                'Error should be near zero when data matches model');
        end

        function testNonZeroError(testCase)
            % When data does not match model, error should be positive
            params = [0 10 1 5 2];
            x = linspace(-3, 3, 20);
            y = ones(size(x)); % constant data, won't match DoG
            err = vis.frequency.dog_error(params, x, y);
            testCase.verifyGreaterThan(err, 0, ...
                'Error should be positive for mismatched data');
        end

        function testWithSte(testCase)
            % Test that the 4-argument form (with STE) runs without error
            params = [1 10 1 5 2];
            x = linspace(-3, 3, 10);
            par = abs(params);
            y = par(1) + par(2).*exp(-x.^2/2./par(3)^2) ...
                - par(4).*exp(-x.^2/2./par(5)^2);
            ste = 0.1 * ones(size(x));
            err = vis.frequency.dog_error(params, x, y, ste);
            testCase.verifyLessThan(abs(err), 1e-6, ...
                'Error with STE should be near zero for matching data');
        end

    end
end
