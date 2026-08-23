classdef test_generate_STA < matlab.unittest.TestCase
    % TEST_GENERATE_STA - unit tests for vis.revcorr.generate_STA
    %
    % These tests pin down the time convention of the reverse correlation:
    % which reconstruction time each plane of the returned array corresponds
    % to. That convention is easy to get wrong and easy to miss, because an
    % error in it leaves the images intact and only relabels the time axis.
    %
    % The convention, stated once here and asserted below:
    %
    %   T_coords is a set of lags, measured backwards from each event. A lag
    %   of +0.100 s means the stimulus that was on the screen 0.100 s BEFORE
    %   the event; a lag of -0.100 s means the stimulus 0.100 s after it.
    %
    %   vis.revcorr.generate_STA returns the reconstruction in reverse-lag
    %   order, so that plane s holds the stimulus at
    %
    %       t_event - T_coords(TMAX+1-s).
    %
    %   ndi.calc.vis.hartley reverses the third dimension of this array
    %   before storing it, after which plane k holds the stimulus at
    %   t_event - T_coords(k), matching T_coords element for element. The
    %   tests below therefore reverse the array the same way and compare
    %   against T_coords directly, which is the form a reader of the stored
    %   document sees.
    %
    % A test that used a T_coords symmetric about zero would pass whether or
    % not the convention were implemented correctly, since reversing the
    % order of a symmetric set of lags is the same as negating it. The lag
    % ranges used here are deliberately asymmetric, as the calculator's own
    % defaults are.

    methods (Test)

        function testLagConvention(testCase)
            % Each plane must hold the stimulus that was on the screen at the
            % lag that labels it. With a single event the average is just the
            % stimulus movie, so every plane can be checked exactly.

            M = 4;
            frameDuration = 0.05;
            [s, kx, ky, frameTimes] = testCase.makeStimulus(41, frameDuration);

            % Lags: asymmetric about zero, as the calculator's default is.
            T_coords = (-0.10:0.05:0.25)';
            tmax = round((T_coords(end)-T_coords(1))/mean(diff(T_coords))) + 1;

            % Offset from the frame boundaries so that no sample falls on the
            % instant a frame changes, where either frame would be defensible.
            tEvent = 0.62;

            expectedTimes = tEvent - T_coords;
            expected = testCase.imagesAt(s, kx, ky, frameTimes, M, expectedTimes);
            testCase.assertPairwiseDistinct(expected);

            sta = vis.revcorr.generate_STA(s, kx, ky, frameTimes, tEvent, ...
                T_coords(end)-T_coords(1), T_coords, tmax, M, @(fraction) []);

            % Reverse as ndi.calc.vis.hartley does, so plane k means lag T(k).
            sta = sta(:,:,end:-1:1);

            for k = 1:tmax
                testCase.verifyEqual(sta(:,:,k), expected(:,:,k), 'AbsTol', 1e-10, ...
                    sprintf(['Plane %d should hold the stimulus at lag %.3f s, ' ...
                    'which is the frame shown at t = %.3f s.'], ...
                    k, T_coords(k), expectedTimes(k)));
            end
        end

        function testCausalLagsExcludeFutureStimulus(testCase)
            % A lag range that is entirely non-negative asks only for stimuli
            % at or before the event. Nothing shown after the event may
            % appear in the reconstruction. This is the same defect as
            % testLagConvention seen from a direction that needs no
            % arithmetic: if the window is taken forwards in time instead of
            % backwards, every plane is drawn from the future stimulus.

            M = 4;
            frameDuration = 0.05;
            nFrames = 41;
            frameTimes = (0:nFrames-1)' * frameDuration;

            tEvent = 1.02;

            % One pattern before the event, a different one after it. The
            % change happens strictly after the event, so the frame on the
            % screen at the moment of the event is still the "past" pattern.
            past = frameTimes <= tEvent;
            s  = ones(nFrames,1);
            kx = zeros(nFrames,1);
            ky = zeros(nFrames,1);
            kx(past)  = 1;   ky(past)  = 0;   s(past)  = 1;
            kx(~past) = 0;   ky(~past) = 2;   s(~past) = -1;

            pastImage   = vlt.neuro.reverse_correlation.hartley.hartley_image(1, 1, 0, M);
            futureImage = vlt.neuro.reverse_correlation.hartley.hartley_image(-1, 0, 2, M);
            testCase.assertNotEqual(pastImage, futureImage, ...
                'Test setup error: the two patterns must differ.');

            T_coords = (0:0.05:0.20)';
            tmax = round((T_coords(end)-T_coords(1))/mean(diff(T_coords))) + 1;

            sta = vis.revcorr.generate_STA(s, kx, ky, frameTimes, tEvent, ...
                T_coords(end)-T_coords(1), T_coords, tmax, M, @(fraction) []);

            for k = 1:tmax
                testCase.verifyEqual(sta(:,:,k), pastImage, 'AbsTol', 1e-10, ...
                    sprintf(['Plane %d was reconstructed from the stimulus shown ' ...
                    'after the event; non-negative lags must use only the ' ...
                    'stimulus at or before it.'], k));
            end
        end

        function testAverageOverEvents(testCase)
            % The reconstruction is a mean over events: computing it for two
            % events must give the mean of the two single-event
            % reconstructions. This checks the 1/N normalization and is
            % independent of the lag convention.

            M = 4;
            frameDuration = 0.05;
            [s, kx, ky, frameTimes] = testCase.makeStimulus(41, frameDuration);

            T_coords = (-0.10:0.05:0.25)';
            tmax = round((T_coords(end)-T_coords(1))/mean(diff(T_coords))) + 1;
            rfRange = T_coords(end)-T_coords(1);

            tEvents = [0.62; 1.13];

            both = vis.revcorr.generate_STA(s, kx, ky, frameTimes, tEvents, ...
                rfRange, T_coords, tmax, M, @(fraction) []);
            first = vis.revcorr.generate_STA(s, kx, ky, frameTimes, tEvents(1), ...
                rfRange, T_coords, tmax, M, @(fraction) []);
            second = vis.revcorr.generate_STA(s, kx, ky, frameTimes, tEvents(2), ...
                rfRange, T_coords, tmax, M, @(fraction) []);

            testCase.verifyEqual(both, (first+second)/2, 'AbsTol', 1e-10, ...
                'The reconstruction must be the mean over events.');
        end

        function testOutputSize(testCase)
            % The reconstruction is M by M by TMAX.

            M = 6;
            [s, kx, ky, frameTimes] = testCase.makeStimulus(41, 0.05);

            T_coords = (-0.10:0.05:0.25)';
            tmax = round((T_coords(end)-T_coords(1))/mean(diff(T_coords))) + 1;

            sta = vis.revcorr.generate_STA(s, kx, ky, frameTimes, 0.62, ...
                T_coords(end)-T_coords(1), T_coords, tmax, M, @(fraction) []);

            testCase.verifySize(sta, [M M tmax]);
        end

    end

    methods (Access = private)

        function [s, kx, ky, frameTimes] = makeStimulus(~, nFrames, frameDuration)
            % A Hartley sequence in which no two frames within any window
            % shorter than 32 frames show the same image, so that a plane of
            % the reconstruction identifies the frame it came from.

            i = (0:nFrames-1)';
            frameTimes = i * frameDuration;
            kx = mod(i, 4);
            ky = mod(floor(i/4), 4);
            s  = 1 - 2*mod(floor(i/16), 2);
        end

        function IM = imagesAt(~, s, kx, ky, frameTimes, M, times)
            % The images on the screen at each of TIMES, as a stack.

            IM = zeros(M, M, numel(times));
            for j = 1:numel(times)
                i = find(frameTimes <= times(j), 1, 'last');
                IM(:,:,j) = vlt.neuro.reverse_correlation.hartley.hartley_image(...
                    s(i), kx(i), ky(i), M);
            end
        end

        function assertPairwiseDistinct(testCase, stack)
            % Guard the test itself: comparing planes against expected images
            % proves nothing if two of those images are identical.

            n = size(stack,3);
            for a = 1:n
                for b = a+1:n
                    testCase.assertGreaterThan(max(max(abs(stack(:,:,a)-stack(:,:,b)))), 1e-6, ...
                        sprintf(['Test setup error: expected images %d and %d are ' ...
                        'identical, so the test could not tell them apart.'], a, b));
                end
            end
        end

    end

end
