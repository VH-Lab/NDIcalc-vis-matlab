classdef testSelectRowsByModulationIndex < matlab.unittest.TestCase
% TESTSELECTROWSBYMODULATIONINDEX - Unit tests for ndi.viz.selectRowsByModulationIndex
%
% Verifies the modulation-index-based row selection and the optional F0/F1
% row-set outputs.

    methods (Static)
        function T = makeTable()
            % Build a small synthetic table with three cells:
            %   id 1/2: simple cell  (MI = 1.5, paired mean+F1)
            %   id 3/4: complex cell (MI = 0.5, paired mean+F1)
            %   id 5  : unpaired mean row (MI = NaN)
            id = (1:5)';
            mi = [1.5; 1.5; 0.5; 0.5; NaN];
            rt = {'mean'; 'F1'; 'mean'; 'F1'; 'mean'};
            T = table(id, mi, rt, 'VariableNames', ...
                {'id', 'x.TC.modulationIndex', 'x.properties.response_type'});
        end
    end

    methods (Test)

        function testSelectedRows(testCase)
            % T_out keeps F1 for simple cells (MI>=1) and mean for complex cells (MI<1),
            % and excludes unpaired (NaN-MI) rows.
            T = ndi.unittest.viz.testSelectRowsByModulationIndex.makeTable();
            T_out = ndi.viz.selectRowsByModulationIndex(T);

            testCase.verifyEqual(height(T_out), 2, ...
                'Expected one selected row per validly-paired cell.');
            % id 2 is the simple cell's F1 row; id 3 is the complex cell's mean row.
            testCase.verifyEqual(sort(T_out.id(:)), [2; 3]);
        end

        function testF0AndF1RowSets(testCase)
            % T_F0 and T_F1 contain one row per validly-paired cell, regardless of MI.
            T = ndi.unittest.viz.testSelectRowsByModulationIndex.makeTable();
            [T_out, T_F0, T_F1] = ndi.viz.selectRowsByModulationIndex(T);

            % T_F0: the 'mean' rows of the two paired cells (ids 1 and 3), not id 5 (NaN MI).
            testCase.verifyEqual(sort(T_F0.id(:)), [1; 3]);
            testCase.verifyTrue(all(strcmpi(T_F0.('x.properties.response_type'), 'mean')));

            % T_F1: the 'F1' rows of the two paired cells (ids 2 and 4).
            testCase.verifyEqual(sort(T_F1.id(:)), [2; 4]);
            testCase.verifyTrue(all(strcmpi(T_F1.('x.properties.response_type'), 'F1')));

            % T_out is unchanged by requesting the extra outputs.
            testCase.verifyEqual(sort(T_out.id(:)), [2; 3]);

            % F0 and F1 sets must be disjoint and exclude the unpaired row.
            testCase.verifyEmpty(intersect(T_F0.id, T_F1.id));
            testCase.verifyFalse(ismember(5, T_F0.id));
        end

        function testSingleFilterMode(testCase)
            % The single-filter code path produces the same results.
            T = ndi.unittest.viz.testSelectRowsByModulationIndex.makeTable();
            [T_out, T_F0, T_F1] = ndi.viz.selectRowsByModulationIndex(T, ...
                'allowMultipleModulationIndexFilters', false);

            testCase.verifyEqual(sort(T_out.id(:)), [2; 3]);
            testCase.verifyEqual(sort(T_F0.id(:)), [1; 3]);
            testCase.verifyEqual(sort(T_F1.id(:)), [2; 4]);
        end

        function testBackwardCompatibleSingleOutput(testCase)
            % Single-output callers continue to receive the selected table.
            T = ndi.unittest.viz.testSelectRowsByModulationIndex.makeTable();
            T_out = ndi.viz.selectRowsByModulationIndex(T);
            testCase.verifyEqual(height(T_out), 2);
        end

    end
end
