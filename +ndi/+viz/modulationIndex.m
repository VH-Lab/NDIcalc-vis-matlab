function T = modulationIndex(T, NameValueArgs)
%MODULATIONINDEX Adds a modulation index column to a table based on paired responses.
%   T_OUT = NDI.VIZ.MODULATIONINDEX(T, NAMEVALUEARGS) calculates a
%   modulation index and adds it as a new column to the input table T.
%
%   DESCRIPTION:
%   This function identifies pairs of rows corresponding to 'mean' (unmodulated)
%   and 'F1' (modulated) responses and computes an index (2 * X.TC.F1/(X.TC.F1+X.TC.mean)).
%   The result is added as a new column to the input table. The column name is derived
%   from the tuning curve data column (e.g., if data is from 'X.TC.mean',
%   the new column will be 'X.TC.modulationIndex').
%
%   For each calculated index, the value is assigned to BOTH rows in the
%   'mean'/'F1' pair. All other rows not part of a successful pair calculation
%   will have a NaN value in this new column.
%
%   The function performs several steps:
%   1.  Identifies columns for response type ('mean'/'F1') and tuning curve data.
%   2.  Groups rows based on identifying columns (e.g., element ID, session ID).
%   3.  Pairs 'mean' and 'F1' rows within each group using the closest timestamp.
%   4.  Calculates the modulation index for each pair and populates the new column.
%
%   INPUTS:
%   T - A MATLAB table containing the data. This table will be modified.
%
%   NAME-VALUE PAIRS:
%   'restrictToZeroTwo' - (Optional) Logical. If true (default), the calculated
%       modulation index is restricted to the range [0, 2].
%   'responseTypeColumn' - (Optional) Name of the column that indicates the
%       response type ('mean' or 'F1'). Searched automatically if not provided.
%   'TuningCurveMeanColumn' - (Optional) Name of the column with tuning curve data.
%       Searched automatically if not provided.
%   'rowPairConstantColumns' - (Optional) Cell array of columns that must have
%       identical values to form a group.
%       Default: {'depends_on_element_id', 'base.session_id', ...
%                 'document_class.class_name', 'document_class.class_version'}
%
%   OUTPUTS:
%   T - The input table with one new column, `X.TC.modulationIndex`, added.
%
%   EXAMPLE:
%   % Assume 'myData' is a table structured like the example in the prompt.
%   myData_with_MI = ndi.viz.modulationIndex(myData);
%

arguments
    T table
    NameValueArgs.restrictToZeroTwo (1,1) logical = true
    NameValueArgs.responseTypeColumn (1,:) char = ''
    NameValueArgs.TuningCurveMeanColumn (1,:) char = ''
    NameValueArgs.rowPairConstantColumns (1,:) cell = ...
        {'depends_on_element_id', 'base.session_id', 'document_class.class_name', 'document_class.class_version'}
end

% --- Task 1 & 2: Find necessary columns ---
% FIX: The default search suffix is now 'properties.response_type'
responseTypeColumn = find_column_by_suffix(T, NameValueArgs.responseTypeColumn, 'properties.response_type');
TuningCurveMeanColumn = find_column_by_suffix(T, NameValueArgs.TuningCurveMeanColumn, 'TC.mean');

% --- Initialize the new output column ---
% Derive the new column name from the tuning curve data column name
new_col_name = strrep(TuningCurveMeanColumn, '.mean', '.modulationIndex');
% Add the column to the table and initialize all values to NaN
T.(new_col_name) = nan(height(T), 1);

% Ensure the datestamp column for pairing exists
if ~ismember('base.datestamp', T.Properties.VariableNames)
    error(['ndi.viz:modulationIndex:MissingColumn  ' ...
        'The required column "base.datestamp" is not present in the table.']);
end

% --- Task 3: Find sets of row pairs by grouping ---
missing_cols = setdiff(NameValueArgs.rowPairConstantColumns, T.Properties.VariableNames);
if ~isempty(missing_cols)
    error('ndi:viz:modulationIndex:MissingGroupingColumn', ...
        'The following specified grouping columns are missing from the table: %s', strjoin(missing_cols, ', '));
end

% Create a temporary table for grouping, casting all variables to string for robustness
grouping_vars_table = table();
for i = 1:numel(NameValueArgs.rowPairConstantColumns)
    col_name = NameValueArgs.rowPairConstantColumns{i};
    grouping_vars_table.(col_name) = string(T.(col_name));
end

[G, ~] = findgroups(grouping_vars_table);
num_groups = max(G);

% --- Task 4 & 5: Process each group to find pairs and calculate index ---
for i = 1:num_groups
    group_indices = find(G == i);
    group_table = T(group_indices, :);

    is_mean = strcmpi(group_table.(responseTypeColumn), 'mean');
    is_F1 = strcmpi(group_table.(responseTypeColumn), 'F1');

    mean_indices_in_group = find(is_mean);
    F1_indices_in_group = find(is_F1);

    if isempty(mean_indices_in_group) || isempty(F1_indices_in_group)
        continue; % Skip group if it doesn't have at least one of each type
    end

    % Identify the best pair based on closest datestamp
    best_mean_idx_in_group = -1;
    best_F1_idx_in_group = -1;

    if isscalar(mean_indices_in_group) && isscalar(F1_indices_in_group)
        best_mean_idx_in_group = mean_indices_in_group(1);
        best_F1_idx_in_group = F1_indices_in_group(1);
    else
        min_time_diff = inf;
        mean_dates = datetime(group_table.('base.datestamp')(mean_indices_in_group), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSS''Z''', 'TimeZone', 'UTC');
        F1_dates = datetime(group_table.('base.datestamp')(F1_indices_in_group), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSS''Z''', 'TimeZone', 'UTC');
        
        for m = 1:numel(mean_indices_in_group)
            for f = 1:numel(F1_indices_in_group)
                time_diff = abs(mean_dates(m) - F1_dates(f));
                if time_diff < min_time_diff
                    min_time_diff = time_diff;
                    best_mean_idx_in_group = mean_indices_in_group(m);
                    best_F1_idx_in_group = F1_indices_in_group(f);
                end
            end
        end
    end
    
    if best_mean_idx_in_group == -1
        continue; % Could not form a pair
    end
    
    % --- Calculate Modulation Index for the pair ---
    mean_data_cell = group_table.(TuningCurveMeanColumn)(best_mean_idx_in_group);
    F1_data_cell = group_table.(TuningCurveMeanColumn)(best_F1_idx_in_group);

    if isempty(mean_data_cell{1}) || isempty(F1_data_cell{1})
        continue; % Skip if data is missing
    end
    
    mean_data = mean_data_cell{1};
    F1_data = F1_data_cell{1};

    [mean_max_val, mean_max_idx] = max(mean_data);
    [F1_max_val, F1_max_idx] = max(F1_data);
    
    mi = NaN; % Default to NaN

    if F1_max_val > mean_max_val
        numerator = 2 * F1_max_val;
        denominator = F1_max_val + mean_data(F1_max_idx);
    else % mean value is greater or equal
        numerator = 2 * F1_data(mean_max_idx);
        denominator = F1_data(mean_max_idx) + mean_max_val;
    end
    
    if denominator ~= 0
        mi = numerator / denominator;
    elseif mean_max_val == 0 && F1_max_val == 0
        mi = 0;
    end
    
    if NameValueArgs.restrictToZeroTwo
        if mi < 0
            mi = 0;
        elseif mi > 2
            mi = 2;
        end
    end

    % --- Update the main table with the calculated value for BOTH rows ---
    % Get the original row indices from the full table T
    original_mean_idx = group_indices(best_mean_idx_in_group);
    original_F1_idx = group_indices(best_F1_idx_in_group);
    
    % Assign the same MI value to both rows in the pair
    T.(new_col_name)(original_mean_idx) = mi;
    T.(new_col_name)(original_F1_idx) = mi;
end

end

% --- Helper Function to find a column name ---
function col_name = find_column_by_suffix(T, specified_name, suffix)
    if ~isempty(specified_name)
        if ~ismember(specified_name, T.Properties.VariableNames)
            error('ndi:viz:modulationIndex:ColumnNotFound', ...
                'The specified column "%s" was not found in the table.', specified_name);
        end
        col_name = specified_name;
        return;
    end

    % Convert both names and suffix to lowercase for case-insensitive matching.
    candidate_cols = T.Properties.VariableNames(endsWith(lower(T.Properties.VariableNames), lower(suffix)));

    if isscalar(candidate_cols)
        col_name = candidate_cols{1};
    elseif numel(candidate_cols) > 1
        error('ndi:viz:modulationIndex:AmbiguousColumn', ...
            'Found multiple columns ending in "%s": %s. Please specify one using the appropriate Name-Value argument.', ...
            suffix, strjoin(candidate_cols, ', '));
    else % 0 candidates
        error('ndi:viz:modulationIndex:ColumnNotFound', ...
            'Could not automatically find a column ending in "%s". Please specify the column name manually.', suffix);
    end
end