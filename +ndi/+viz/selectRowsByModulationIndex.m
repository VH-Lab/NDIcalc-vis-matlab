function T_out = selectRowsByModulationIndex(T, NameValueArgs)
% SELECTROWSBYMODULATIONINDEX - Filter rows based on the modulation index
%
%   T_OUT = NDI.VIZ.SELECTROWSBYMODULATIONINDEX(T, NAMEVALUEARGS)
%
%   Filters a table that has been produced by ndi.viz.modulationIndex.
%
%   It takes as input a table with a column that ends in the string
%   'X.TC.modulationIndex'. The user can specify the column name as an
%   optional input.
%
%   The function examines each entry with an X.TC.modulationIndex value.
%   If the modulation index is greater than or equal to 1, then the row
%   is kept in the output table if the row X.properties.response_type is
%   equal to 'F1'. If the modulation index is less than 1, then the row
%   is kept in the output table if the row X.properties.response_type is
%   equal to 'mean'.
%
%   This function can also be specified using name/value pairs.
%   'modulationIndexColumn'      : The column with the modulation index.
%                                    If not specified, the function will look for a
%                                    column that ends in '.TC.modulationIndex'.
%   'responseTypeColumn'         : The column with the response type.
%                                    If not specified, the function will look for a
%                                    column that ends in '.properties.response_type'.
%   'allowMultipleModulationIndexFilters' : (true/false) If true, then the function
%                                    can search for many columns X.TC.modulationIndex
%                                    and run the filter for each type of modulationIndex.
%                                    Default is false.
%
%   Example:
%      T_filtered = ndi.viz.selectRowsByModulationIndex(T_with_MI)
%

arguments
    T table
    NameValueArgs.modulationIndexColumn (1,:) char = ''
    NameValueArgs.responseTypeColumn (1,:) char = ''
    NameValueArgs.allowMultipleModulationIndexFilters (1,1) logical = true
end

if ~NameValueArgs.allowMultipleModulationIndexFilters
    % single filter mode
    modulationIndexColumn = find_column_by_suffix(T, NameValueArgs.modulationIndexColumn, '.TC.modulationIndex');

    if ~isempty(NameValueArgs.responseTypeColumn)
        % user provided a column name, let's use it but verify it exists
        responseTypeColumn = NameValueArgs.responseTypeColumn;
        if ~ismember(responseTypeColumn, T.Properties.VariableNames)
             error('ndi:viz:selectRowsByModulationIndex:ColumnNotFound', ...
                'The specified responseTypeColumn "%s" was not found in the table.', responseTypeColumn);
        end
    else
        % user did not provide a column name, so we have to find it.
        % First, let's try to derive it from the modulation index column name
        responseTypeColumn_derived = strrep(modulationIndexColumn, '.TC.modulationIndex', '.properties.response_type');
        if ismember(responseTypeColumn_derived, T.Properties.VariableNames)
             responseTypeColumn = responseTypeColumn_derived;
        else
             % derivation didn't work, so let's search for a column with the suffix
             responseTypeColumn = find_column_by_suffix(T, '', '.properties.response_type');
        end
    end

    MI = T.(modulationIndexColumn);
	response_type = T.(responseTypeColumn);

	keep_rows = find(...
			(MI>=1 & strcmpi(response_type,'F1')) | ...
			(MI<1 & strcmpi(response_type,'mean'))  ...
		);
	T_out = T(keep_rows,:);
else
    % multiple filter mode
    if ~isempty(NameValueArgs.modulationIndexColumn)
        warning('ndi:viz:selectRowsByModulationIndex:argignored', 'modulationIndexColumn is ignored when allowMultipleModulationIndexFilters is true.');
    end
    if ~isempty(NameValueArgs.responseTypeColumn)
        warning('ndi:viz:selectRowsByModulationIndex:argignored', 'responseTypeColumn is ignored when allowMultipleModulationIndexFilters is true.');
    end

	modulationIndexColumns = find_column_by_suffix(T, '', '.TC.modulationIndex', 'allow_multiple', true);
    keep_rows = [];
    for i=1:numel(modulationIndexColumns)
        modulationIndexColumn = modulationIndexColumns{i};
        responseTypeColumn = strrep(modulationIndexColumn,'.TC.modulationIndex','.properties.response_type');
        if ~ismember(responseTypeColumn, T.Properties.VariableNames)
            error(['Could not find response column ' responseTypeColumn ' for modulation index column ' modulationIndexColumns{i} '.']);
        end;
        MI = T.(modulationIndexColumn);
        response_type = T.(responseTypeColumn);
        keep_rows_here = find(...
                (MI>=1 & strcmpi(response_type,'F1')) | ...
                (MI<1 & strcmpi(response_type,'mean'))  ...
            );
        keep_rows = unique([keep_rows(:); keep_rows_here(:)]);
    end
    T_out = T(keep_rows,:);
end

end

function col_name = find_column_by_suffix(T, specified_name, suffix, NameValueArgs)
% FIND_COLUMN_BY_SUFFIX - helper function to find a column name by suffix
%
%  COL_NAME = FIND_COLUMN_BY_SUFFIX(T, SPECIFIED_NAME, SUFFIX, NAMEVALUEARGS)
%
%  Finds a column name in a table T that ends in SUFFIX.
%  If SPECIFIED_NAME is not empty, it is checked to see if it exists and is returned.
%  Otherwise, the function searches for a column ending in SUFFIX.
%  If not exactly 1 match is found, an error is generated.
%
%  This function accepts name/value pairs that modify its behavior:
%  'allow_multiple' (logical)   : if true, multiple matches are allowed and a cell
%                                 array of strings is returned. Default is false.
%
	arguments
		T table
		specified_name (1,:) char
		suffix (1,:) char
		NameValueArgs.allow_multiple (1,1) logical = false
	end

    if ~isempty(specified_name)
        if ~ismember(specified_name, T.Properties.VariableNames)
            error('ndi:viz:selectRowsByModulationIndex:ColumnNotFound', ...
                'The specified column "%s" was not found in the table.', specified_name);
        end
        col_name = specified_name;
        return;
    end

    candidate_cols = T.Properties.VariableNames(endsWith(lower(T.Properties.VariableNames), lower(suffix)));

	if NameValueArgs.allow_multiple,
		col_name = candidate_cols;
		return;
	end;

    if isscalar(candidate_cols)
        col_name = candidate_cols{1};
    elseif numel(candidate_cols) > 1
        error('ndi:viz:selectRowsByModulationIndex:AmbiguousColumn', ...
            ['Found multiple columns ending in "%s": %s. Please specify one using the ' ...
             'appropriate Name-Value argument or set allowMultipleModulationIndexFilters to true.'], ...
            suffix, strjoin(candidate_cols, ', '));
    else % 0 candidates
        error('ndi:viz:selectRowsByModulationIndex:ColumnNotFound', ...
            'Could not automatically find a column ending in "%s". Please specify the column name manually.', suffix);
    end
end