# Documentation: NDIc-vis-matlab Add-on for NDI

## 1. Introduction

**NDIc-vis-matlab** is an add-on package for the Neuroscience Data Interface (NDI) framework, implemented in Matlab. It provides a suite of specialized analysis tools, packaged as **NDI Calculators** (`ndi.calculator`), primarily focused on common analyses in visual neuroscience and electrophysiology.

These calculators automate the process of analyzing data stored within an NDI session (like spike times, stimulus presentations, and basic responses) and store the detailed results, including statistical significance and model fits, back into the NDI database as new, queryable NDI documents. The add-on also typically provides built-in methods for visualizing the results of these calculations.

## 2. Core Concept: NDI Calculators

The tools within NDIc-vis-matlab are built upon the core NDI concept of the `ndi.calculator` class. Key features of this approach include:

* **Automation:** Calculators encapsulate specific analysis workflows.
* **Data Provenance:** Output documents automatically store dependencies on the input documents used (e.g., which neuron element, which stimulus presentation, which tuning curve document).
* **Standardization:** Results are stored in defined NDI document formats (specified by JSON schemas), making them consistent and queryable.
* **Modularity:** Calculators can be chained together, with the output of one serving as input for another.
* **Integration:** They operate directly on data within an NDI session object.

## 3. General Workflow with NDIc-vis-matlab Calculators

Using a calculator from this add-on typically involves the following steps in Matlab:

1.  **Instantiate the Calculator:** Create an instance of the desired calculator class, passing your NDI session object (`S`).
    ```matlab
    % Example: Instantiate the orientation/direction tuning calculator
    oridir_calc = ndi.calc.vis.oridir_tuning(S);
    ```
2.  **Prepare Input Parameters:** Determine the input parameters needed for the `calculate` method. This always includes a `depends_on` structure specifying the IDs of the input NDI documents. You might also need `input_parameters` specific to the analysis.
    * You can often use the calculator's `default_search_for_input_parameters()` and `search_for_input_parameters()` methods to automatically find suitable input documents based on default queries (e.g., finding the right type of `stimulus_tuningcurve` document).
    ```matlab
    % Example: Find potential input parameters using default query
    default_params = oridir_calc.default_search_for_input_parameters();
    % Refine query if needed, e.g., restrict to a specific element_id
    element_doc = S.database_search(ndi.query('element.name','exact_string','neuron_001','')); % Find your element
    if ~isempty(element_doc)
       default_params.query.query = default_params.query.query & ndi.query('','depends_on','element_id',element_doc{1}.id());
    end
    % Search for specific inputs matching the (potentially refined) query
    input_params_list = oridir_calc.search_for_input_parameters(default_params);
    ```
3.  **Run the Calculator:** Execute the analysis using the `run` method. You typically pass the specific input parameters structure found in the previous step. The `'Replace'` option ensures existing results for the same inputs are replaced.
    ```matlab
    % Example: Run the calculator for the first found input parameter set
    if ~isempty(input_params_list)
        output_docs = oridir_calc.run('Replace', input_params_list{1});
        disp(['Calculator generated ' int2str(numel(output_docs)) ' output document(s).']);
    else
        disp('No suitable input documents found to run the calculator.');
    end
    ```
4.  **Find Output Documents:** Locate the newly created output document(s). You can search by type (`_calc` document name) and/or by dependency on the input documents, or use the list returned directly by the `run` command. Alternatively, use the calculator's `search_for_calculator_docs` method.
    ```matlab
    % Example: Find output docs based on the input stimulus_tuningcurve_id
    stimulus_tuningcurve_id = did.db.struct_name_value_search(input_params_list{1}.depends_on,'stimulus_tuningcurve_id');
    q_out = ndi.query('','isa','oridirtuning_calc') & ndi.query('','depends_on','stimulus_tuningcurve_id',stimulus_tuningcurve_id);
    found_output_docs = S.database_search(q_out);
    ```
5.  **Examine and Visualize Results:** Load the output document and access its properties (as shown in Tutorial 5 for `oridirtuning_calc`). Use the calculator's built-in `plot` method for visualization.
    ```matlab
    % Example: Plot the result from the first output document found
    if exist('found_output_docs','var') && ~isempty(found_output_docs)
        figure; % Create a new figure
        axes_h = gca; % Get axes handle
        oridir_calc.plot(found_output_docs{1}, 'axes_handle', axes_h);
    end
    ```

## 4. Available Calculators in NDIc-vis-matlab

Based on the provided code, the following calculator classes are available (assuming they reside in an `ndi.calc.vis` namespace):

| Calculator Class Name                | Purpose                                     | Primary Input Dependency(ies)                                | Output Document Type                | Plot Method? |
| :----------------------------------- | :------------------------------------------ | :----------------------------------------------------------- | :---------------------------------- | :----------- |
| `ndi.calc.vis.contrast_sensitivity`  | Calculates contrast sensitivity curves      | `element_id`, `stimulus_presentation_id` (finds `stimulus_response_scalar`, `contrast_tuning` internally) | `contrastsensitivity_calc`        | Yes          |
| `ndi.calc.vis.contrast_tuning`       | Fits Naka-Rushton models to contrast data | `stimulus_tuningcurve_id` (contrast)                       | `contrasttuning_calc`             | Yes          |
| `ndi.calc.vis.hartley_calc`          | Performs Hartley reverse correlation        | `element_id`, `stimulus_presentation_id` (Hartley)         | `hartley_calc`                    | Yes          |
| `ndi.calc.vis.oridir_tuning`         | Analyzes orientation/direction tuning     | `stimulus_tuningcurve_id` (angle/direction)              | `oridirtuning_calc`               | Yes          |
| `ndi.calc.vis.spatial_frequency_tuning`| Analyzes spatial frequency tuning       | `stimulus_tuningcurve_id` (spatial frequency)            | `spatial_frequency_tuning_calc`   | Yes          |
| `ndi.calc.vis.speed_tuning`          | Analyzes speed tuning (SF+TF interaction) | `stimulus_tuningcurve_id` (SF & TF)                        | `speedtuning_calc`                | Yes          |
| `ndi.calc.vis.temporal_frequency_tuning`| Analyzes temporal frequency tuning      | `stimulus_tuningcurve_id` (temporal frequency)           | `temporal_frequency_tuning_calc`  | Yes          |

## 5. Key Output Document Fields (Examples)

The `_calc` documents store rich information. Refer back to the documentation snippets provided earlier for full details. Here are highlights for a few common types:

**`oridirtuning_calc` (via `orientation_direction_tuning` property):**

* **`tuning_curve`**: Mean, stderr, individual responses per direction.
* **`significance`**: p-values for visual responsiveness and tuning significance (ANOVA).
* **`vector`**: Fitless measures like preferred angles, circular variance, dot product significance.
* **`fit`**: Double Gaussian fit parameters (`Rsp`, `Rp`, `theta_pref`, `sigma`, `Rn`), derived preferred angles, tuning width (HWHH), orientation/direction indices.

**`contrasttuning_calc` (via `contrast_tuning` property):**

* **`tuning_curve`**: Mean, stderr, individual responses per contrast level.
* **`significance`**: p-values for visual responsiveness and tuning significance (ANOVA).
* **`fitless`**: Interpolated C50 value.
* **`fit`**: Parameters, fit values, R^2, empirical C50, sensitivity, etc., for 3 different Naka-Rushton model variations (RB, RBN, RBNS).

**`spatial_frequency_tuning_calc` / `temporal_frequency_tuning_calc` (via `spatial_frequency_tuning` / `temporal_frequency_tuning` properties):**

* **`tuning_curve`**: Mean, stderr, individual responses per frequency.
* **`significance`**: p-values for visual responsiveness and tuning significance (ANOVA).
* **`fitless`**: Measures like preferred frequency (Pref), low/high cutoff frequencies (L50/H50), bandwidth, low/high pass indices.
* **`fit_dog` / `fit_movshon` / `fit_movshon_c`**: Parameters, fit values, R^2, and derived Pref/L50/H50/bandwidth for different fit models (Difference-of-Gaussians, Movshon).
* **`abs`**: Repeats fitless and fit measures for the absolute value of the responses.

## 6. Associated Document Types Defined

This add-on defines the primary calculator output documents (e.g., `contrastsensitivity_calc`, `contrasttuning_calc`, `hartley_calc`, `oridirtuning_calc`, `spatial_frequency_tuning_calc`, `speedtuning_calc`, `temporal_frequency_tuning_calc`).

It also defines several intermediate or base document types that are likely used as superclasses or components within the calculator outputs (e.g., `contrast_tuning`, `hartley_reverse_correlation`, `reverse_correlation`, `spatial_frequency_tuning`, `speed_tuning`, `temporal_frequency_tuning`). Understanding these base types can help in interpreting the structure of the final `_calc` documents.

## 7. Illustrative Example (Orientation/Direction Tuning)

```matlab
% --- Assumes S is a valid NDI session object ---

% 1. Instantiate calculator
oridir_calc = ndi.calc.vis.oridir_tuning(S);

% 2. Find input parameters (e.g., for a specific element)
element_doc = S.database_search(ndi.query('element.name','exact_string','neuron_001',''));
if isempty(element_doc), error('Element not found'); end;
element_id = element_doc{1}.id();

params_query = oridir_calc.default_search_for_input_parameters();
params_query.query.query = params_query.query.query & ndi.query('','depends_on','element_id',element_id);
input_params_list = oridir_calc.search_for_input_parameters(params_query);

if isempty(input_params_list)
    error('No suitable stimulus_tuningcurve document found for this element.');
end
input_parameters = input_params_list{1}; % Use the first found input set

% 3. Run the calculator
disp('Running oridir_tuning calculator...');
output_docs = oridir_calc.run('Replace', input_parameters);
disp('Calculation complete.');

if isempty(output_docs)
    error('Calculator did not produce an output document.');
end
oridir_result_doc = output_docs{1};

% 4. (Optional) Find the output document again via search
% q_out = ndi.query('','isa','oridirtuning_calc') & ...
%         ndi.query('','depends_on','stimulus_tuningcurve_id', input_parameters.depends_on(1).value);
% found_docs = S.database_search(q_out);
% oridir_result_doc = found_docs{1};

% 5. Visualize the result
disp('Plotting results...');
figure;
axes_h = gca;
oridir_calc.plot(oridir_result_doc, 'axes_handle', axes_h);
title('Orientation/Direction Tuning Result'); % Add a custom title if needed

% 6. Examine specific fields (optional)
ot_data = oridir_result_doc.document_properties.orientation_direction_tuning;
disp(['Fit Direction Preference: ' num2str(ot_data.fit.direction_angle_preference)]);
disp(['Fit HWHH: ' num2str(ot_data.fit.hwhh)]);
