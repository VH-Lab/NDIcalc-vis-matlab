# Notes on `generate_mock_docs` in `ndi.calc.vis.*` objects

The `generate_mock_docs` method is a testing utility present in NDI visualization calculator classes (e.g., `contrast_tuning`, `oridir_tuning`, `spatial_frequency_tuning`, `speed_tuning`, `temporal_frequency_tuning`). It generates synthetic input data (mock documents) and runs the calculator to produce actual outputs, which can then be compared against expected outputs.

## Function Signature

```matlab
[docs, doc_output, doc_expected_output] = generate_mock_docs(obj, scope, number_of_tests, Name, Value)
```

## Inputs

| Input Argument | Type | Description |
| :--- | :--- | :--- |
| `obj` | `ndi.calculator` | The instance of the calculator object (e.g., `contrast_tuning_obj`). |
| `scope` | String | Defines the testing conditions, primarily affecting noise levels and repetition counts. Allowed values:<br>- `'highSNR'`: High signal-to-noise ratio (standard conditions).<br>- `'lowSNR'`: Low signal-to-noise ratio (high noise). |
| `number_of_tests` | Integer | The number of distinct test cases to generate. The method will loop from 1 to `number_of_tests`, generating different parameters for each iteration via `generate_mock_parameters`. |
| `Name, Value` | Key/Value pairs | Optional parameters to modify behavior. These are defined in an `arguments` block. |

### Optional Name-Value Arguments

| Argument | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `generate_expected_docs` | logical | `false` | If `true`, the method saves the current output as the "expected" output for future tests. Use this when updating the calculator logic or creating new tests. |
| `specific_test_inds` | double vector | `[]` | Allows specifying a subset of test indices to run (e.g., `[1 5]`). If empty, all `number_of_tests` are run. |

## Outputs

| Output Argument | Type | Description |
| :--- | :--- | :--- |
| `docs` | Cell Array | A 1xN cell array (where N is number of tests). `docs{i}` is a cell array containing the **input documents** generated for the `i`-th test. These mock documents (typically stimulus response documents) simulate the raw data in the database. |
| `doc_output` | Cell Array | A 1xN cell array. `doc_output{i}` contains the **actual output document** produced by the calculator when processing the inputs in `docs{i}`. |
| `doc_expected_output` | Cell Array | A 1xN cell array. `doc_expected_output{i}` contains the **expected output document**. If `generate_expected_docs` is `false`, this is loaded from a stored file. If `true`, it mirrors `doc_output{i}` (and is written to disk). |

## General Workflow

1.  **Parameter Generation**: For each test index, the method calls `generate_mock_parameters(scope, i)` to determine the underlying properties of the synthetic data (e.g., optimal contrast, tuning width, orientation preference).
2.  **Data Simulation**: It generates synthetic response data (e.g., firing rates) based on these parameters and the noise level defined by `scope`.
3.  **Input Creation**: It creates mock NDI documents (e.g., using `ndi.mock.fun.stimulus_response`) that encapsulate this synthetic data.
4.  **Execution**: It configures the calculator to use these mock inputs and runs it (`obj.run`).
5.  **Output & Validation**: It returns the inputs, the actual calculator output, and the expected output (for comparison using `compare_mock_docs`).
