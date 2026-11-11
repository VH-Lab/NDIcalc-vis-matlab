# NDIcalc-vis-matlab

NDI calculator objects and documents for vision research.

## Getting started

- [Installation instructions](docs/installation.md)
- [Documentation](https://vh-lab.github.io/NDIcalc-vis-matlab/)

## Manuals

- [NDI Pipeline Manager & Calculator GUI User Manual](manuals/NDI_Pipeline_Manager_User_Manual.md)
- [NDI Calculator Manual (in NDI-matlab)](https://github.com/VH-Lab/NDI-matlab/blob/main/src/ndi/docs/NDI-matlab/manual/developer/NDI_Calculator_Manual.md) — canonical guide to writing `ndi.calculator` objects.

## Example dataset and demonstrations

The paper that describes these *calculators* (Lepsky et al., "A motif for
domain-specific analysis applets that are easy to learn, reuse, test, and to
compose into pipelines: application to vision science") is accompanied by an
example dataset in NDI Cloud,
[doi:10.63884/ndic.2026.i1iqjm6a](https://doi.org/10.63884/ndic.2026.i1iqjm6a).
It holds two recording sessions, one from Griswold and Van Hooser (2025) and
one from the study of Suárez Casanova et al., together with the output
documents the *calculators* wrote for them.

The package `ndi.calc.paperFigures` downloads that dataset and shows the
*calculators* at work on it:

```matlab
% download the dataset into userpath (or open it, if it is already there)
D = ndi.calc.paperFigures.openDataset();

% plot the stored results of the example cells, for every calculator
ndi.calc.paperFigures.plotResults(D);

% or only some of them
ndi.calc.paperFigures.plotResults(D, 'doContrast', false, 'doHartley', false);

% run a calculator again on the example cells' recordings and plot the result;
% nothing is saved. Name each calculator to recompute (all are off by default);
% Hartley takes minutes per cell.
ndi.calc.paperFigures.recomputeResults(D, 'doOridir', true);
ndi.calc.paperFigures.recomputeResults(D, 'doSpeed', true, 'doHartley', true);
```

| Function | What it does |
|---|---|
| `openDataset` | Downloads the dataset (documents only; data files are fetched when first needed), or opens it if `<userpath>/<NDI Cloud id>` already holds it. |
| `plotResults` | Plots the stored output of each *calculator* for the example cells: orientation and direction, contrast, contrast sensitivity, and spatial and temporal frequency tuning from the Griswold and Van Hooser session (16 plots per window), and speed tuning and Hartley receptive fields from the Suárez Casanova session (one window each, the last Hartley epoch per cell). |
| `recomputeResults` | Runs the chosen *calculators* again on the same cells and plots the new results in windows named "Recomputed: …", without saving them. |
| `exampleCells` | The example cells, by element name and reference, for each session. |

A recomputed result need not match the stored one digit for digit: the fits
start from random points. That each *calculator* recovers known answers is shown
by its self-tests, which run on synthetic data with known parameters.

## License

[MIT](LICENSE). Note that NDI-matlab, which these *calculators* require, is
distributed under its own license (CC BY-NC-SA 4.0).
