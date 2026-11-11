# ndicalc-vis-matlab (Python packaging)

This folder lets NDI-python install the document definitions of
[NDIcalc-vis-matlab](https://github.com/VH-Lab/NDIcalc-vis-matlab) with pip or
uv. It contains no calculators; the calculators are MATLAB code and MATLAB
users do not need anything here.

The wheel contains the `ndicalc_vis_matlab` package and a copy of the
repository's `ndi_common/database_documents` and `ndi_common/schema_documents`
JSON trees.

```python
import ndicalc_vis_matlab
ndicalc_vis_matlab.toolbox_dir()   # directory that contains ndi_common/
```

It is also advertised as the entry point `NDIcalc-vis-matlab` in the group
`ndi.calc`, so NDI-python can find it with
`importlib.metadata.entry_points(group="ndi.calc")`.

Install from GitHub without git:

```
pip install https://github.com/VH-Lab/NDIcalc-vis-matlab/archive/refs/heads/main.zip
```

Run the tests from a checkout:

```
pip install . pytest
pytest python/tests
```
