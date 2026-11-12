# Self-test expectations for `ndi.calc.vis.speed_tuning_bootstrap`

The calculator declares 4 self-tests (`numberOfSelfTests = 4`), one per mock
cell defined in `generate_mock_parameters`:

| Test | Mock cell |
|---|---|
| 1 | separable, xi = 0 |
| 2 | speed tuned, xi = 1 |
| 3 | negative ("anti-speed") coupling, xi = -0.5 |
| 4 | preferred SF (0.03 cyc/deg) below the lowest stimulus SF (0.05), xi = 0.5 |

`mock.N.compare.json` (in this folder) lists the fields compared for test N and
their tolerances. The comparisons concentrate on what the calculator is for:
xi and its interval, the p-values (compared by the decision they support,
`pValueConsistent`), and the 0/1 classifications (exact). For test 4 the
passband-derived quantities are not compared, because with the tuning peak
outside the sampled grid they are not stable; its low-SF edge flag is.

The self-tests run at 'highSNR' (noise SD 0.001 x response), where the
bootstrap interval of xi is extremely narrow. For a cell whose true xi is
exactly 0 (test 1) or exactly 1 (test 2), whether that interval excludes 0 or
includes 1 is then decided by fitting and spike-timing error and can differ
between runs, so those exact-boundary checks (p0/differentFrom0 for test 1;
p1/priebe for test 2) are not compared. The rules with a margin
(minimumEffect, equivalence, passband) are compared for those cells.

`mock.N.json`, the expected output documents, must be generated in a licensed
MATLAB session (as for every calculator here; see AGENTS.md). Until they exist,
`tests/+ndi/+unittest/+calc/+vis/testSpeedTuningBootstrap.m` skips the
self-tests (reported as an incomplete test, not a failure). To generate them,
with NDI and this package on the path:

```matlab
p = [tempname '_mock_session']; mkdir(p);
S = ndi.session.dir('test_session', p);   % an empty session, as in the unit tests
c = ndi.calc.vis.speed_tuning_bootstrap(S);
c.generate_mock_docs('highSNR', c.numberOfSelfTests, 'generate_expected_docs', true);
```

`write_mock_expected_output` refuses to overwrite an existing `mock.N.json`;
delete the old files first when the calculator's output changes on purpose.
Then commit the four `mock.N.json` files.
