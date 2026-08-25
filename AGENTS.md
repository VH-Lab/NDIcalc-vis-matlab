# Instructions for AI agents working in this repository

## Primary language

MATLAB. Running or testing anything here needs a licensed MATLAB. If you do not
have one, do not try to execute `.m` files: make the change, push it to a
branch, and say plainly that it has not been run. CI runs the real suite on
every pull request, so a push is how a change here gets verified.

## Self-tests, seeds, and stored expectations

The seven *calculator* self-tests build synthetic inputs from known parameters
and compare the result against a stored expected document, field by field,
within tolerances declared in `+ndi/+calc/+vis/mock/<calculator>/mock.N.compare.json`.

**The synthetic inputs and the fit start points are deliberately not seeded.**
They draw from a private stream (so they do not disturb a caller's own random
numbers) but a fresh draw on every run. A test compared against one pinned draw
would show only that the calculation is deterministic; unpinned, passing it
shows the *calculator* recovering its answer from data it has not seen before.
Do not reintroduce a fixed seed to make a test pass.

**A changed seed is not a reason to regenerate the stored expectations.** The
comparison is `|actual - expected| < tolerance`, and the expectation is just a
value. If a quantity is recovered to well inside its tolerance on any draw,
then any draw's stored value serves as well as any other. Regeneration is
needed when the *code* changes what a calculator computes -- not when the
randomness changes.

So when a self-test goes out of tolerance, the question is never "which seed
produced the expectation". It is whether that field's spread across draws
exceeds its declared tolerance, and if so, whether the tolerance or the
comparison method is wrong.

### Tolerances on the constrained fits

`speed_tuning` runs three fits for a nested F test: `fit` with the speed
exponent xi free, `fit_fullspeed` with xi forced to 1, and `fit_no_speed` with
xi forced to 0. On any given mock at least one constrained arm is fitting a
model that is wrong for the data -- that is what makes the F test informative.

A mis-specified fit has no interior optimum. On self-test 7, whose cell is not
speed tuned, `fit_fullspeed` drives the preferred temporal frequency to exactly
32, the upper bound, while the free-xi fit recovers the true 16 with a sum of
squared errors of zero. A parameter pinned to a bound is not an estimate; it is
the fit saying "as high as you will let me", and different optimizer starts stop
at different points along that flat ridge.

So the constrained arms' `Priebe_fit_values` and frequency preferences carry
deliberately wide tolerances, set to span the parameter's allowed range
(sf0 in [0.05, 1.2], tf0 in [0.5, 32], responses in [0, 2*max(R)]). What is
still compared tightly is what those fits are *for*: `r_squared`, `partial_r2`,
and the nested F test's p-value, compared by the decision it supports rather
than its digits. Those reproduce fine. The unconstrained `fit` is well specified
and keeps its tight tolerances -- do not loosen those.

## Related repositories

`VH-Lab/NDI-matlab` holds the framework, including the mock generation these
self-tests use (`ndi.calc.tuning_fit`, `ndi.mock.fun.stimulus_response`). A
change to how mocks are built usually belongs there, not here. The paper that
describes these calculators is in `VH-Lab-Biz/2026-Calculators-Paper`.
