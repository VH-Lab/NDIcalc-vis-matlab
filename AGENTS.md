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

### The spline fits

`spatial_frequency_tuning` and `temporal_frequency_tuning` each report three
fits, and one of them, `fit_spline`, is not a fit at all: a spline interpolates
its inputs, so it passes exactly through the data and has no free parameters.
It therefore cannot be more reproducible than the data are, and with the mock
inputs drawn fresh on every run its curve moves on every run by construction.

It is also badly conditioned. The knots are 100 points of
`logspace(-2, log10(60))`, whose narrowest gap is 0.0009 and whose widest is
5.05 -- a ratio of about 5500 to 1. A cubic spline's second-derivative system
at that spacing amplifies a small change in the data into a large excursion
between knots, so the 0.1 percent noise the mocks carry does not stay 0.1
percent in the interpolant. `spatial_frequency_tuning.m` already says as much,
in the comment that keeps the spline out of the plot: *the spline fits are
terrible*.

So `fit_spline.fit` carries a deliberately wide tolerance, 100, against
response amplitudes of at most 20. It is a presence-and-sanity check, not an
accuracy check: it still fails on a missing field, an empty one, or an `Inf`,
and it no longer fails because a different noise draw moved an interpolant.

It is a wide tolerance rather than `none` on purpose. In
`ndi.database.doctools.docComparison`, `'none'` is the first branch of the
comparison loop and `continue`s before the value is read, so it skips the
field-missing check as well -- a field set to `none` is not loosely checked,
it is invisible. `fit_spline.values`, `L50`, `Pref`, `H50` and `bandwidth` are
all `none` already, and `R2` is identically 1 for an interpolant, so a wide
number on `fit` is the only thing keeping that group wired up at all.

## Related repositories

`VH-Lab/NDI-matlab` holds the framework, including the mock generation these
self-tests use (`ndi.calc.tuning_fit`, `ndi.mock.fun.stimulus_response`). A
change to how mocks are built usually belongs there, not here. The paper that
describes these calculators is in `VH-Lab-Biz/2026-Calculators-Paper`.
