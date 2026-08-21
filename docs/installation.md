# Installation:

1. Install NDI

See the [NDI Installation Guide](https://vh-lab.github.io/NDI-matlab/installation/) and check out its requirements.

2. In a terminal or DOS shell, navigate to the **same parent directory** where NDI-matlab is installed. NDIcalc-vis-matlab must be a sibling of NDI-matlab so that NDI can automatically discover its calculator documents and schemas at runtime. For example, if NDI-matlab is at `/Users/yourusername/Documents/MATLAB/tools/NDI-matlab`, then you should `cd` to `/Users/yourusername/Documents/MATLAB/tools/`.

3. Clone this repository using `git clone https://github.com/VH-Lab/NDIcalc-vis-matlab`.

Your directory structure should look like:

```
tools/
├── NDI-matlab/
├── NDIcalc-vis-matlab/
└── ...
```



## If the self-tests cannot find the calculator schemas

NDI discovers this package by scanning the directory that holds NDI-matlab for
folders whose names match `NDIcalc*-matlab`. If the checkout ends up somewhere
else, or under a different name — a downloaded zip expands to
`NDIcalc-vis-matlab-main`, and continuous-integration checkouts pick their own
names — the calculators load but their self-tests fail while reading the JSON
schema that defines each one's output type:

```
DID:Document:readjsonfilelocation: found no match for {calculator}
```

Either move the checkout beside NDI-matlab under the name
`NDIcalc-vis-matlab`, as above, or run

```matlab
ndi.fun.linkCalcDirectory();
```

which creates a link of that name pointing at the current checkout. It does
nothing if the checkout can already be found, and the test suite calls it
automatically, so the self-tests now run from a checkout of any name.
