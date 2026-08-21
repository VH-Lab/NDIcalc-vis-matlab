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
else, or under a different name, the calculators load but their self-tests fail
while reading the JSON schema that defines each one's output type:

```
DID:Document:readjsonfilelocation: found no match for {calculator}
```

Names that do not match are easier to end up with than renaming suggests:
GitHub's **Download ZIP** expands to `NDIcalc-vis-matlab-main`, and an archived
release expands to `NDIcalc-vis-matlab-1.0.0`. Neither matches.

To find out whether the copy you have can be discovered, and what to do if it
cannot, run

```matlab
ndi.fun.checkCalcDirectory();
```

It reports the checkout, the directory NDI searches, and the exact command to
link the two. It changes nothing on disk — move or link the folder yourself. The
test suite calls it first, so a checkout in the wrong place fails with that
message rather than with the schema error above.
