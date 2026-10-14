"""NDIcalc-vis document definitions for NDI-python.

This package carries a copy of the repository's ``ndi_common`` JSON trees
(``database_documents`` and ``schema_documents``) so NDI-python can find them
the same way MATLAB's ``ndi.fun.find_calc_directories`` does: it is given a
toolbox directory and appends ``ndi_common/<key>``.

NDI-python discovers this package through the ``ndi.calc`` entry point, whose
value is :func:`toolbox_dir`.
"""

from __future__ import annotations

from pathlib import Path

__all__ = ["toolbox_dir"]

_PACKAGE_DIR = Path(__file__).resolve().parent


def toolbox_dir() -> Path:
    """Return the directory that contains ``ndi_common/``.

    For a normal (wheel) install that is this package's own directory, which
    holds the packaged copy of ``ndi_common``. For an editable install from a
    checkout the packaged copy is absent, and the repository root -- two
    levels above this package, at ``<repo>/python/ndicalc_vis_matlab`` -- is
    returned instead.

    Raises:
        FileNotFoundError: if neither location contains ``ndi_common``.
    """
    candidates = (_PACKAGE_DIR, _PACKAGE_DIR.parent.parent)
    for candidate in candidates:
        if (candidate / "ndi_common").is_dir():
            return candidate
    searched = ", ".join(str(c / "ndi_common") for c in candidates)
    raise FileNotFoundError(
        "ndicalc_vis_matlab could not find its ndi_common directory; "
        f"looked in: {searched}. The installation appears to be incomplete."
    )
