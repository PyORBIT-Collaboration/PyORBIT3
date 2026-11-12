try:
    from importlib import metadata
except ImportError:
    # Running on pre-3.8 Python; use importlib-metadata package
    import importlib_metadata as metadata  # type: ignore

try:
    __version__ = metadata.version("ess-python-tools")  # type: ignore
except ModuleNotFoundError:
    # package is not installed
    pass

from . import aperture
from . import bumps
from . import bunch_generators
from . import bunch_utils
from . import collimation
from . import core
from . import diagnostics
from . import errors
from . import fieldtracker
from . import foils
from . import impedances
from . import injection
from . import kickernodes
from . import lattice
from . import matching
from . import matrix_lattice
from . import orbit_correction
from . import parsers
from . import py_linac
from . import rf_cavities
from . import space_charge
from . import teapot
from . import teapot_base
from . import time_dep
from . import utils

__all__ = [
    "aperture",
    "bumps",
    "bunch_generators",
    "bunch_utils",
    "collimation",
    "core",
    "diagnostics",
    "errors",
    "fieldtracker",
    "foils",
    "impedances",
    "injection",
    "kickernodes",
    "lattice",
    "matching",
    "matrix_lattice",
    "orbit_correction",
    "parsers",
    "py_linac",
    "rf_cavities",
    "space_charge",
    "teapot",
    "teapot_base",
    "time_dep",
    "utils",
]
