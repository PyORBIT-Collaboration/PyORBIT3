import importlib.util
import sys

import _orbit

pkg_path = _orbit.__file__


def _load_module(mod_name, path):
    spec = importlib.util.spec_from_file_location(mod_name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[mod_name] = module
    spec.loader.exec_module(module)

    return module


__all__ = ["orbit_mpi", "bunch", "teapot_base", 
           "linac", "spacecharge", "orbit_utils",
           "aperture", "foil", "collimator"]

for mod_name in __all__:
    locals()[mod_name] = _load_module(mod_name, pkg_path)
del mod_name