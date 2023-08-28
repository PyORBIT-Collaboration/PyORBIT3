import importlib.machinery
import importlib.util

from orbit.core import _orbit


def load_module(mod_name: str):
    # The name of the module in the C++ shared library starts with an underscore.
    mod_name = "_" + mod_name

    pkg_path = _orbit.__file__

    loader = importlib.machinery.ExtensionFileLoader(mod_name, pkg_path)
    spec = importlib.util.spec_from_file_location(mod_name, pkg_path)
    module = importlib.util.module_from_spec(spec)
    loader.exec_module(module)

    return module
