import os
import sys
import importlib.util

"""
This script reproduces the behabior of the legacy pyORBIT executable.
It loads the pyORBIT extension modules and makes them loadable using their short names.

When called with a script as an argument, it will execute the script.

If called with python -i flag, it will open an interactive Python interpreter
where pyORBIT modules can be directly loaded.
"""


def load_extension_modules():
    # Find submodules
    submodules = [
        f.name for f in os.scandir(os.path.dirname(__file__)) if f.is_dir() and not f.name.startswith(".") and not f.name.startswith("_")
    ]
    # Add a link to the modules without the leading underscore to be able to load them in legacy code
    for module in submodules:
        sys.modules[module] = sys.modules[f"_{module}"]


if __name__ == "__main__":
    load_extension_modules()
    if len(sys.argv) > 1:
        # Running the script as a module
        spec = importlib.util.spec_from_file_location("script", sys.argv[1])
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
