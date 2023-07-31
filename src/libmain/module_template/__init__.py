"""
This is not a module, but a template used by setup.py to generate boilerplate code for each one of the pyORBIT C++ extension modules.

All C++ extension modules MUST be listed in the setup.py file to be able to use them.

This init file load a single extension module and makes it available for the parent package, i.e., orbit.core.
"""

from orbit.core._module_loader import load_module

mod_name = __name__.split(".")[-1]

locals().update(load_module(mod_name).__dict__)

del mod_name
