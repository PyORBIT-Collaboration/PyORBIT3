## Dynamically import all submodules
import os
import importlib

submodules = [f.name for f in os.scandir(os.path.dirname(__file__)) if f.is_dir() and not f.name.startswith(".")]

for module in submodules:
    _mod = importlib.import_module(f"orbit.core.{module}", package=None)
    locals().update({module: _mod})

del _mod
