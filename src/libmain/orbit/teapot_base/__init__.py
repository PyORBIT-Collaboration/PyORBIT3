from orbit.core._module_loader import load_module

mod_name = __name__.split(".")[-1]

locals().update(load_module(mod_name).__dict__)

del mod_name
