from orbit.core._module_loader import load_module

mod_name = __name__.split(".")[-1]

locals()[mod_name] = load_module(mod_name)

del mod_name
