## \namespace orbit
## \brief The base module for all pyORBIT core python
## modules.
##
## Some of the modules:
## - lattice     - the base lattice implementation.
## - parsers     - collection of parsers.
## - teapot      - the TEAPOT-like ORBIT lattice elements.
## - teapot_base - fundamental bricks of teapot elements.
## - utils       - the collection of useful Python utilites.

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
