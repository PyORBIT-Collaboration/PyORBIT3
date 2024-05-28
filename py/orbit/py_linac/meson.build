
subdir('linac_parsers')
subdir('lattice_modifications')
subdir('materials')
subdir('lattice')
subdir('overlapping_fields')
subdir('orbit_correction')
subdir('errors')
subdir('rf_field_readers')


py_sources = files([
    '__init__.py'
])

python.install_sources(
    py_sources,
    subdir: 'orbit/py_linac',
    # pure: true,
)        
