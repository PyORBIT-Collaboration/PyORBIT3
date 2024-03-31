
subdir('fitting')
subdir('xml')
subdir('orbit_mpi_utils')


py_sources = files([
    'phaseOperations.py',
	'TypedObject.py',
	'__init__.py',
	'ParamsDictObject.py',
	'consts.py',
	'multiDimArray.py',
	'NamedObject.py',
	'orbitFinalize.py'
])

python.install_sources(
    py_sources,
    subdir: 'orbit/utils',
    # pure: true,
)        
