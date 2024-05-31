
subdir('sc1d')
subdir('directforce2p5d')
subdir('sc2dslicebyslice')
subdir('sc3d')
subdir('sc2p5d')


py_sources = files([
    'scLatticeModifications.py',
	'scAccNodes.py',
	'__init__.py',
	'sccenteredLatticeModifications.py'
])

python.install_sources(
    py_sources,
    subdir: 'orbit/space_charge',
    # pure: true,
)        
