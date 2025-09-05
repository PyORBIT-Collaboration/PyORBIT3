ext_sources = files([
    'DanilovEnvelopeSolver20.cc',
    'DanilovEnvelopeSolver22.cc',
    'wrap_DanilovEnvelopeSolver20.cc',
    'wrap_DanilovEnvelopeSolver22.cc',
    'wrap_danilov_envelope.cc',
    'danilov_envelope.cc',
])


python.extension_module('danilov_envelope',
    sources: ext_sources,
    include_directories: inc,
    cpp_args: ['-fPIC', '-std=c++11'],
    dependencies: [core_dep],
    install: true,
    subdir: 'orbit/ext',
)
