project('PyORBIT',
  ['cpp'],
  version: run_command('python', '-m',  'setuptools_scm').stdout().strip(),
  default_options: [
    'buildtype=release',
    'b_lto=true',
    'cpp_std=c++11',
    'warning_level=3',
  ],
)

subdir('src')
subdir('py/orbit')
