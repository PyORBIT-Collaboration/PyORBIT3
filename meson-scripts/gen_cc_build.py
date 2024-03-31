# This file is outdated, left here for reference.

from pathlib import Path

sources = Path(__file__).parent.parent.parent / 'src'
base_src = '../../src/'

files = []
headers = set()

headers.add('main')
headers.add('utils/ellipticalint')

for file_path in sources.glob('**/*.cc'):
    f = file_path.relative_to(sources)
    files.append(str(f))
    headers.add(str(f.parent))

files_merged = f"',\n\t'{base_src}".join(files)
headers_merged = f"',\n\t'{base_src}".join(headers)

buffer = \
f"""
sources = files([
    '{base_src}{files_merged}'
]) 
headers = include_directories([
    '{base_src}{headers_merged}'
    
])

inc += headers

core_lib = library('core',
    sources: sources,
    include_directories: inc,
    cpp_args: ['-fPIC', '-std=c++11'],    
    install: false,

)

core_dep = declare_dependency(link_with : core_lib)
    
"""


print(buffer)



