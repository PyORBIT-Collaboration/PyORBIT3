import os

from pathlib import Path
from queue import Queue

root = Path('py')
f_q = Queue[Path]()
f_q.put(root)

src_split = ',\n\t'
while not f_q.empty():
    path = f_q.get()
    pys = [f"'{p.name}'" for p in path.iterdir() if p.is_file() and p.name.endswith('.py')]
    sub = [p for p in path.iterdir() if p.is_dir()]
    subdirs = [f"subdir('{p.name}')\n" for p in sub]
    [f_q.put(p) for p in sub]
    contents = f"""
{''.join(subdirs)}
"""
    if pys:
        contents += f"""
py_sources = files([
    {src_split.join(pys)}
])

python.install_sources(
    py_sources,
    subdir: '{path.relative_to(root)}',
    # pure: true,
)        
"""
    print(contents)
    with open(path / 'meson.build', 'w') as f:
        f.write(contents)
