import os
import sys
import subprocess
import argparse
from pathlib import Path
from distutils.ccompiler import new_compiler


import platform
p = platform.uname().system
LIB_TYPE = 'shared'
DEFAULT_LIBS = ['m']
SEARCH = ['/usr']

if p == 'Darwin':
    LIB_TYPE = 'dylib'
    DEFAULT_LIBS = []
    SEARCH = ['/opt/homebrew']


def parse_options(line):
    library_dirs = []
    libraries = []
    include_dirs = []
    compile_options = []
    for option in line.split(' ')[1:]:
        if option.startswith('-L'):
            library_dirs.append(option[2:])
        elif option.startswith('-l'):
            libraries.append(option[2:])
        elif option.startswith('-I'):
            include_dirs.append(option[2:])
        else:
            compile_options.append(option[2:])

    return library_dirs, libraries, include_dirs, compile_options


def find_mpi(places_to_look=SEARCH):
    if not places_to_look:
        places_to_look = []
    for place in places_to_look:
        for p in Path(place).glob('**/mpicc'):
            mpi_compiler = p
            print(f"Found MPI compiler at {mpi_compiler}")
            return mpi_compiler
    return None

def check_library(name, include_dir=None, lib_dir=None, test_file=None):
    old_path = Path.cwd()
    comp_dir = test_file if test_file else Path(__file__).parent
    os.chdir(comp_dir)

    f_name = f'test_{name}'
    src = f'{f_name}.c'
    obj = f'{f_name}.o'
    try:
        compiler = new_compiler()
        if include_dir:
            compiler.add_include_dir(str(include_dir))
        if lib_dir:
            compiler.add_library_dir(str(lib_dir))
        [compiler.add_library(n) for n in DEFAULT_LIBS]
        compiler.add_library(name)

        compiler.compile([src])
        compiler.link_executable([obj], f_name)
        return True
    except Exception as e:
        print(f"Library {name} fails compilation.")
        return False
    finally:
        if Path(f_name).exists():
            os.remove(f_name)
        if Path(obj).exists():
            os.remove(obj)
        os.chdir(old_path)


def find_library(name, places_to_look=SEARCH, test_file=None, lib_type=LIB_TYPE):
    if not places_to_look:
        places_to_look = []
    compiler = new_compiler()
    library = compiler.library_filename(name, lib_type=lib_type)
    for place in places_to_look:
        for p in Path(place).glob(f'**/{library}'):
            lib_dir = p.parent
            include_dir = lib_dir.parent / 'include'
            if not include_dir.exists():
                continue
            print(f"Found {name} library at {p}")
            if check_library(name, include_dir, lib_dir):
                return str(lib_dir), str(include_dir)

    return None


def main():


    epilog_text = '''Examples: 
           python create_env.py .venv
           # will create a VENV using default python in .venv
           # update the pip.conf and patch activation script for use in ICS                   
           '''
    parser = argparse.ArgumentParser(epilog=epilog_text, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--python', default=None, type=str,
                        help='Python interpreter to use')
    parser.add_argument('--mpi', default=None, type=str,
                        help='Path to MPI installation')

    parser.add_argument('venv_name', default='.venv', type=str, nargs='?',
                        help='Name of virtual environment')

    args = parser.parse_args()

    if args.mpi:
        SEARCH.insert(0, args.mpi)

    venv_name = args.venv_name
    python = args.python if args.python else sys.executable

    path = Path(venv_name)
    if path.exists():
        print(f'{venv_name} already exists.')
        return

    #  create virtual environment
    try:

        result = subprocess.run([python, '-m', 'venv', venv_name])
        if result.returncode != 0:
            print(f'Failed to create {venv_name}.')
            return
        print(f'Virtual environment {venv_name} created.')

        mpi = find_mpi()
        fftw3 = find_library('fftw3')

        def insert_after(origin, marker, new_lines):
            for i, l in enumerate(origin):
                if l.startswith(marker):
                    break
            for l in reversed(new_lines):
                origin.insert(i+1, l)

            return origin

        p = path / 'bin' / 'activate'
        with open(p, "r") as f:
            contents = f.readlines()

        export = []
        unexport = []
        if mpi:
            export.append(f'export MPICC={mpi}\n')
            export.append(f'export PATH={str(mpi.parent)}:$PATH\n')
            unexport.append('\tunset MPICC\n')

        if fftw3:
            export.append(f'export FFTW3_LIB_DIR={fftw3[0]}\n')
            export.append(f'export FFTW3_INCLUDE_DIR={fftw3[1]}\n')
            unexport.append('\tunset FFTW3_LIB_DIR\n')
            unexport.append('\tunset FFTW3_INCLUDE_DIR\n')

        export.insert(0, '\n# Start PyORBIT patch\n')
        export.append('# End PyORBIT patch\n\n')
        unexport.insert(0, '\n\t# Start PyORBIT patch\n')
        unexport.append('\t# End PyORBIT patch\n\n')
        contents = insert_after(contents, 'PATH="$VIRTUAL_ENV/bin:$PATH"', export)
        contents = insert_after(contents, 'deactivate () {', unexport)



        with open(p, "w") as f:
            f.write("".join(contents))
        print(f'Bash activation script patched.')

        # show python version
        python = path / 'bin' / 'python'
        result = subprocess.run([python, '--version'], stdout=subprocess.PIPE)
        if result.returncode != 0:
            print(f'Failed to launch Python interpreter.')
            return
        print(f'Interpreter version is {str(result.stdout, "utf-8")}', end='')

        #  update pip
        pip = path / 'bin' / 'pip'
        result = subprocess.run([pip, 'install', '-U', 'pip', 'setuptools', 'setuptools_scm', 'pytest'], stdout=subprocess.PIPE)
        if result.returncode != 0:
            print(f'Failed to update pip .')
            return
        print(f'pip updated.')


        result = subprocess.run([pip, 'list'], stdout=subprocess.PIPE)
        if result.returncode != 0:
            print(f'Failed to list packages.')
            return

        print(f'Installed packages:\n{str(result.stdout, "utf-8")}')
        print(f'You can activate it by running "source {venv_name}/bin/activate"')

    except Exception as e:
        print(f'Failed to create virtual environment: {e}')


if __name__ == "__main__":
    main()
