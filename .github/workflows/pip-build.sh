python3 -m venv .po3
. .po3/bin/activate
pip install -U pip
pip install -r requirements.txt
pip install -U setuptools
pip install --no-build-isolation --editable .