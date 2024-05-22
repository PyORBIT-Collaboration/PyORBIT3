. /opt/conda/etc/profile.d/conda.sh
conda env create -n po3 --file environment.yml
conda activate po3
pip install -U meson-python setuptools setuptools-scm
pip install --no-build-isolation --editable .