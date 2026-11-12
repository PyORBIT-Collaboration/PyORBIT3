. /opt/conda/etc/profile.d/conda.sh
conda env create -n po3 --file environment-linux.yml
conda activate po3
pip install --no-build-isolation --editable .
