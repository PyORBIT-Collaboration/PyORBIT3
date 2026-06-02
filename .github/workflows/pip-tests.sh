. .po3/bin/activate
pip install pytest
python -m pytest tests/py/ -v
cd examples/SNS_Linac/pyorbit3_linac_model/
python pyorbit3_sns_linac_mebt_hebt2.py
