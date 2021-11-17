# Bioplatforms Australia Submission Generator

Generate metadata suitable for submissions to international data
repositories, from the metadata contained in the BPA Metadata Portal.

Example setup with python virtual environment in linux home dir:
```
python -m venv ~/.virtual/bpa-submission-generator
. ~/.virtual/bpa-submission-generator/bin/activate
python setup.py install
pip install -r requirements.txt
```

Example usage:
bpa-submit -k <ckan-api-key> -u https://data.bioplatforms.com base-ncbi 2>&1 | tee base.log
