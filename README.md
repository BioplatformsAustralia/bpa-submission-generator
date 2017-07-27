# Bioplatforms Australia Submission Generator

Generate metadata suitable for submissions to international data
repositories, from the metadata contained in the BPA Metadata Portal.

Example usage:

bpa-submit -k <ckan-api-key> -u https://data.bioplatforms.com base-ncbi 2>&1 | tee base.log
