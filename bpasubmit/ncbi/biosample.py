
class NCBIBioSampleMetagenomeEnvironmental(object):
    # TSV format template, source: https://submit.ncbi.nlm.nih.gov/biosample/template/
    template = """
# This is a submission template for batch deposit of 'Metagenome or environmental; version 1.0' samples to the NCBI BioSample database (http://www.ncbi.nlm.nih.gov/biosample/).
# Fields with an asterisk (*) are mandatory. Your submission will fail if any mandatory fields are not completed. If information is unavailable for any mandatory field, please enter 'not collected', 'not applicable' or 'missing' as appropriate.
# All other fields are optional. Leave optional fields empty if no information is available.
# You can add any number of custom fields to fully describe your BioSamples, simply include them in the table.
# CAUTION: Be aware that Excel may automatically apply formatting to your data. In particular, take care with dates, incrementing autofills and special characters like / or -. Doublecheck that your text file is accurate before uploading to BioSample.
# TO MAKE A SUBMISSION:
#     1. Complete the template table (typically in Excel, or another spreadsheet application)
#     2. Save the worksheet as a Text (Tab-delimited) file - use 'File, Save as, Save as type: Text (Tab-delimited)'
#     3. Upload the file on the 'Attributes' tab of the BioSample Submission Portal at https://submit.ncbi.nlm.nih.gov/subs/biosample/.
#     4. If you have any questions, please contact us at biosamplehelp@ncbi.nlm.nih.gov.
*sample_name	sample_title	bioproject_accession	*organism	host	isolation_source	*collection_date	*geo_loc_name	*lat_lon	ref_biomaterial	rel_to_oxygen	samp_collect_device	samp_mat_process	samp_size	source_material_id	description
"""
    pass

