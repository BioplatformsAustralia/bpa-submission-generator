
import csv
import itertools
from ..util import make_logger

logger = make_logger(__name__)


class NCBIBioSampleMetagenomeEnvironmental(object):
    # TSV format template, source: https://submit.ncbi.nlm.nih.gov/biosample/template/
    ncbi_template = """\
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
"""
    ncbi_field_header = '*sample_name	sample_title	bioproject_accession	*organism	host	isolation_source	*collection_date	*geo_loc_name	*lat_lon	ref_biomaterial	rel_to_oxygen	samp_collect_device	samp_mat_process	samp_size	source_material_id	description'
    fields = (
        'sample_name',
        'sample_title',
        'bioproject_accession',
        'organism',
        'host',
        'isolation_source',
        'collection_date',
        'geo_loc_name',
        'lat_lon',
        'ref_biomaterial',
        'rel_to_oxygen',
        'samp_collect_device',
        'samp_mat_process',
        'samp_size',
        'source_material_id',
        'description',
    )

    chunk_size = 1000

    @classmethod
    def chunk_write(cls, custom_fields, base_filename, rows):
        """
        write out n files with cls.chunk_size rows
        returns a mapping from sample_name to output file number
        """
        # TODO there is probably a more pythonic way of doing this
        chunk = 0
        rows_chunk = list(itertools.islice(rows, cls.chunk_size))
        info = {}
        while rows_chunk:
            chunk += 1
            filename = '{0}-{1}.tsv'.format(base_filename, chunk)
            info.update(dict((t['sample_name'], chunk) for t in rows_chunk))
            with open(filename, 'w') as fd:
                cls.write(custom_fields, fd, rows_chunk)
            rows_chunk = list(itertools.islice(rows, cls.chunk_size))
        return info

    @classmethod
    def write(cls, custom_fields, fd, rows):
        """
        write NCBI BioSample Metagenome or Environmental; version 1.0 submission sheet to `fd`
        each row in `rows` must be a dictionary with keys corresponding to the fields member of
        this class, plus any `custom_fields` provided
        """
        # note: the NCBI template uses DOS linefeeds
        fd.write(cls.ncbi_template.replace('\n', '\r\n'))
        fd.write('\t'.join((cls.ncbi_field_header,) + custom_fields))
        fd.write('\r\n')
        writer = csv.DictWriter(fd, cls.fields + custom_fields, dialect='excel-tab')
        writer.writerows(rows)
