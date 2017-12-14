
import csv
from ..util import make_logger

logger = make_logger(__name__)


class NCBISRASubtemplate(object):
    fields = (
        'bioproject_accession',
        'biosample_accession',
        'sample_name',
        'library_ID',
        'title',
        'library_strategy',
        'library_source',
        'library_selection',
        'library_layout',
        'platform',
        'instrument_model',
        'design_description',
        'reference_genome_assembly',
        'alignment_software',
        'forward_read_length',
        'reverse_read_length')
    file_header = (
        'filetype',
        'filename',
        'MD5_checksum')

    chunk_size = 500

    @classmethod
    def numbered_file_header(cls, count):
        """
        return a tuple of count file headers
        """
        rval = ('filetype', 'filename', 'MD5_checksum')
        for i in range(2, count + 1):
            rval = rval + ('filetype' + str(i), 'filename' + str(i), 'MD5_checksum' + str(i))
        return rval

    @classmethod
    def write(cls, custom_fields, fd, rows):
        """
        write NCBI SRA Subtemplate v2.8
        """
        # note: the NCBI template uses DOS linefeeds
        writer = csv.writer(fd)
        # writer.writerow(cls.fields + cls.file_header * 4)
        writer.writerow(cls.fields + cls.numbered_file_header(4))
        for row_obj, file_objs in rows:
            if not file_objs:
                continue
            row = [row_obj[t] for t in cls.fields]
            for file_obj in sorted(file_objs):
                row += file_obj
            writer.writerow(row)
