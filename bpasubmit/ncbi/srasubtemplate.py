
import csv
import itertools
from ..util import make_logger

logger = make_logger(__name__)


class NCBISRASubtemplate(object):
    fields = (
        'bioproject_accession',
        'biosample_accession',
        'sample_name',
        'library_ID',
        'title/short description',
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
    def chunk_write(cls, custom_fields, base_filename, rows):
        """
        write out n files with cls.chunk_size rows
        """
        # TODO there is probably a more pythonic way of doing this
        chunk = 0
        rows_chunk = list(itertools.islice(rows, cls.chunk_size))
        while rows_chunk:
            logger.info(len(rows_chunk))
            chunk += 1
            with open('{0}-{1}.csv'.format(base_filename, chunk), 'w') as fd:
                cls.write(custom_fields, fd, rows_chunk)
            rows_chunk = list(itertools.islice(rows, cls.chunk_size))

    @classmethod
    def write(cls, custom_fields, fd, rows):
        """
        write NCBI SRA Subtemplate v2.8
        """
        # note: the NCBI template uses DOS linefeeds
        writer = csv.writer(fd)
        writer.writerow(cls.fields + cls.file_header * 4)
        for row_obj, file_objs in rows:
            if not file_objs:
                continue
            row = [row_obj[t] for t in cls.fields]
            for file_obj in sorted(file_objs):
                row += file_obj
            writer.writerow(row)
