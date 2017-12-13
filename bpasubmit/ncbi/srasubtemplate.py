
import csv
import itertools
from ..util import make_logger
from collections import defaultdict

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
    def chunk_write(cls, custom_fields, base_filename, rows, biosample_chunks):
        """
        write out n files with cls.chunk_size rows.
        output files are set up to correspond with any biosample submission
        files, if relevant.
        """

        # bin rows by biosample chunk; there'll be a bin for None which is fine,
        # they just don't depend on biosample submissions at all
        bins = defaultdict(list)
        for row_obj, file_objs in rows:
            bins[biosample_chunks.get(row_obj['sample_name'], 'NA')].append((row_obj, file_objs))
        
        for biosample_filenum, bin_rows in bins.items():
            # TODO there is probably a more pythonic way of doing this
            chunk = 0
            it = iter(bin_rows)
            rows_chunk = list(itertools.islice(it, cls.chunk_size))
            while rows_chunk:
                chunk += 1
                with open('{}-{}.{}.csv'.format(base_filename, biosample_filenum, chunk), 'w') as fd:
                    cls.write(custom_fields, fd, rows_chunk)
                rows_chunk = list(itertools.islice(it, cls.chunk_size))

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
