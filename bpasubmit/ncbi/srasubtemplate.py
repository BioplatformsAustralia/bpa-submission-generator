
import csv


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
