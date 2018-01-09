
from collections import Counter
from .srasubtemplate import NCBISRASubtemplate
from .biosample import NCBIBioSampleMetagenomeEnvironmental
import itertools


# from itertools docs
def grouper(iterable, n, fillvalue=None):
    args = [iter(iterable)] * n
    return itertools.zip_longest(*args, fillvalue=fillvalue)


def write_sra_biosample(biosample_custom_fields, biosample_base, biosample_rows, sra_custom_fields, sra_base, sra_rows):
    #
    # write out the BioSample and SRA submission files, with a one to one link between each BioSample
    # file and a corresponding SRA file. In practice this means that BioSample files will tend to be
    # short.
    #

    # coalesce so we can slice and dice
    biosample_rows = list(biosample_rows)
    sra_rows = list(sra_rows)

    sample_nsrarows = Counter(row['sample_name'] for row, _ in sra_rows if row['sample_name'])

    # bin samples into SRA template files where new samples are being uploaded
    sra_chunks = []
    chunk = []
    counter = 0
    for sample_id, srarows in sorted(sample_nsrarows.items(), key=lambda kv: int(kv[0].split('/', 1)[-1])):
        if counter + srarows > NCBISRASubtemplate.chunk_size:
            sra_chunks.append(chunk)
            counter = 0
            chunk = []
        chunk.append(sample_id)
        counter += srarows

    # For each chunk, write out the BioSample and SRA templates
    for output_filenum, sample_ids in enumerate(sra_chunks, start=1):
        br = [row for row in biosample_rows if row['sample_name'] in sample_ids]
        biosample_filename = '{}-{}.tsv'.format(biosample_base, output_filenum)
        with open(biosample_filename, 'w') as fd:
            NCBIBioSampleMetagenomeEnvironmental.write(biosample_custom_fields, fd, br)

        sr = [row_res for row_res in sra_rows if row_res[0]['sample_name'] in sample_ids]
        sra_filename = '{}-{}.tsv'.format(sra_base, output_filenum)
        with open(sra_filename, 'w') as fd:
            NCBISRASubtemplate.write(sra_custom_fields, fd, sr)

    # Spit out the file uploads for the existing samples
    sra_existing = [t for t in sra_rows if not t[0]['sample_name']]
    for output_filenum, sr in enumerate(grouper(sra_existing, NCBISRASubtemplate.chunk_size), start=1):
        sr = [t for t in sr if t]
        sra_filename = '{}-SA{}.tsv'.format(sra_base, output_filenum)
        with open(sra_filename, 'w') as fd:
            NCBISRASubtemplate.write(sra_custom_fields, fd, sr)
