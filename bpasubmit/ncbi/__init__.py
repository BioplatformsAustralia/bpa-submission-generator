
from collections import Counter, defaultdict
from ..util import bpa_id_short
from .srasubtemplate import NCBISRASubtemplate
from .biosample import NCBIBioSampleMetagenomeEnvironmental


def write_sra_biosample(biosample_custom_fields, biosample_base, biosample_rows, sra_custom_fields, sra_base, sra_rows):
    #
    # write out the BioSample and SRA submission files, with a one to one link between each BioSample
    # file and a corresponding SRA file. In practice this means that BioSample files will tend to be
    # short.
    #

    sample_srarows = Counter(row['sample_name'] for row, _ in sra_rows if row['sample_name'])

    # bin samples into SRA template files
    sra_chunks = []
    chunk = []
    counter = 0
    for sample_id, srarows in sorted(sample_srarows.items(), key=lambda kv: int(kv[0].split('/', 1)[-1])):
        if counter + srarows > NCBISRASubtemplate.chunk_size:
            sra_chunks.append(chunk)
            counter = 0
            chunk = []
        chunk.append(sample_id)
        counter += srarows
    
    print(len(sra_chunks))
    print(sra_chunks)

