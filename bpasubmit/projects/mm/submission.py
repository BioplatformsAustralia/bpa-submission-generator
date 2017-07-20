from collections import defaultdict
from ...util import make_logger, ckan_packages_of_type, common_values, ckan_spatial_to_ncbi_lat_lon
from ...ncbi.biosample import NCBIBioSampleMetagenomeEnvironmental
from ...ncbi.srasubtemplate import NCBISRASubtemplate

logger = make_logger(__name__)


class MarineMicrobes(object):
    def __init__(self, ckan, args):
        self.ckan = ckan
        self.amplicons = ckan_packages_of_type(ckan, 'mm-genomics-amplicon')
        self.metagenomics = ckan_packages_of_type(ckan, 'mm-metagenomics')
        self.metatranscriptome = ckan_packages_of_type(ckan, 'mm-metatranscriptome')
        self.packages = self.metagenomics + self.amplicons + self.metatranscriptome
        self.write_ncbi()

    @classmethod
    def _build_id_depth_metadata(cls, packages):
        # group together by (bpa_id, depth), then take the common values
        by_bpaid_depth = defaultdict(list)
        for package in packages:
            # cooerce to string, so lists and nested objects don't blow up common_values
            by_bpaid_depth[(package['bpa_id'], package.get('depth', ''))].append(
                dict((t, unicode(u)) for (t, u) in package.items()))

        uniqued = []
        for (bpa_id, depth), packages in by_bpaid_depth.items():
            uniqued.append(common_values(packages))
        return uniqued

    @classmethod
    def packages_to_submit(cls, packages):
        for obj in sorted(packages, key=lambda obj: int(obj['bpa_id'].rsplit('.', 1)[1])):
            yield obj

    def ncbi_metagenome_objects(self):

        def represent_depth(depth):
            # some are floating point values, but we need to integer-f
            try:
                return int(float(depth))
            except ValueError:
                # not cast-able to a float, just return as a string
                # example: "10_20"
                return depth

        def generate_isolate(bpa_id, depth):
            if not bpa_id or not depth:
                return ''
            return '%s_%s' % (bpa_id, represent_depth(depth))

        id_depth_metadata = self._build_id_depth_metadata(self.packages)
        for obj in self.packages_to_submit(id_depth_metadata):
            bpa_id_slash = '/'.join(obj['bpa_id'].rsplit('.', 1))
            depth = obj.get('depth', '')
            yield {
                'sample_name': bpa_id_slash,
                'collection_date': obj.get('date_sampled', '2015'),
                'geo_loc_name': '%s: %s' % (obj.get('geo_loc_name', 'Australia'), obj.get('location_description', '')),
                'lat_lon': ckan_spatial_to_ncbi_lat_lon(obj),
                'bioproject_accession': 'PRJNA385736',
                'depth': depth,
                'isolate': generate_isolate(bpa_id_slash, depth),
                # constant values: FIXME, put these in CKAN once confirmed correct
                'organism': 'marine metagenome',
                'isolation_source': 'Marine',
            }

    def ncbi_sra_objects(self):
        base_obj = {
            'bioproject_accession': 'PRJNA385736',
            'design_description': 'http://www.bioplatforms.com/marine-microbes/',
            'reference_genome_assembly': '',
            'alignment_software': '',
            'library_selection': 'PCR',
            'library_layout': 'paired',
            'platform': 'ILLUMINA',
        }
        # genomics amplicons: each row is a unique (bpa_id, amplicon, flow_cell_id): which happens
        # to be how we modelled things in CKAN
        for obj in self.packages_to_submit(self.amplicons):
            bpa_id_slash = '/'.join(obj['bpa_id'].rsplit('.', 1))
            file_info = []
            for resource_obj in obj['resources']:
                if resource_obj['read'] not in ('R1', 'R2'):
                    continue
                file_info.append(['fastq', resource_obj['url'].rsplit('/', 1)[-1], resource_obj['md5']])
            row_obj = base_obj.copy()
            row_obj.update({
                'biosample_accession': obj.get('ncbi_biosample_accession', ''),
                'sample_name': bpa_id_slash,
                'library_ID': '%s_%s_%s' % (obj['bpa_id'].split('.')[-1], obj['amplicon'].upper(), obj['flow_id']),
                'title/short description': 'Marine_amplicon',
                'library_strategy': 'AMPLICON',
                'library_source': 'GENOMIC',
                'instrument_model': 'Illumina MiSeq',
                'forward_read_length': obj['read_length'],
                'reverse_read_length': obj['read_length'],
                'filetype': 'fastq',
            })
            yield row_obj, file_info
        # metagenomics
        for obj in self.packages_to_submit(self.metagenomics):
            bpa_id_slash = '/'.join(obj['bpa_id'].rsplit('.', 1))
            file_info = []
            for resource_obj in obj['resources']:
                if resource_obj['read'] not in ('R1', 'R2'):
                    continue
                if resource_obj.get('ncbi_file_uploaded') is True:
                    logger.debug('skipped an uploaded file')
                    continue
                file_info.append(['fastq', resource_obj['url'].rsplit('/', 1)[-1], resource_obj['md5']])
            row_obj = base_obj.copy()
            row_obj.update({
                'biosample_accession': obj.get('ncbi_biosample_accession', ''),
                'sample_name': bpa_id_slash,
                'library_ID': '%s' % (obj['bpa_id'].split('.')[-1]),
                'title/short description': 'Marine_metagenomics',
                'library_strategy': 'METAGENOMICS',
                'library_source': 'METAGENOMICS',
                'instrument_model': obj.get('library_construction_protocol', ''),
                'forward_read_length': obj['read_length'],
                'reverse_read_length': obj['read_length'],
                'filetype': 'fastq',
            })
            yield row_obj, file_info
        # metatranscriptome
        for obj in self.packages_to_submit(self.metagenomics):
            bpa_id_slash = '/'.join(obj['bpa_id'].rsplit('.', 1))
            file_info = []
            for resource_obj in obj['resources']:
                if resource_obj['read'] not in ('R1', 'R2'):
                    continue
                if resource_obj.get('ncbi_file_uploaded') is True:
                    logger.debug('skipped an uploaded file')
                    continue
                file_info.append(['fastq', resource_obj['url'].rsplit('/', 1)[-1], resource_obj['md5']])
            row_obj = base_obj.copy()
            row_obj.update({
                'biosample_accession': obj.get('ncbi_biosample_accession', ''),
                'sample_name': bpa_id_slash,
                'library_ID': '%s' % (obj['bpa_id'].split('.')[-1]),
                'title/short description': 'Marine_metatranscriptome',
                'library_strategy': 'METATRANSCRIPTOME',
                'library_source': 'METATRANSCRIPTOME',
                'instrument_model': obj.get('library_construction_protocol', ''),
                'forward_read_length': obj['read_length'],
                'reverse_read_length': obj['read_length'],
                'filetype': 'fastq',
            })
            yield row_obj, file_info

    def write_ncbi(self):
        with open('Metagenome.environmental.1.0-MM.tsv', 'w') as fd:
            NCBIBioSampleMetagenomeEnvironmental.write(('depth', 'isolate'), fd, self.ncbi_metagenome_objects())
        with open('SRA_subtemplate_v2-8-MM.csv', 'w') as fd:
            NCBISRASubtemplate.write(('depth', 'isolate'), fd, self.ncbi_sra_objects())
