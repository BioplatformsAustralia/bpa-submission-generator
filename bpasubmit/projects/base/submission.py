import json
from collections import defaultdict
from ...util import make_logger, ckan_packages_of_type, common_values
from ...ncbi.biosample import NCBIBioSampleMetagenomeEnvironmental
from ...ncbi.srasubtemplate import NCBISRASubtemplate

logger = make_logger(__name__)


class BASE(object):
    def __init__(self, ckan, args):
        self.ckan = ckan
        self.metagenomics = ckan_packages_of_type(ckan, 'base-metagenomics')
        self.amplicons = ckan_packages_of_type(ckan, 'base-genomics-amplicon')
        self.packages = self.metagenomics + self.amplicons
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

    def ncbi_metagenome_objects(self):
        def ncbi_lat_lon(obj):
            spatial_json = obj.get('spatial')
            if not spatial_json:
                return ''
            spatial = json.loads(spatial_json)
            lng, lat = spatial['coordinates']
            n_s = 'N'
            if lat < 0:
                lat = abs(lat)
                n_s = 'S'
            e_w = 'E'
            if lng < 0:
                lng = abs(lng)
                e_w = 'W'
            return '%f %s %f %s' % (lat, n_s, lng, e_w)

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
        for obj in sorted(id_depth_metadata, key=lambda obj: int(obj['bpa_id'].rsplit('.', 1)[1])):
            bpa_id_slash = '/'.join(obj['bpa_id'].rsplit('.', 1))
            depth = obj.get('depth', '')
            yield {
                'sample_name': bpa_id_slash,
                'collection_date': obj.get('date_sampled', '2015'),
                'geo_loc_name': '%s: %s' % (obj.get('geo_loc_name', 'Australia'), obj.get('location_description', '')),
                'lat_lon': ncbi_lat_lon(obj),
                'bioproject_accession': obj.get('ncbi_bioproject_accession', ''),
                'depth': depth,
                'isolate': generate_isolate(bpa_id_slash, depth),
                # constant values: FIXME, put these in CKAN once confirmed correct
                'organism': 'soil metagenome',
                'isolation_source': 'Soil',
            }

    def ncbi_sra_objects(self):
        # genomics amplicons: each row is a unique (bpa_id, amplicon, flow_cell_id): which happens
        # to be how we modelled things in CKAN
        for obj in sorted(self.amplicons, key=lambda obj: int(obj['bpa_id'].rsplit('.', 1)[1])):
            bpa_id_slash = '/'.join(obj['bpa_id'].rsplit('.', 1))
            file_info = []
            for resource_obj in obj['resources']:
                if resource_obj['read'] not in ('R1', 'R2'):
                    continue
                file_info.append(['fastq', resource_obj['url'].rsplit('/', 1)[-1], resource_obj['md5']])
            yield ({
                'bioproject_accession': obj.get('ncbi_bioproject_accession', ''),
                'biosample_accession': obj.get('ncbi_biosample_accession', ''),
                'sample_name': bpa_id_slash,
                'library_ID': '%s_%s_%s' % (obj['bpa_id'].split('.')[-1], obj['amplicon'].upper(), obj['flow_id']),
                'title/short description': 'Soil_amplicon',
                'library_strategy': 'AMPLICON',
                'library_source': 'GENOMIC',
                'library_selection': 'PCR',
                'library_layout': 'paired',
                'platform': 'ILLUMINA',
                'instrument_model': 'Illumina MiSeq',
                'design_description': 'http://www.bioplatforms.com/soil-biodiversity/',
                'reference_genome_assembly': '',
                'alignment_software': '',
                'forward_read_length': '',  # FIXME
                'reverse_read_length': '',  # FIXME
                'filetype': 'fastq',
            }, file_info)
        # metagenomics
        for obj in sorted(self.metagenomics, key=lambda obj: int(obj['bpa_id'].rsplit('.', 1)[1])):
            bpa_id_slash = '/'.join(obj['bpa_id'].rsplit('.', 1))
            file_info = []
            for resource_obj in obj['resources']:
                if resource_obj['read'] not in ('R1', 'R2'):
                    continue
                file_info.append(['fastq', resource_obj['url'].rsplit('/', 1)[-1], resource_obj['md5']])
            yield ({
                'bioproject_accession': obj.get('ncbi_bioproject_accession', ''),
                'biosample_accession': obj.get('ncbi_biosample_accession', ''),
                'sample_name': bpa_id_slash,
                'library_ID': '%s_%s' % (obj['bpa_id'].split('.')[-1], obj['flow_id']),
                'title/short description': 'Soil_metagenomics',
                'library_strategy': 'METAGENOMICS',
                'library_source': 'METAGENOMICS',
                'library_selection': 'PCR',
                'library_layout': 'paired',
                'platform': 'ILLUMINA',
                'instrument_model': obj.get('library_construction_protocol', ''),
                'design_description': 'http://www.bioplatforms.com/soil-biodiversity/',
                'reference_genome_assembly': '',
                'alignment_software': '',
                'forward_read_length': '',  # FIXME
                'reverse_read_length': '',  # FIXME
                'filetype': 'fastq',
            }, file_info)

    def write_ncbi(self):
        with open('Metagenome.environmental.1.0-BASE.tsv', 'w') as fd:
            NCBIBioSampleMetagenomeEnvironmental.write(('depth', 'isolate'), fd, self.ncbi_metagenome_objects())
        with open('SRA_subtemplate_v2-8-BASE.csv', 'w') as fd:
            NCBISRASubtemplate.write(('depth', 'isolate'), fd, self.ncbi_sra_objects())
