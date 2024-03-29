from collections import defaultdict
from ...util import sample_id_short, sample_id_slash, make_logger, ckan_packages_of_type, common_values, ckan_spatial_to_ncbi_lat_lon, apply_embargo, fix_instrument_hiseq_model
from ...ncbi import write_sra_biosample

logger = make_logger(__name__)


class MarineMicrobes(object):
    def __init__(self, ckan, args):
        self.ckan = ckan

        def with_embargo(typ):
            return apply_embargo(ckan_packages_of_type(ckan, typ), months=3)

        mandatory_fields = ('utc_date_sampled', 'geo_loc_name', 'spatial')

        def with_mandatory(packages):
            r = []
            for obj in packages:
                missing = [t for t in mandatory_fields if t not in obj]
                if missing:
                    logger.warn('Skipping package (missing_mandatory {}) package_id: {}'.format(
                        str(mandatory_fields), obj.get('id')))
                    continue
                r.append(obj)
            return r

        self.amplicons = with_mandatory(with_embargo('mm-genomics-amplicon'))
        self.metagenomics = with_mandatory(with_embargo('mm-metagenomics'))
        self.metatranscriptome = with_mandatory(
            with_embargo('mm-metatranscriptome'))
        self.packages = self.metagenomics + self.amplicons + self.metatranscriptome
        self.write_ncbi()

    @classmethod
    def _build_id_depth_metadata(cls, packages):
        # group together by (sample_id, depth), then take the common values
        by_bpaid_depth = defaultdict(list)
        for package in packages:
            # cooerce to string, so lists and nested objects don't blow up common_values
            by_bpaid_depth[(package['sample_id'], package.get('depth', ''))].append(
                dict((t, str(u)) for (t, u) in list(package.items())))

        uniqued = []
        for (sample_id, depth), packages in list(by_bpaid_depth.items()):
            uniqued.append(common_values(packages))
        return uniqued

    @classmethod
    def packages_to_submit(cls, packages):
        for obj in sorted(packages, key=lambda obj: int(sample_id_short(obj['sample_id']))):
            # TODO hardcoded filter
            if not obj.get('sample_type'):
                logger.warn('Skipping package (sample_type) sample_id: {0} id: {1} sample_type: {2} has-resources: {3}'.format(
                    obj.get('sample_id'), obj.get('id'), obj.get('sample_type'), 'resources' in obj))
                continue
            yield obj

    @classmethod
    def resources_to_submit(cls, resources):
        for resource_obj in resources:

            # TODO hard coded filter on ncbi_file_uploaded
            if resource_obj.get('ncbi_file_uploaded') == 'True':
                logger.info('Skipping resource (ncbi_file_uploaded) package_id: {0} id: {1}'.format(
                    resource_obj.get('package_id'), resource_obj['id']))
                continue

            if not resource_obj.get('read'):
                logger.warn('Skipping resource (read missing) package_id: {0} id: {1} read: {2}'.format(
                    resource_obj.get('package_id'), resource_obj['id'], resource_obj.get('read')))
                continue

            # TODO hardcoded filter on read
            if resource_obj.get('read') not in ('R1', 'R2', 'I1', 'I2'):
                logger.info('Skipping resource (read) package_id: {0} id: {1} read: {2}'.format(
                    resource_obj.get('package_id'), resource_obj['id'], resource_obj.get('read')))
                continue
            yield resource_obj

    def ncbi_metagenome_objects(self):

        def represent_depth(depth):
            # some are floating point values, but we need to integer-f
            try:
                return int(float(depth))
            except ValueError:
                # not cast-able to a float, just return as a string
                # example: "10_20"
                return depth

        def generate_isolate(sample_id, depth):
            if not sample_id or not depth:
                return ''
            return '%s_%s' % (sample_id_slash(sample_id), represent_depth(depth))

        id_depth_metadata = self._build_id_depth_metadata(self.packages)
        for obj in self.packages_to_submit(id_depth_metadata):

            # Request NOT to include biosample entries where a biosample_accession already exists
            biosample_accession = obj.get('ncbi_biosample_accession', '')
            if biosample_accession:
                # logger.info('Skipping (ncbi_biosample_accession) package_id: {0} id: {1} biosample_accession: {2}'.format(obj.get('package_id'), obj['id'], biosample_accession))
                logger.info('Skipping package (ncbi_biosample_accession) package_id: {0} biosample_accession: {1}'.format(
                    obj.get('id'), biosample_accession))
                continue

            yield {
                'sample_name': sample_id_slash(obj['sample_id'], 'MANDATORY'),
                'collection_date': obj.get('utc_date_sampled', 'MANDATORY'),
                'geo_loc_name': obj.get('geo_loc_name', 'MANDATORY'),
                'lat_lon': ckan_spatial_to_ncbi_lat_lon(obj, 'MANDATORY'),
                # TODO hard coded bioproject_accession
                'bioproject_accession': 'PRJNA385736',
                'depth': obj.get('depth', ''),
                'isolate': generate_isolate(obj['sample_id'], obj.get('depth', '')),
                # TODO hard coded values: FIXME, put these in CKAN once confirmed correct
                'organism': 'marine metagenome',
                'isolation_source': obj.get('sample_type', ''),
            }

    def ncbi_sra_objects(self):

        def amplicon_specific(obj):
            # genomics amplicons: each row is a unique (sample_id, amplicon, flow_cell_id): which happens
            # to be how we modelled things in CKAN
            return {
                'library_ID': '%s_%s_%s' % (sample_id_short(obj['sample_id']), obj['amplicon'].upper(), obj['mm_amplicon_linkage']),
                # TODO hard coded values
                'title': 'Marine_amplicon',
                'library_strategy': 'AMPLICON',
                'library_source': 'GENOMIC',
                'instrument_model': 'Illumina MiSeq',
            }

        def metagenomic_specific(obj):
            # TODO hard coded instrument model field. The code is slightly redundant to allow for us to log the
            # specific issues with the data
            instrument_model = fix_instrument_hiseq_model(obj)

            return {
                'library_ID': sample_id_slash(obj['sample_id']),
                # TODO hard coded values
                'title': 'Marine_metagenomics',
                'library_strategy': 'WGS',
                'library_source': 'METAGENOMIC',
                'instrument_model': instrument_model,
            }

        def metatranscriptome_specific(obj):
            # TODO hard coded instrument model field. The code is slightly redundant to allow for us to log the
            # specific issues with the data
            instrument_model = fix_instrument_hiseq_model(obj)
            return {
                'library_ID': sample_id_slash(obj['sample_id']),
                # TODO hard coded values
                'title': 'Marine_metatranscriptome',
                'library_strategy': 'RNA-Seq',
                'library_source': 'METATRANSCRIPTOMIC',
                'instrument_model': instrument_model,
            }

        def resource_file_info(resources):
            rval = []
            for resource_obj in resources:
                # TODO hard coded values: resource_obj['Format']?
                rval.append(['fastq', resource_obj['url'].rsplit(
                    '/', 1)[-1], resource_obj['md5']])
            return rval

        # TODO hard coded values
        base_obj = {
            'bioproject_accession': 'PRJNA385736',
            'design_description': 'http://www.bioplatforms.com/marine-microbes/',
            'reference_genome_assembly': '',
            'alignment_software': '',
            'library_selection': 'PCR',
            'library_layout': 'paired',
            'platform': 'ILLUMINA',
            'filetype': 'fastq',
        }

        for obj in self.packages_to_submit(self.packages):
            file_info = resource_file_info(
                self.resources_to_submit(obj['resources']))
            row_obj = base_obj.copy()

            # biosample_accession and sample_name cannot both be set
            biosample_accession = obj.get('ncbi_biosample_accession', '')
            sample_name = sample_id_slash(obj['sample_id'])
            if biosample_accession:
                sample_name = None

            row_obj.update({
                'biosample_accession': biosample_accession,
                'sample_name': sample_name,
                # TODO received feedback that the values coming from CKAN are not always correct, so we setting them in the _specific methods
                # 'instrument_model': obj.get('sequencer', ''),
                'forward_read_length': obj['read_length'],
                'reverse_read_length': obj['read_length'],
            })

            if obj['type'] == 'mm-genomics-amplicon':
                row_obj.update(amplicon_specific(obj))
            elif obj['type'] == 'mm-metagenomics':
                row_obj.update(metagenomic_specific(obj))
            elif obj['type'] == 'mm-metatranscriptome':
                row_obj.update(metatranscriptome_specific(obj))
            else:
                logger.error('Skipping package (type) sample_id: {0} id: {1} has-resources: {2}'.format(
                    obj.get('sample_id'), obj. get('id'), 'resources' in obj))
                continue

            # TODO do we need to yield if there is no file info???
            if file_info:
                yield row_obj, file_info

    def write_ncbi(self):
        write_sra_biosample(
            biosample_custom_fields=('depth', 'isolate'),
            biosample_base='Metagenome.environmental.1.0-MM',
            biosample_rows=self.ncbi_metagenome_objects(),
            sra_custom_fields=('depth', 'isolate'),
            sra_base='SRA_subtemplate_v2-8-MM',
            sra_rows=self.ncbi_sra_objects())
