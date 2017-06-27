import json
from collections import defaultdict
from ...util import make_logger, ckan_packages_of_type, common_values
from ...ncbi.biosample import NCBIBioSampleMetagenomeEnvironmental

logger = make_logger(__name__)


class BASE(object):
    def __init__(self, ckan, args):
        self.ckan = ckan
        self.packages = self._get_packages(ckan)
        self.write_ncbi()

    def _get_packages(self, ckan):
        # group together by (bpa_id, depth), then take the common values
        packages = ckan_packages_of_type(ckan, 'base-metagenomics') + ckan_packages_of_type(ckan, 'base-genomics-amplicon')
        by_bpaid_depth = defaultdict(list)
        for package in packages:
            # cooerce to string, so lists and nested objects don't blow up common_values
            by_bpaid_depth[(package['bpa_id'], package.get('depth', ''))].append(
                dict((t, unicode(u)) for (t, u) in package.items()))

        uniqued = []
        for (bpa_id, depth), packages in by_bpaid_depth.items():
            uniqued.append(common_values(packages))
        return uniqued

    def ncbi_objects(self):
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

        logger.debug(self.packages[0])
        logger.debug(list(sorted(self.packages[0].keys())))
        logger.debug(self.packages[0])
        for obj in sorted(self.packages, key=lambda obj: int(obj['bpa_id'].rsplit('.', 1)[1])):
            bpa_id_slash = '/'.join(obj['bpa_id'].rsplit('.', 1))
            depth = obj.get('depth', '')
            yield {
                'sample_name': bpa_id_slash,
                'collection_date': obj.get('date_sampled', '2015'),
                'geo_loc_name': 'Australia: %s' % obj.get('location_description', ''),
                'lat_lon': ncbi_lat_lon(obj),
                'bioproject_accession': obj.get('ncbi_bioproject_accession', ''),
                'depth': depth,
                'isolate': generate_isolate(bpa_id_slash, depth),
                # constant values: FIXME, put these in CKAN once confirmed correct
                'organism': 'soil metagenome',
                'isolation_source': 'Soil',
            }

    def write_ncbi(self):
        with open('Metagenome.environmental.1.0-BASE.tsv', 'w') as fd:
            NCBIBioSampleMetagenomeEnvironmental.write(('depth', 'isolate'), fd, self.ncbi_objects())
