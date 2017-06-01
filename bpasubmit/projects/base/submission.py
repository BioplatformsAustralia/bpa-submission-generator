import json
from ...util import make_logger, ckan_packages_of_type
from ...ncbi.biosample import NCBIBioSampleMetagenomeEnvironmental

logger = make_logger(__name__)


class BASE(object):
    def __init__(self, ckan, args):
        self.ckan = ckan
        self.metagenome = ckan_packages_of_type(ckan, 'base-metagenomics')
        # self.amplicon = ckan_packages_of_type(ckan, 'base-genomics-amplicon')
        self.write_ncbi()

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

        def generate_isolate(bpa_id, depth):
            if not bpa_id or not depth:
                return ''
            return '%s_%d' % (bpa_id, int(float(depth)))

        logger.debug(list(sorted(self.metagenome[0].keys())))
        logger.debug(self.metagenome[0])
        for obj in sorted(self.metagenome, key=lambda obj: int(obj['bpa_id'].rsplit('.', 1)[1])):
            bpa_id_slash = '/'.join(obj['bpa_id'].rsplit('.', 1))
            depth = obj.get('depth', '')
            yield {
                'sample_name': bpa_id_slash,
                'collection_date': obj.get('date_sampled', ''),
                'geo_loc_name': obj.get('description', ''),
                'lat_lon': ncbi_lat_lon(obj),
                'depth': depth,
                'isolate': generate_isolate(bpa_id_slash, depth),
                # constant values: FIXME, put these in CKAN once confirmed correct
                'organism': 'soil metagenome',
                'isolation_source': 'Soil',
            }

    def write_ncbi(self):
        with open('Metagenome.environmental.1.0-BASE.tsv', 'w') as fd:
            NCBIBioSampleMetagenomeEnvironmental.write(('depth', 'isolate'), fd, self.ncbi_objects())
