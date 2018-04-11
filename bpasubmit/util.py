import logging
import json
import os

import ckanapi
import requests
from dateutil.relativedelta import relativedelta
import datetime


def make_logger(name):
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)
    handler = logging.StreamHandler()
    fmt = logging.Formatter("%(asctime)s [%(levelname)-7s] [%(threadName)s]  %(message)s")
    handler.setFormatter(fmt)
    logger.addHandler(handler)
    return logger


logger = make_logger(__name__)


def make_ckan_api(args):
    ckan = ckanapi.RemoteCKAN(args.ckan_url, apikey=args.api_key)
    return ckan


CKAN_AUTH = {
    'login': 'CKAN_USERNAME',
    'password': 'CKAN_PASSWORD'
}


# http://stackoverflow.com/questions/38271351/download-resources-from-private-ckan-datasets
def authenticated_ckan_session(ckan):
    s = requests.Session()
    data = dict((k, os.environ.get(v)) for k, v in CKAN_AUTH.items())
    if any(t is None for t in data.values()):
        raise Exception('please set %s' % (', '.join(CKAN_AUTH.values())))
    url = ckan.address + '/login_generic'
    r = s.post(url, data=data)
    if 'field-login' in r.text:
        raise RuntimeError('Login failed.')
    return s


def ckan_packages_of_type(ckan, typ, limit=10000):
    # 10,000 is hard-coded in the BPA version of CKAN (upped from default limit of 1,000
    # cache for local dev
    cache_filename = 'cache/{}.json'.format(typ)
    try:
        with open(cache_filename) as fd:
            return json.load(fd)
    except IOError:
        data = ckan.action.package_search(q='type:%s' % typ, include_private=True, rows=limit)['results']
        with open(cache_filename, 'w') as fd:
            json.dump(data, fd, indent=2, sort_keys=True)
        return data


def common_values(dicts):
    """
    given a list of dicts, return a dict with only the values shared
    in common between those dicts
    """
    all_keys = set()
    for d in dicts:
        all_keys = all_keys.union(set(d.keys()))
    r = {}
    for k in all_keys:
        vals = set([d.get(k) for d in dicts])
        if len(vals) == 1:
            r[k] = dicts[0][k]
    return r


def ckan_spatial_to_ncbi_lat_lon(obj, default=''):
    spatial_json = obj.get('spatial')
    if not spatial_json:
        return default
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


def bpa_id_slash(bpa_id, default=None):
    """
    replace the last '.' in bpa_id with a '/'
    """
    if not bpa_id:
        return default
    return '/'.join(bpa_id.rsplit('.', 1))


def bpa_id_short(bpa_id, default=None):
    """
    short version of a bpa_id, the number after the last '.'
    """
    if not bpa_id:
        return default
    return bpa_id.split('.')[-1]


def apply_embargo(ckan_packages, months):
    def within_embargo(package):
        ingest_date = package.get('archive_ingestion_date')
        if not ingest_date:
            logger.error('Skipping package (no archive_ingestion_date) bpa_id: {} id: {} archive_ingestion_date: {} has-resources: {} ticket: {}'.format(package.get('bpa_id'), package.get('id'), package.get('archive_ingestion_date'), 'resources' in package, package.get('ticket')))
            return False

        ingest_date = datetime.datetime.strptime(ingest_date, "%Y-%m-%d").date()
        if ingest_date + relativedelta(months=3) > datetime.date.today():
            logger.error('Skipping package (embargoed) bpa_id: {0} id: {1} archive_ingestion_date: {2} has-resources: {3}'.format(package.get('bpa_id'), package.get('id'), package.get('archive_ingestion_date'), 'resources' in package))
            return False

        return True

    return [t for t in ckan_packages if within_embargo(t)]


def _build_hiseq_fix_map():
    m = {}
    _hiseq_prefixes = 'HiSeq', 'HiSeq ', 'Illumina HiSeq', 'Illumina HiSeq '
    for i in range(2500, 2599):
        for p in _hiseq_prefixes:
            m["%s%d" % (p, i)] = ILLUMINA_HISEQ_DEFAULT
    return m


ILLUMINA_HISEQ_DEFAULT = 'Illumina HiSeq 2500'
_hiseq_fix_map = _build_hiseq_fix_map()


def fix_instrument_hiseq_model(obj):
    original = obj.get('sequencer', '').strip()
    if not original:
        renamed = ILLUMINA_HISEQ_DEFAULT
    else:
        renamed = _hiseq_fix_map.get(original)
    if renamed is not None and renamed != original:
        logger.warn('Rename (instrument_model) bpa_id: {} id: {} ({} -> {})'.format(obj.get('bpa_id'), obj.get('id'), original, renamed))
        return renamed
    return original
