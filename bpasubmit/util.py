import logging
import json
import os

import ckanapi
import requests


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
    return ckan.action.package_search(q='type:%s' % typ, include_private=True, rows=limit)['results']


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
