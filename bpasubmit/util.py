import logging
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
