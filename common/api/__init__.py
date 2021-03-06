import jwt
import requests
from datetime import datetime, timedelta
from common.basedir import PERSIST
from selfdrive.version import version
import base64
import json

class Api():
  def __init__(self, dongle_id):
    self.dongle_id = dongle_id
    with open(PERSIST+'/comma/id_rsa') as f:
      self.private_key = f.read()

  def get(self, *args, **kwargs):
    return self.request('GET', *args, **kwargs)

  def post(self, *args, **kwargs):
    return self.request('POST', *args, **kwargs)

  def request(self, method, endpoint, timeout=None, access_token=None, **params):
    return api_get(endpoint, method=method, timeout=timeout, access_token=access_token, **params)

  def get_token(self):
    now = datetime.utcnow()
    payload = {
      'identity': self.dongle_id,
      'nbf': now.isoformat(),
      'iat': now.isoformat(),
      'exp': (now + timedelta(hours=1)).isoformat()
    }

    return base64.b64encode(json.dumps(payload).encode('ascii')).decode('ascii')
    
    #return jwt.encode(payload, self.private_key, algorithm='RS256').decode('utf8')

def api_get(endpoint, method='GET', timeout=None, access_token=None, **params):
  backend = "http://openpilot.api.codemotive.io/"

  headers = {}
  if access_token is not None:
    headers['Authorization'] = "NONE "+access_token

  headers['User-Agent'] = "openpilot-" + version

  return requests.request(method, backend+endpoint, timeout=timeout, headers = headers, params=params)

