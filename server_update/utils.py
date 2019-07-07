import os
import gzip
import json
import shutil
import requests

from abc import ABC, abstractmethod
from datetime import datetime
from collections import OrderedDict
from typing import List
from requests.models import Response


class RemoteResource(ABC):

    info_file_schema = {
        'domain': None,
        'filename': None,
        'source': None,
        'title': None,
        'tags': [],
        'size': None,
        'datetime': None,
        # used only if files are compressed
        'uncompressed': None,
        'compression': None,
    }

    @abstractmethod
    def handle_response(self, response: Response) -> str:
        pass

    @abstractmethod
    def prepare_data(self):
        pass

    def fetch_url(self, url: str) -> Response:
        # TODO: properly handle exceptions
        print(f'Fetching data from {url}.')

        try:
            response = requests.get(url, timeout=10)
            response.raise_for_status()
        except requests.exceptions.HTTPError as err:
            raise err
        except requests.exceptions.RequestException as err:
            raise err
        else:
            return response

    @staticmethod
    def last_modified(url: str) -> datetime:
        """ Check Last-Modified header tag """
        response = requests.head(url)
        last_modified = response.headers.get('Last-Modified')
        return datetime.strptime(last_modified, '%a, %d %b %Y %H:%M:%S GMT')

    @abstractmethod
    def __init__(self):
        self.domain: str = ''
        self.file_name: str = ''
        self.title: str = ''
        self.tags: List[str] = []
        self.compression: str = 'gz'
        self.download_url: str = ''

    def create_info_file(self, **kwargs):
        info_dict = OrderedDict(self.info_file_schema)

        info_dict.update(**kwargs)
        info_dict['datetime'] = '{0:%Y-%m-%d %H:%M:%S.%f}'.format(datetime.today())
        info_dict['size'] = os.stat(self.file_name).st_size
        info_dict['source'] = 'server_file'

        with open(self.file_name + '.info', 'wt') as f:
            json.dump(info_dict, f)

    def to_serverfile_format(self):
        uncompresed_size = os.stat(self.file_name).st_size

        with open(self.file_name, 'rb') as f_in:
            with gzip.open(self.file_name + '.temp', 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

        os.remove(self.file_name)
        os.rename(self.file_name + '.temp', self.file_name)

        self.create_info_file(domain=self.domain,
                              filename=self.file_name,
                              title=self.title,
                              tags=self.tags,
                              uncompressed=uncompresed_size,
                              compression=self.compression)
