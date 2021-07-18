import configparser
from pathlib import Path
import shutil


def read_config(path):
    config = configparser.ConfigParser()
    config.read(path)

    return config


def reset_config():
    shutil.copy(current_dir / 'config.ini', config_file_path)


class Singleton(type):
    #https://stackoverflow.com/questions/6760685/creating-a-singleton-in-python
    _instances = {}

    def __call__(cls, *args, **kwargs):
        if cls not in cls._instances:
            cls._instances[cls] = super(Singleton, cls).__call__(*args, **kwargs)
        return cls._instances[cls]


class ConfigurationSettings(metaclass=Singleton):

    def __init__(self, file_path=None):
        """
        Parameters
        ----------
        file_path : :obj:`str`  (or pathlib.Path)
            pass any *.ini file as argument to load and update, by default reads config.ini in pyhdx/panel
        """

        pth = file_path or config_file_path
        self._config = read_config(pth)

    def load_config(self, pth):
        self._config = read_config(pth)

    def __getitem__(self, item):
        return self._config.__getitem__(item)

    def get(self, *args, **kwargs):
        """configparser get"""
        return self._config.get(*args, **kwargs)

    def set(self, *args, **kwargs):
        """configparser set"""
        self._config.set(*args, **kwargs)

    @property
    def cluster(self):
        """Returns cluster address"""

        return f"{self.get('cluster', 'ip')}:{self.get('cluster', 'port')}"

    @cluster.setter
    def cluster(self, address):
        ip, port = address.split(':')
        self.set('cluster', 'ip', ip)
        self.set('cluster', 'port', port)

    def write_config(self, path=None):
        """
        This method is used to update the configuration file.
        """

        pth = path or config_file_path

        with open(pth, 'w') as config_file:
            self._config.write(config_file)


home_dir = Path.home()
config_dir = home_dir / '.pyhdx'
config_dir.mkdir(parents=False, exist_ok=True)
current_dir = Path(__file__).parent

config_file_path = config_dir / 'config.ini'
if not config_file_path.exists():
    reset_config()