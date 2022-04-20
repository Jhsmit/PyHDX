import configparser
from pathlib import Path
from pyhdx._version import get_versions
from packaging import version
import torch
import warnings


__version__ = get_versions()["version"]
del get_versions


def read_config(path):
    """read .ini config file at path, return configparser.ConfigParser object"""
    config = configparser.ConfigParser()
    config.read(path)

    return config


def write_config(path, config):
    """write a config .ini file from a configparse.Configparser object"""
    with open(path, "w") as config_file:
        version_string = "; pyhdx configuration file " + __version__ + "\n\n"
        config_file.write(version_string)
        config.write(config_file)


def reset_config():
    """create a new config.ini file in the user home dir/.pyhdx folder"""

    with open(config_file_path, "w") as target:
        version_string = "; pyhdx configuration file " + __version__ + "\n\n"
        target.write(version_string)

        with open(current_dir / "config.ini") as source:
            for line in source:
                target.write(line)

    # shutil.copy(current_dir / 'config.ini', config_file_path)


class Singleton(type):
    # https://stackoverflow.com/questions/6760685/creating-a-singleton-in-python
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
        if pth.exists():
            self._config = read_config(pth)
        else:
            self._config = read_config(current_dir / "config.ini")

    def load_config(self, pth):
        """load a new configuration from pth"""
        self._config = read_config(pth)

    def __getitem__(self, item):
        return self._config.__getitem__(item)

    def get(self, *args, **kwargs):
        """configparser get"""
        return self._config.get(*args, **kwargs)

    def getint(self, *args, **kwargs):
        return self._config.getint(*args, **kwargs)

    def getfloat(self, *args, **kwargs):
        return self._config.getfloat(*args, **kwargs)

    def getboolean(self, *args, **kwargs):
        return self._config.getboolean(*args, **kwargs)

    def set(self, *args, **kwargs):
        """configparser set"""
        self._config.set(*args, **kwargs)

    def write_config(self, path=None):
        """
        This method is used to update the configuration file.
        """

        warnings.warn(
            "write_config method is deprecation, use the module level function instead",
            DeprecationWarning,
        )

        pth = path or config_file_path

        with open(pth, "w") as config_file:
            self._config.write(config_file)

    @property
    def assets_dir(self):
        spec_path = self.get("server", "assets_dir")
        assets_dir = Path(spec_path.replace("~", str(home_dir)))

        return assets_dir

    @property
    def TORCH_DTYPE(self):
        dtype = self.get("fitting", "dtype")
        if dtype in ["float64", "double"]:
            return torch.float64
        elif dtype in ["float32", "float"]:
            return torch.float32
        else:
            raise ValueError(f"Unsupported data type: {dtype}")

    @property
    def TORCH_DEVICE(self):
        device = self.get("fitting", "device")
        return torch.device(device)


def valid_config():
    """Checks if the current config file in the user home directory is a valid config
    file for the current pyhdx version

    """
    if not config_file_path.exists():
        return False
    else:
        with open(config_file_path, "r") as f:
            version_string = f.readline().strip("; ").split(" ")[-1]

        pyhdx_version = version.parse(__version__)
        cfg_version = version.parse(version_string)

        return pyhdx_version.public == cfg_version.public


home_dir = Path.home()
config_dir = home_dir / ".pyhdx"
config_dir.mkdir(parents=False, exist_ok=True)
current_dir = Path(__file__).parent

config_file_path = config_dir / "config.ini"
if not valid_config():
    try:
        reset_config()
    except FileNotFoundError:
        # This will happen on conda-forge docker build.
        # When no config.ini file is in home_dir / '.pyhdx',
        # ConfigurationSettings will use the hardcoded version (pyhdx/config.ini)
        pass
        # ( this is run twice due to import but should be OK since cfg is singleton)

cfg = ConfigurationSettings()
