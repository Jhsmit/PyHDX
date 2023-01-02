from __future__ import annotations

from os import PathLike
from pathlib import Path
from typing import Union, Dict, Any

import torch
from omegaconf import OmegaConf, DictConfig, DictKeyType
from packaging import version

from pyhdx._version import get_versions

__version__ = get_versions()["version"]
del get_versions


def reset_config():
    """create a new config.yaml file in the user home dir/.pyhdx folder"""

    with open(config_file_path, "w") as target:
        version_string = "# pyhdx configuration file " + __version__ + "\n\n"
        target.write(version_string)

        with open(current_dir / "config.yaml") as source:
            for line in source:
                target.write(line)

    # shutil.copy(current_dir / 'config.ini', config_file_path)


class Singleton(type):
    _instances: Dict[type, "Singleton"] = {}

    def __call__(cls, *args: Any, **kwargs: Any) -> Any:
        if cls not in cls._instances:
            cls._instances[cls] = super(Singleton, cls).__call__(*args, **kwargs)
        return cls._instances[cls]

    def instance(cls: Any, *args: Any, **kwargs: Any) -> Any:
        return cls(*args, **kwargs)


class PyHDXConfig(metaclass=Singleton):
    __slots__ = ["conf"]

    def __init__(self) -> None:
        self.conf = None

    def __getattr__(self, item: str) -> Any:
        return getattr(self.conf, item)

    def load_config(self, config_file: PathLike[str]):
        conf = OmegaConf.create(Path(config_file).read_text())
        self.set_config(conf)

    def set_config(self, conf: DictConfig) -> None:
        self.conf = conf

    def get(self, key: DictKeyType, default_value: Any = None) -> Any:
        return self.conf.get(key, default_value)

    @property
    def assets_dir(self) -> Path:
        spec_path = self.conf.server.assets_dir
        assets_dir = Path(spec_path.replace("~", str(Path().home())))

        return assets_dir

    @property
    def log_dir(self) -> Path:
        spec_path = self.conf.server.log_dir
        log_dir = Path(spec_path.replace("~", str(Path().home())))

        return log_dir

    @property
    def TORCH_DTYPE(self) -> Union[torch.float64, torch.float32]:
        dtype = self.conf.fitting.dtype
        if dtype in ["float64", "double"]:
            return torch.float64
        elif dtype in ["float32", "float"]:
            return torch.float32
        else:
            raise ValueError(f"Unsupported data type: {dtype}")

    @property
    def TORCH_DEVICE(self) -> torch.device:
        device = self.conf.fitting.device
        return torch.device(device)


def valid_config() -> bool:
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
config_file_path = config_dir / "config.yaml"

current_dir = Path(__file__).parent

# Current config version is outdated
if not valid_config():
    try:
        reset_config()
        conf = OmegaConf.load(config_file_path)
    except FileNotFoundError:
        # This will happen on conda-forge docker build.
        # When no config.yaml file is in home_dir / '.pyhdx',
        # ConfigurationSettings will use the hardcoded version (pyhdx/config.ini)
        conf = OmegaConf.load(current_dir / "config.yaml")
        # (this is run twice due to import but should be OK since conf is singleton)
else:
    conf = OmegaConf.load(config_file_path)


cfg = PyHDXConfig()
cfg.set_config(conf)
