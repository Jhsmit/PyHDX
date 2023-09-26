from __future__ import annotations

from contextlib import contextmanager
from os import PathLike
from pathlib import Path
from typing import Union, Dict, Any, Optional, Generator

import torch
from omegaconf import OmegaConf, DictConfig, DictKeyType
from packaging import version


PACKAGE_NAME = "pyhdx"


def reset_config() -> None:
    """create a new config.yaml file in the `~home/.pyhdx` folder"""

    with open(config_file_path, "w") as target:
        from pyhdx.__version__ import __version__

        version_string = f"# {PACKAGE_NAME} configuration file " + __version__ + "\n\n"
        target.write(version_string)

        with open(current_dir / "config.yaml") as source:
            for line in source:
                target.write(line)


class Singleton(type):
    _instances: dict[type, Singleton] = {}

    def __call__(cls, *args: Any, **kwargs: Any) -> Any:
        if cls not in cls._instances:
            cls._instances[cls] = super(Singleton, cls).__call__(*args, **kwargs)
        return cls._instances[cls]

    def instance(cls: Any, *args: Any, **kwargs: Any) -> Any:
        return cls(*args, **kwargs)


class PyHDXConfig(metaclass=Singleton):
    """PyHDX configuration class.

    This object is a singleton, thus always the same instance is returned upon creation.

    Attributes:
        conf: OmegaConf DictConfig object.
    """

    __slots__ = ["conf"]

    def __init__(self) -> None:
        self.conf: Optional[DictConfig] = None

    def __getattr__(self, item: str) -> Any:
        return getattr(self.conf, item)

    def __setattr__(self, key: str, value: Any) -> None:
        if key in self.__slots__:
            super().__setattr__(key, value)
        elif key in self.conf.keys():
            setattr(self.conf, key, value)
        else:
            raise AttributeError(f"Config has no attribute {key}")

    def load_config(self, config_file: PathLike[str]) -> None:
        """Load a config file and merge with the current config.

        Args:
            config_file: Path to the config file to load.

        """
        conf_new = OmegaConf.create(Path(config_file).read_text())
        self.merge_config(conf_new)

    def merge_config(self, conf: Union[DictConfig, dict]) -> None:
        """Merge a new config with the current config.

        Args:
            conf: New config to merge with current config.

        """
        if isinstance(conf, dict):
            conf = OmegaConf.create(conf)
        conf = OmegaConf.merge(self.conf, conf)
        self.set_config(conf)

    def set_config(self, conf: DictConfig) -> None:
        self.conf = conf

    def get(self, key: DictKeyType, default_value: Any = None) -> Any:
        return self.conf.get(key, default_value)

    @property
    def assets_dir(self) -> Path:
        """PyHDX server assets directory"""
        spec_path = self.conf.server.assets_dir
        assets_dir = Path(spec_path.replace("~", str(Path().home())))

        return assets_dir

    @property
    def log_dir(self) -> Path:
        """PyHDX server log directory"""
        spec_path = self.conf.server.log_dir
        log_dir = Path(spec_path.replace("~", str(Path().home())))

        return log_dir

    @property
    def database_dir(self) -> Path:
        """HDXMS-datasets database directory"""
        spec_path = self.conf.server.database_dir
        database_dir = Path(spec_path.replace("~", str(Path().home())))

        return database_dir

    @property
    def TORCH_DTYPE(self) -> Union[torch.float64, torch.float32]:
        """PyTorch dtype used for ΔG calculations"""
        dtype = self.conf.fitting.dtype
        if dtype in ["float64", "double"]:
            return torch.float64
        elif dtype in ["float32", "float"]:
            return torch.float32
        else:
            raise ValueError(f"Unsupported data type: {dtype}")

    @property
    def TORCH_DEVICE(self) -> torch.device:
        """PyTorch device used for ΔG calculations"""
        device = self.conf.fitting.device
        return torch.device(device)

    @contextmanager
    def context(self, settings: dict) -> Generator[PyHDXConfig, None, None]:
        from pyhdx.support import rsetattr

        original_config = self.conf.copy()

        try:
            for attr, value in settings.items():
                rsetattr(cfg, attr, value)
            yield cfg
        finally:
            cfg.conf = original_config


def valid_config() -> bool:
    """Checks if the current config file in the user home directory is a valid config
    file for the current pyhdx version.

    """
    if not config_file_path.exists():
        return False
    else:
        with open(config_file_path, "r") as f:
            version_string = f.readline().strip("; ").split(" ")[-1]

        from pyhdx.__version__ import __version__

        local_version = version.parse(__version__)
        cfg_version = version.parse(version_string)

        return local_version.public == cfg_version.public


home_dir = Path.home()
config_dir = home_dir / f".{PACKAGE_NAME}"
config_dir.mkdir(parents=False, exist_ok=True)
config_file_path = config_dir / "config.yaml"

current_dir = Path(__file__).parent
conf_src_pth = current_dir / "config.yaml"


# Current config version is outdated
if not valid_config():
    try:
        reset_config()
        conf = OmegaConf.load(config_file_path)
    except FileNotFoundError:
        # This will happen on conda-forge docker build.
        # When no config.yaml file is in home_dir / '.{PACKAGE_NAME}>',
        # ConfigurationSettings will use the hardcoded version (pyhdx/config.yaml)
        conf = OmegaConf.load(conf_src_pth)
        # (this is run twice due to import but should be OK since conf is singleton)
else:
    conf = OmegaConf.load(config_file_path)


cfg = PyHDXConfig()
cfg.set_config(conf)
