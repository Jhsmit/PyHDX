from pathlib import Path
from pyhdx.config import (
    ConfigurationSettings,
    read_config,
    config_file_path,
    reset_config,
    write_config,
)
import pytest

directory = Path(__file__).parent


class TestConfig(object):
    def test_cfg_singleton(self, tmp_path):
        cfg = ConfigurationSettings()
        scheduler_address = "127.0.0.1:00000"
        cfg.set("cluster", "scheduler_address", scheduler_address)

        assert cfg.get("cluster", "scheduler_address") == scheduler_address

        cfg2 = ConfigurationSettings()
        assert id(cfg) == id(cfg2)

        fpath = Path(tmp_path) / "config.ini"
        write_config(fpath, cfg._config)

        cp_config = read_config(fpath)
        assert cp_config["cluster"]["scheduler_address"] == "127.0.0.1:00000"

        reset_config()  # the effect of this is currently not tested, needs more extensive tests
        cp_config = read_config(config_file_path)
        assert cp_config["cluster"]["scheduler_address"] == "127.0.0.1:52123"
