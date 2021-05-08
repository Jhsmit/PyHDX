from pathlib import Path
from pyhdx.panel.config import ConfigurationSettings, read_config, config_file_path, reset_config
import pytest

directory = Path(__file__).parent


class TestConfig(object):
    def test_cfg_singleton(self):
        cfg = ConfigurationSettings()
        cluster = '127.0.0.1:00000'
        cfg.cluster = cluster

        assert cfg.get('cluster', 'port') == '00000'

        cfg2 = ConfigurationSettings()
        assert id(cfg) == id(cfg2)

        cluster = 'asdfsadf'
        with pytest.raises(ValueError):
            cfg.cluster = cluster

        cfg.write_config()

        cp_config = read_config(config_file_path)
        assert cp_config['cluster']['port'] == '00000'

        reset_config()
        cp_config = read_config(config_file_path)
        assert cp_config['cluster']['port'] == '52123'



