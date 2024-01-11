from pathlib import Path
import yaml
from pyhdx.config import (
    PyHDXConfig,
)
from omegaconf import OmegaConf

directory = Path(__file__).parent


class TestConfig(object):
    def test_cfg_singleton(self, tmp_path):
        cfg = PyHDXConfig()
        scheduler_address = "127.0.0.1:00000"
        cfg.cluster.scheduler_address = scheduler_address

        assert cfg.cluster.scheduler_address == scheduler_address
        assert cfg.conf.cluster.scheduler_address == scheduler_address

        cfg2 = PyHDXConfig()
        assert id(cfg) == id(cfg2)
        assert cfg.cluster.scheduler_address == scheduler_address

        s = """
        foo: 
         - bar
         - apple
         - banana
        value: 3
        """

        Path(tmp_path / "configfile.yaml").write_text(yaml.dump(s))
        conf = OmegaConf.load(tmp_path / "configfile.yaml")

        merged = OmegaConf.merge(cfg.conf, conf)

        cfg.set_config(merged)

        assert "bar" in cfg.foo
        assert cfg.value == 3
        assert cfg2.value == 3
