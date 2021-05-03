from pathlib import Path
from pyhdx.panel.config import ConfigurationSettings
import pytest

directory = Path(__file__).parent


class TestConfig(object):

    def test_update_cluster(self):
        cluster = ('127.0.0.1', 52125)
        with pytest.raises(TypeError):
            ConfigurationSettings().update_cluster(*cluster)

        cluster = ('127.0.0.1', '52125')
        ConfigurationSettings().update_cluster(*cluster)

        set_value = ConfigurationSettings().load_cluster()
        assert set_value == ':'.join(cluster)



