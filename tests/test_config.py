from pyhdx.fileIO import csv_to_protein, txt_to_np, read_dynamx, fmt_export, csv_to_np
from pyhdx.models import Protein
from pathlib import Path
from pyhdx.panel.configurations import ConfigurationSettings
from io import StringIO
import numpy as np
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



