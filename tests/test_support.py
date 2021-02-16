import numpy as np
import matplotlib as mpl
from pyhdx.support import rgb_to_hex


class TestSupportFunctions(object):

    def test_color_converions(self):
        all_rbg = np.arange(256 ** 3, dtype=np.uint32).view(np.uint8).reshape(256 ** 3, -1)

        selected_rgb = all_rbg[np.random.randint(0, 1000)::1000]
        hex_mpl = np.array([mpl.colors.to_hex(rgba).upper() for rgba in selected_rgb / 255])

        hex_pyhdx = rgb_to_hex(selected_rgb)
        assert np.all(hex_pyhdx == hex_mpl)
