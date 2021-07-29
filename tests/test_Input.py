# Standard imports
from pathlib import Path
import unittest

# Third-party imports
from netCDF4 import Dataset
import numpy as np

# Local imports
from src.Input import Input

class TestInput(unittest.TestCase):
    """Tests Input class methods."""

    def test___get_gb_data(self):
        """Tests __get_gb_data method."""

        basin_json = Path(__file__).parent / "test_data" / "reaches.json"
        input = Input(None, basin_json, None)

        gb_file = Path(__file__).parent / "test_data" / "gb_data_test.nc"
        gb = Dataset(gb_file, 'r')
        actual = input._Input__get_gb_data(gb, "logQ", True)
        gb.close()

        expected = np.exp(np.array([7.052769661, np.nan, 7.637092749, 7.037340959, np.nan]))
        np.testing.assert_array_almost_equal(actual, expected, decimal=2)