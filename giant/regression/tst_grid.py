import unittest

import numpy
from scitbx.array_family import flex

from giant.grid.utils import *


class TestGridUtils(unittest.TestCase):


    def setUp(self):
        self.grid_sizes = [(1,2,3),(8,3,6),(4,1,9)]

    def test_idx_to_grid(self):
        for gs in self.grid_sizes:
            g = flex.grid(gs)
            for i, gp in enumerate(flex.nested_loop(gs)):
                gp_calc = idx_to_grid(idx=i, grid_size=gs)
                self.assertEqual(gp,gp_calc)
                self.assertEqual(i, g(gp))

