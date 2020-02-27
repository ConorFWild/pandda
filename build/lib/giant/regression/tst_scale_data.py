import unittest

from scitbx.array_family import flex
from giant.stats.optimisation import LinearFitting

class TestLinearFitting(unittest.TestCase):

    def setUp(self):

        self.v0 = 16.0
        self.v1 = 15.0

        self.t1_x_values = flex.double(range(1000))

    def test_linear_scaling(self):

        t1_ref_values = self.v1*self.t1_x_values + self.v0
        t1_mov_values = flex.double(self.t1_x_values)

        ls = LinearFitting(x_values     = self.t1_x_values,
                           ref_y_values = t1_ref_values,
                           mov_y_values = t1_mov_values)

        self.assertAlmostEqual(round(ls.optimised_values[0],3), self.v0)
        self.assertAlmostEqual(round(ls.optimised_values[1],3), self.v1)

        # Check scaling is working correctly
        scaled = ls.transform(t1_mov_values)
        self.assertAlmostEqual(flex.sum(t1_ref_values-scaled), 0.0)

    def test_exponential_scaling(self):
        pass
