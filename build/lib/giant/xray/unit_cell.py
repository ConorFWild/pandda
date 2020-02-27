
import numpy
import pandas
from math import cos, pi, sqrt
from scitbx.array_family import flex

from bamboo.common import Info


class UnitCellVariation(Info):
    def __init__(self, unit_cells):
        unit_cell_params = numpy.array([uc.parameters() for uc in unit_cells])
        self.mean  = tuple(unit_cell_params.mean(axis=0))
        self.std   = tuple(unit_cell_params.std(axis=0))
        self.max   = tuple(unit_cell_params.max(axis=0))
        self.min   = tuple(unit_cell_params.min(axis=0))
        self.range = tuple(numpy.array(self.max)-numpy.array(self.min))
        self.std_perc = tuple(100.0*numpy.array(self.std)/numpy.array(self.mean))
        self.max_perc = tuple(100.0*numpy.array(self.min)/numpy.array(self.mean))
        self.min_perc = tuple(100.0*numpy.array(self.max)/numpy.array(self.mean))
        self._initialized = True
    def as_pandas_table(self):
        """Return the variation as a pandas table"""
        import pandas
        pd = pandas.DataFrame(index=['a','b','c','alpha','beta','gamma'])
        pd['min']   = self.min
        pd['mean']  = self.mean
        pd['max']   = self.max
        pd['std']   = self.std
        pd['range'] = self.range
        pd['%std']  = self.std_perc
        pd['%min']  = self.min_perc
        pd['%max']  = self.max_perc
        return pd


#############################################################################

def unit_cell_diagonals(a, b, c, alpha=90, beta=90, gamma=90):
    """Calculate the lengths of the diagonals in the unit cell (across the three faces from the origin)"""
    d_ab = sqrt(a**2 + b**2 - 2*a*b*cos((1-gamma/180.0)*pi) )
    d_bc = sqrt(b**2 + c**2 - 2*b*c*cos((1-alpha/180.0)*pi) )
    d_ca = sqrt(c**2 + a**2 - 2*c*a*cos((1-beta/180.0)*pi)  )
    return d_ab, d_bc, d_ca

def fractional_change_for_diagonals(uc1_diags, uc2_diags):
    """Calculate the fractional change between two sets of diagonals"""
    d_ab_1, d_bc_1, d_ca_1 = uc1_diags
    d_ab_2, d_bc_2, d_ca_2 = uc2_diags
    m_ab = abs(d_ab_1-d_ab_2)/min(d_ab_1, d_ab_2)
    m_bc = abs(d_bc_1-d_bc_2)/min(d_bc_1, d_bc_2)
    m_ca = abs(d_ca_1-d_ca_2)/min(d_ca_1, d_ca_2)
    return m_ab, m_bc, m_ca

#############################################################################

def lcv_from_diagonals(uc1_diags, uc2_diags):
    chg_diags = fractional_change_for_diagonals(uc1_diags=uc1_diags, uc2_diags=uc2_diags)
    return max(chg_diags)

def lcv_from_unit_cell_parameters(uc1_params, uc2_params):
    """Calculate the Linear Cell Variation for two sets of unit cell parameters (a,b,c,alpha,beta,gamma)"""
    uc1_diags = unit_cell_diagonals(*uc1_params)
    uc2_diags = unit_cell_diagonals(*uc2.params)
    return lcv_from_diagonals(uc1_diags=uc1_diags, uc2_diags=uc2_diags)

def lcv_from_unit_cells(uc1, uc2):
    """Calculate the Linear Cell Variation for two unit cells"""
    return lcv_from_unit_cell_parameters(uc1_params=uc1.parameters(), uc2_params=uc2.parameters())

#############################################################################

def pairwise_lcv(unit_cells=None, parameters=None):
    """Calculate the LCV for all pairs of unit cells/parameters"""
    assert [unit_cells, parameters].count(None) == 1, "Provide either unit_cells or parameters"
    if unit_cells:  parameters = [uc.parameters() for uc in unit_cells]
    num_objs = len(parameters)
    diags = numpy.empty([num_objs], dtype=object)
    dist  = numpy.empty([num_objs]*2, dtype=object)
    for i_uc1, uc1 in enumerate(parameters):
        diags[i_uc1] = unit_cell_diagonals(*uc1)
    for i_uc1, i_uc2 in flex.nested_loop((num_objs, num_objs)):
        dist[i_uc1, i_uc2] = lcv_from_diagonals(uc1_diags=diags[i_uc1], uc2_diags=diags[i_uc2])
    return dist

def pairwise_ecv(unit_cells=None, parameters=None):
    """Calculate the Euclidean Cell Variation (ECV) for all pairs of unit cells/parameters"""
    assert [unit_cells, parameters].count(None) == 1, "Provide either unit_cells or parameters"
    if unit_cells:  parameters = [uc.parameters() for uc in unit_cells]
    num_objs = len(parameters)
    diags = numpy.empty([num_objs], dtype=object)
    dist  = numpy.empty([num_objs]*2, dtype=object)
    for i_uc1, uc1 in enumerate(parameters):
        diags[i_uc1] = unit_cell_diagonals(*uc1)
    for i_uc1, i_uc2 in flex.nested_loop((num_objs, num_objs)):
        dist[i_uc1, i_uc2] = flex.double(fractional_change_for_diagonals(diags[i_uc1],diags[i_uc2])).norm()
    return dist
