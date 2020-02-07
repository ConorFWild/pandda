import numpy
from scitbx.array_family import flex

def pairwise_dists(xyz1, xyz2):
    """Calculate pairwise distances between xyz1 and xyz2"""
    return numpy.array([flex.sqrt((xyz2 - c1).dot()) for c1 in xyz1])

def is_within(dist, coords_1, coords_2):
    """Checks if any of coords_1 is within dist of coords_2"""

    dist_sq = dist**2
    if len(coords_2) < len(coords_1):
        i1 = coords_2; i2 = coords_1
    else:
        i1 = coords_1; i2 = coords_2

    for c1 in i1:
        diffs = i2 - c1
        min_d_sq = flex.min(diffs.dot())
        if min_d_sq < dist_sq:
            return True
    return False
