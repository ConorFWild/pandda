import itertools
import scipy.spatial
import scipy.cluster
import numpy

def generate_group_idxs(group_vals):
    num_vals = len(group_vals)
    # Sort the list idxs by group val
    sorted_idx = sorted(range(num_vals), key=lambda i: group_vals[i])
    # Reorder the group vals
    sorted_val = [group_vals[i] for i in sorted_idx]
    # For each group val, return idxs for this group
    for val, sel_idxs in itertools.groupby(range(num_vals), key=lambda i: sorted_val[i]):
        yield val, [sorted_idx[idx] for idx in sel_idxs]

def find_connected_groups(connection_matrix):
    """Take a matrix of connected elements and output groups"""

    assert not set(connection_matrix.flatten()).difference(set([0,1])), 'CONNECTION MATRIX MUST CONSIST OF 1s AND 0s ONLY'
    assert connection_matrix.ndim == 2, 'CONNECTION MATRIX MUST BE OF DIMENSION 2'
    assert connection_matrix.shape[0] > 1, 'MATRIX MUST BE LARGER THAN 1x1'
    assert connection_matrix.shape[0] == connection_matrix.shape[1], 'CONNECTION MATRIX MUST BE SQUARE'

    # Make it symmetrical - if 1 is connected to 2, 2 is connected to 1
    sym_connection_matrix = ((connection_matrix + connection_matrix.T) != 0).astype(int)
    # Convert from a connection array to a distance array (0 for connected, 1 for unconnected)
    dist_mx = 1 - sym_connection_matrix
    # Convert to condensed distance matrix
    cd_dist_mx = scipy.spatial.distance.squareform(dist_mx)
    # Calculate the linkage matrix
    l = scipy.cluster.hierarchy.linkage(cd_dist_mx, method='single', metric='euclidean')
    # Cluster with very low cutoff
    connected_groups = scipy.cluster.hierarchy.fcluster(Z=l, t=0.01, criterion='distance').tolist()

    return connected_groups

