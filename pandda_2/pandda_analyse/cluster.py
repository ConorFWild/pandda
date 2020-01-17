import time

import numpy
import scipy.cluster
from hdbscan import HDBSCAN

from scitbx.array_family import flex

from bamboo.stats.cluster import find_connected_groups, generate_group_idxs

from giant.grid.utils import idx_to_grid


def cluster_z_scores(z_map_data, point_mask_idx, grid_clustering_cutoff, log=None, negative_values=False,
                     cluster_method="hierarchical", contour_level=0):
    """

    Parameters
    ----------
    z_map_data : flex (p,q,r)
    point_mask_idx : flex (p*q*r)
    grid_clustering_cutoff
    log
    negative_values
    cluster_method
    contour_level

    Returns
    -------

    """
    """Finds all the points in the z-map above `z_cutoff`, points will then be clustered into groups of cutoff `clustering_cutoff` angstroms"""

    # Truncate map
    above_val, above_idx, above_gps = truncate_map(z_map_data,
                                                   point_mask_idx,
                                                   negative_values,
                                                   contour_level)

    # Check valid
    # No Cluster points found
    if len(above_val) == 0:
        return 0, []
    # One Cluster point found
    elif len(above_val) == 1:
        return 1, [(above_gps, above_val)]

    # Cluster
    cluter_time_start = time.time()
    if cluster_method == "hdbscan":
        cluster_ids = cluster_hdbscan(above_gps)
    elif cluster_method == "hierarchical":
        cluster_ids = cluster_hierarcical(above_gps, grid_clustering_cutoff)
    else:
        raise Exception("No valid clustering method specified!")
    cluter_time_finish = time.time()

    # Get number of clusters
    num_clusters = get_num_clusters(cluster_ids)

    # Group clusters
    z_clusters = group_clusters(cluster_ids,
                                num_clusters,
                                above_gps,
                                above_val)

    if log is not None:
        log('> Clustered {!s} Points.'.format(len(above_val)))
        log('> Clustering > Time Taken: {!s} seconds'.format(int(cluter_time_finish - cluter_time_start)))

    return num_clusters, z_clusters


def mask_map(z_map_data, point_mask_idx):
    # Select these values from the map
    point_mask_idx = flex.size_t(point_mask_idx)
    point_mask_val = z_map_data.select(point_mask_idx)
    return point_mask_val


def truncate_map(z_map_data, point_mask_idx, negative_values=False, contour_level=0):

    # Mask map
    z_map_masked = mask_map(z_map_data, point_mask_idx)

    # truncate map
    z_truncated_map_truncated = get_truncated_map_mask(z_map_masked,
                                                       negative_values,
                                                       contour_level)

    # Extract these points
    above_val, above_idx, above_gps = extract_points(z_map_masked,
                                                     point_mask_idx,
                                                     z_truncated_map_truncated,
                                                     z_map_data)

    return above_val, above_idx, above_gps


def get_truncated_map_mask(point_mask_val, negative_values=False, contour_level=0):
    # Find values above cutoff
    if negative_values:
        above_idx = (point_mask_val >= contour_level).iselection()
        below_idx = (point_mask_val <= -1.0 * contour_level).iselection()
        sel_idx = above_idx.concatenate(below_idx)
    else:
        sel_idx = (point_mask_val >= contour_level).iselection()

    return sel_idx


def extract_points(point_mask_val, point_mask_idx, sel_idx, z_map_data):
    # Extract values and grid points for these sites
    above_val = point_mask_val.select(sel_idx)
    above_idx = point_mask_idx.select(sel_idx)
    above_gps = flex.vec3_double([idx_to_grid(i, grid_size=z_map_data.all()) for i in above_idx])
    # above_len = len(above_val)

    return above_val, above_idx, above_gps


def cluster_hdbscan(above_gps):
    sample_by_feature = above_gps.to_numpy()
    print(sample_by_feature.shape)

    clusterer = HDBSCAN()
    clusterer.fit(sample_by_feature)
    cluster_ids = list(clusterer.labels_)

    return cluster_ids


def cluster_hierarcical(above_gps, grid_clustering_cutoff):
    cluster_ids = scipy.cluster.hierarchy.fclusterdata(X=above_gps,
                                                       t=grid_clustering_cutoff,
                                                       criterion='distance',
                                                       metric='euclidean',
                                                       method='single')
    cluster_ids = list(cluster_ids)
    return cluster_ids


def get_num_clusters(cluster_ids):
    # Get the number of clusters
    cluster_id_array = numpy.array(cluster_ids)
    unique_cluster_ids = numpy.unique(cluster_id_array)
    num_clusters = len(unique_cluster_ids)

    return num_clusters


def group_clusters(cluster_ids, num_clusters, above_gps, above_val):
    # Group the values by cluster id
    z_clusters = []
    for c_id, c_idxs in generate_group_idxs(cluster_ids):
        c_idxs = flex.size_t(c_idxs)
        c_gps = above_gps.select(c_idxs)
        c_val = above_val.select(c_idxs)
        z_clusters.append((c_gps, c_val))
    assert num_clusters == len(z_clusters)

    return z_clusters





