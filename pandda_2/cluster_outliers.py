from __future__ import print_function

from collections import OrderedDict

import time

import numpy
import scipy

from scitbx.array_family import flex

from bamboo.stats.cluster import find_connected_groups, generate_group_idxs

from giant.grid.utils import idx_to_grid
from giant.structure.select import protein, non_water, find_nearest_atoms
from giant.grid.masks import GridMask


class PanDDADefaultCluster:

    def __init__(self, grid_clustering_cutoff, negative_values, cluster_method, contour_level, outer_mask,
                 inner_mask_symmetry, grid=None, dataset_mask=None):
        self.grid_clustering_cutoff = grid_clustering_cutoff
        self.negative_values = negative_values
        self.cluster_method = cluster_method
        self.contour_level = contour_level

        self.outer_mask = outer_mask
        self.inner_mask_symmetry = inner_mask_symmetry

    def __call__(self, dataset, z_map, grid):

        z_map_data = z_map.get_map_data(sparse=False)

        # ============================================================================>
        # Extract the global mask object from the grid
        # ============================================================================>
        dset_total_temp = grid.global_mask().total_mask_binary().copy()

        # ============================================================================>
        # Generate symmetry masks for this dataset
        # ============================================================================>
        # Generate symmetry contacts for this dataset and align to reference frame
        dataset_sym_copies = dataset.model.crystal_contacts(distance_cutoff=self.outer_mask + 5,
                                                            combine_copies=True)
        dataset_sym_copies.atoms().set_xyz(
            dataset.model.alignment.nat2ref(dataset_sym_copies.atoms().extract_xyz()))

        # Extract protein atoms from the symmetry copies
        dataset_sym_sites_cart = non_water(dataset_sym_copies).atoms().extract_xyz()
        # Generate symmetry contacts grid mask
        dataset_mask = GridMask(parent=grid,
                                sites_cart=dataset_sym_sites_cart,
                                max_dist=self.outer_mask,
                                min_dist=self.inner_mask_symmetry)

        # Combine with the total mask to generate custom mask for this dataset
        dset_total_temp.put(dataset_mask.inner_mask_indices(), 0)
        point_mask_idx = numpy.where(dset_total_temp)[0]

        # Truncate map
        above_val, above_idx, above_gps = self.truncate_map(z_map_data,
                                                            point_mask_idx,
                                                            self.negative_values,
                                                            self.contour_level)

        # Check valid
        # No Cluster points found
        if len(above_val) == 0:
            return 0, []
        # One Cluster point found
        elif len(above_val) == 1:
            return 1, [(above_gps, above_val)]

        # Cluster
        cluter_time_start = time.time()
        cluster_ids = self.cluster_hierarcical(above_gps, self.grid_clustering_cutoff)
        cluter_time_finish = time.time()

        # Get number of clusters
        num_clusters = self.get_num_clusters(cluster_ids)

        # Group clusters
        z_clusters = self.group_clusters(cluster_ids,
                                         num_clusters,
                                         above_gps,
                                         above_val)

        return num_clusters, z_clusters

    def mask_map(self, z_map_data, point_mask_idx):
        # Select these values from the map
        point_mask_idx = flex.size_t(point_mask_idx)
        point_mask_val = z_map_data.select(point_mask_idx)
        return point_mask_val

    def truncate_map(self, z_map_data, point_mask_idx, negative_values=False, contour_level=0):

        # Mask map
        z_map_masked = self.mask_map(z_map_data, point_mask_idx)

        # truncate map
        z_truncated_map_truncated = self.get_truncated_map_mask(z_map_masked,
                                                                negative_values,
                                                                contour_level)

        # Extract these points
        above_val, above_idx, above_gps = self.extract_points(z_map_masked,
                                                              point_mask_idx,
                                                              z_truncated_map_truncated,
                                                              z_map_data)

        return above_val, above_idx, above_gps

    def get_truncated_map_mask(self, point_mask_val, negative_values=False, contour_level=0):
        # Find values above cutoff
        if negative_values:
            above_idx = (point_mask_val >= contour_level).iselection()
            below_idx = (point_mask_val <= -1.0 * contour_level).iselection()
            sel_idx = above_idx.concatenate(below_idx)
        else:
            sel_idx = (point_mask_val >= contour_level).iselection()

        return sel_idx

    def extract_points(self, point_mask_val, point_mask_idx, sel_idx, z_map_data):
        point_mask_idx = flex.size_t(point_mask_idx)
        # Extract values and grid points for these sites
        above_val = point_mask_val.select(sel_idx)
        above_idx = point_mask_idx.select(sel_idx)
        above_gps = flex.vec3_double([idx_to_grid(i, grid_size=z_map_data.all()) for i in above_idx])
        # above_len = len(above_val)

        return above_val, above_idx, above_gps

    # def cluster_hdbscan(self, above_gps):
    #     sample_by_feature = above_gps.to_numpy()
    #     print(sample_by_feature.shape)
    #
    #     clusterer = HDBSCAN()
    #     clusterer.fit(sample_by_feature)
    #     cluster_ids = list(clusterer.labels_)
    #
    #     return cluster_ids

    def cluster_hierarcical(self, above_gps, grid_clustering_cutoff):
        cluster_ids = scipy.cluster.hierarchy.fclusterdata(X=above_gps,
                                                           t=grid_clustering_cutoff,
                                                           criterion='distance',
                                                           metric='euclidean',
                                                           method='single')
        cluster_ids = list(cluster_ids)
        return cluster_ids

    def get_num_clusters(self, cluster_ids):
        # Get the number of clusters
        cluster_id_array = numpy.array(cluster_ids)
        unique_cluster_ids = numpy.unique(cluster_id_array)
        num_clusters = len(unique_cluster_ids)

        return num_clusters

    def group_clusters(self, cluster_ids, num_clusters, above_gps, above_val):
        # Group the values by cluster id
        z_clusters = []
        for c_id, c_idxs in generate_group_idxs(cluster_ids):
            c_idxs = flex.size_t(c_idxs)
            c_gps = above_gps.select(c_idxs)
            c_val = above_val.select(c_idxs)
            z_clusters.append((c_gps, c_val))
        assert num_clusters == len(z_clusters)

        return z_clusters

    def repr(self):
        repr = OrderedDict()
        repr["grid_clustering_cutoff"] = self.grid_clustering_cutoff
        repr["negative_values"] = self.negative_values
        repr["cluster_method"] = self.cluster_method
        repr["contour_level"] = self.contour_level
        repr["outer_mask"] = self.outer_mask
        repr["inner_mask_symmetry"] = self.inner_mask_symmetry

        return repr
