from __future__ import print_function

from collections import OrderedDict

import numpy

from scitbx.array_family import flex

from bamboo.common import Info
from bamboo.stats.cluster import find_connected_groups, generate_group_idxs

from giant.structure.select import protein, non_water, find_nearest_atoms
from giant.xray.symmetry import find_symmetry_equivalent_groups


class PointCluster(object):
    def __init__(self, id, points, values):
        """Class to hold information about an identified cluster of points, with associated values"""

        # Check that both are lists for the moment
        points = list(points)
        values = list(values)
        assert len(points) == len(values)

        self.parent = None
        self.id = id
        # Store as flex arrays
        self.points = flex.vec3_double(points)
        self.values = flex.double(values)

        stats = self.values.min_max_mean()
        self.min = stats.min
        self.max = stats.max
        self.mean = stats.mean
        self.size = len(self.values)

        self.peak = self.points[values.index(max(values))]
        self.centroid = self.points.mean()

    def summary(self):
        out = []
        out.append('Point Cluster {!s} Summary'.format(self.id))
        out.append('{!s} Points'.format(self.size))
        out.append('Min:'.format(self.min))
        out.append('Max:'.format(self.max))
        out.append('Mean:'.format(self.mean))
        return '\n'.join(out)


class Event(object):
    _attributes = ['estimated_pseudo_occupancy', 'estimated_bdc', 'global_correlation', 'local_correlation']

    def __init__(self, id, cluster, info=None):
        """Class to hold information about an event in a dataset"""
        # Process cluster object (points and values for Event)
        assert isinstance(cluster, PointCluster), 'cluster must be of type PointCluster'
        self.cluster = cluster
        self.cluster.parent = self
        # Give it a name
        if id:
            self.id = id
        else:
            self.id = cluster.id
        # Allow a parent
        self.parent = None
        # Add Meta to the object
        if info:
            assert isinstance(info, Info)
            for a in self._attributes: assert hasattr(info, a)
            self.info = info
        else:
            self.info = Info(self._attributes)

    def summary(self):
        out = []
        out.append('Event {!s} Summary'.format(self.id))
        out.append(self.cluster.summary())
        out.append(self.info.summary())
        return '\n'.join(out)


class PanDDADefaultEventFinder:

    def __init__(self, min_blob_z_peak=0, grid_minimum_volume=0, grid_spacing=1, grid=None):

        self.min_blob_z_peak = min_blob_z_peak
        self.grid_minimum_volume = int(float(grid_minimum_volume) / float(grid_spacing) ** 3)
        self.grid_spacing = grid_spacing

        self.grid = None

    def __call__(self, dataset, clusters, grid):

        num_clusters = clusters[0]
        z_clusters = clusters[1]

        # TODO: Not very functional...
        self.grid = grid

        print(self.grid_minimum_volume)

        log_strs = []
        print(num_clusters)
        # ============================================================================>
        # Too many points to cluster -- probably a bad dataset
        # ============================================================================>
        if num_clusters == -1:
            # TODO: changed this hard, make sure upstream changes still consistent
            return None

        # ============================================================================>
        #####
        # FILTER/SELECT CLUSTERS OF Z-SCORES
        #####
        # ============================================================================>
        # Filter the clusters by size and peak height
        # ============================================================================>
        if num_clusters > 0:
            num_clusters, z_clusters = self.filter_z_clusters_1(z_clusters=z_clusters)
            self.validate_clusters(z_clusters)
            if num_clusters == 0: log_strs.append('===> Minimum cluster peak/size not reached.')
            print(num_clusters)
        # ============================================================================>
        # Filter the clusters by distance from protein
        # ============================================================================>
        if num_clusters > 0:
            num_clusters, z_clusters = self.filter_z_clusters_2(z_clusters=z_clusters, dataset=dataset)
            self.validate_clusters(z_clusters)
            if num_clusters == 0: log_strs.append('===> Clusters too far from protein.')
            print(num_clusters)
        # ============================================================================>
        # Group Nearby Clusters Together
        # ============================================================================>
        if num_clusters > 0:
            num_clusters, z_clusters = self.group_clusters(z_clusters=z_clusters)
            self.validate_clusters(z_clusters)
            print(num_clusters)
        # ============================================================================>
        # Filter the clusters by symmetry equivalence
        # ============================================================================>
        if num_clusters > 0:
            num_clusters, z_clusters = self.filter_z_clusters_3(z_clusters=z_clusters, dataset=dataset)
            self.validate_clusters(z_clusters)
            print(num_clusters)

        events = self.make_clusters(dataset, z_clusters)
        print(num_clusters)
        return num_clusters, z_clusters, events

    def validate_clusters(self, z_clusters):
        for i, (gps, vals) in enumerate(z_clusters):
            assert len(gps) == len(vals)

    def filter_z_clusters_1(self, z_clusters):
        """Filter the z-clusters on a variety of criteria (size, peak value)"""
        filt_z_clusters = z_clusters

        # Filter out small clusters - get numbers of clusters satisfying the minimum cluster size
        large_clusters = (
                flex.int([x[1].size() for x in filt_z_clusters]) >= int(self.grid_minimum_volume)).iselection()
        if large_clusters.size() == 0:  return 0, []
        filt_z_clusters = [filt_z_clusters[i] for i in large_clusters]
        # Filter out weak clusters - get numbers of clusters satisfying the minimum z_peak value
        strong_clusters = (flex.double(
            [x[1].min_max_mean().max for x in filt_z_clusters]) >= self.min_blob_z_peak).iselection()
        if strong_clusters.size() == 0: return 0, []
        filt_z_clusters = [filt_z_clusters[i] for i in strong_clusters]

        return len(filt_z_clusters), filt_z_clusters

    def filter_z_clusters_2(self, z_clusters, dataset, min_contact_dist=6):
        """Find and remove clusters more than a minimum distance from the protein"""

        # min_contact_dist - blobs are rejected if they are more than this distance from the protein

        # Extract the protein sites in the reference frame
        ref_sites_cart = dataset.model.alignment.nat2ref(protein(dataset.model.hierarchy).atoms().extract_xyz())
        # Save time - calculate the square of the contact distance
        min_contact_dist_sq = min_contact_dist ** 2

        # Remove any clusters that are more than min_contact_dist from the protein
        filtered_c_idxs = []
        for c_idx, (c_gps, c_val) in enumerate(z_clusters):
            # Extract points in cluster
            cluster_points_cart = self.grid.grid2cart(c_gps)
            # Calculate minimum distance to protein
            for r_site_cart in ref_sites_cart:
                diff_vecs_cart = cluster_points_cart - r_site_cart
                # Keep cluster if minimum distance is less than min_contact_dist
                if min(diff_vecs_cart.dot()) < min_contact_dist_sq:
                    filtered_c_idxs.append(c_idx)
                    break

        filt_z_clusters = [z_clusters[i] for i in filtered_c_idxs]

        return len(filt_z_clusters), filt_z_clusters

    def filter_z_clusters_3(self, z_clusters, dataset, max_contact_dist=8):
        """Find and remove symmetry equivalent clusters"""

        if len(z_clusters) == 1:
            return 1, z_clusters

        # Extract the protein sites in the reference frame
        d_sites_cart = protein(dataset.model.hierarchy).atoms().extract_xyz()
        d_unit_cell = dataset.model.unit_cell
        d_sym_ops = dataset.model.crystal_contact_operators()

        # Cartesianise and fractionalise the points in each of the clusters (in the crystallographic frame)
        points_cart = [None] * len(z_clusters)
        points_frac = [None] * len(z_clusters)
        for c_idx, (c_gps, c_val) in enumerate(z_clusters):
            # Extract points in cluster
            points_cart[c_idx] = dataset.model.alignment.ref2nat(self.grid.grid2cart(c_gps))
            # Fractionalise them to the unit cell of the dataset
            points_frac[c_idx] = d_unit_cell.fractionalize(points_cart[c_idx])
        # Find the sets of clusters that are symmetry related
        sym_equiv_groups = find_symmetry_equivalent_groups(points_frac=points_frac,
                                                           sym_ops=d_sym_ops,
                                                           unit_cell=d_unit_cell,
                                                           cutoff_cart=1.05 * 1.7321 * self.grid_spacing)
        # max_contact_dist - a point contacts an atom if the atoms is within this distance of it
        # Save time - calculate the square of the contact distance
        max_contact_dist_sq = max_contact_dist ** 2
        # Iterate through and chose one from each group to keep
        filt_z_clusters = []
        for g_id, g_idxs in generate_group_idxs(sym_equiv_groups):
            # Count the number of contact for each cluster in the group
            c_contacts = []
            # Iterate through cluster in the group
            for c_idx in g_idxs:
                # Initialise contact counter
                contacts = 0
                # Get the cartesian points for the cluster
                c_points_cart = points_cart[c_idx]
                # Again, use the brute force all-v-all method
                for rp in d_sites_cart:
                    diffs_cart = c_points_cart - rp
                    # Check to see if site closer to cluster than minimum
                    if min(diffs_cart.dot()) < max_contact_dist_sq:
                        contacts += 1
                # Record the number of contacts (over size of cluster)
                c_contacts.append(1.0 * contacts / len(c_points_cart))
            #                if self.log.verbose:
            #                    print('CLUSTER:', c_idx, ', CONTACTS PER POINT:', round(c_contacts[-1],3))

            # Find the cluster with the most contacts
            max_contacts = max(c_contacts)
            if max_contacts == 0:
                raise Exception('MAX CONTACTS IS 0!')
            else:
                cluster_to_keep = g_idxs[c_contacts.index(max_contacts)]
                filt_z_clusters.append(z_clusters[cluster_to_keep])
        #                if self.log.verbose:
        #                    print('KEEPING CLUSTER', cluster_to_keep)
        assert len(filt_z_clusters) == max(
            sym_equiv_groups), 'NUMBER OF UNIQUE GROUPS AND GROUPS TO BE RETURNED NOT THE SAME'

        return len(filt_z_clusters), filt_z_clusters

    def group_clusters(self, z_clusters, separation_cutoff=5):
        """Join clusters that are separated by less than max_separation"""

        if len(z_clusters) == 1:
            return 1, z_clusters

        # Minimum distance between grid points to be joined (squared)
        grid_cutoff_sq = (separation_cutoff / self.grid_spacing) ** 2

        # Record which clusters are to be joined
        connect_array = numpy.zeros((len(z_clusters), len(z_clusters)), dtype=int)
        for i_clust_1, (c_gps_1, c_val_1) in enumerate(z_clusters):
            for i_clust_2, (c_gps_2, c_val_2) in enumerate(z_clusters):
                # Skip if this is the same blob
                if i_clust_1 == i_clust_2:
                    connect_array[(i_clust_1, i_clust_2)] = 1
                    continue
                # Extract the minimum separation of the grid points
                min_dist_sq = min([min((c_gps_2 - gp).dot()) for gp in c_gps_1])
                # Check to see if they should be joined
                if min_dist_sq < grid_cutoff_sq:
                    connect_array[(i_clust_1, i_clust_2)] = 1
        # Cluster the connection array
        cluster_groupings = find_connected_groups(connection_matrix=connect_array)
        # Concatenate smaller clusters into larger clusters
        grouped_clusters = []
        for g_id, g_idxs in generate_group_idxs(cluster_groupings):
            g_gps = []
            [g_gps.extend(z_clusters[i][0]) for i in g_idxs]
            g_gps = flex.vec3_double(g_gps)
            g_val = []
            [g_val.extend(z_clusters[i][1]) for i in g_idxs]
            g_val = flex.double(g_val)
            grouped_clusters.append((g_gps, g_val))

        assert len(grouped_clusters) == max(cluster_groupings)

        return len(grouped_clusters), grouped_clusters

    def make_clusters(self, dataset, z_clusters):

        events = []

        # ============================================================================>
        # Process the identified features
        # ============================================================================>
        for event_idx, (event_points, event_values) in enumerate(z_clusters):
            # Number events from 1
            event_num = event_idx + 1
            # Create a unique identifier for this event
            event_key = (dataset.tag, event_num)
            # ============================================================================>
            # Create a point cluster object
            # ============================================================================>
            point_cluster = PointCluster(id=event_key, points=event_points, values=event_values)
            # ============================================================================>
            # Create an event object
            # ============================================================================>
            event_obj = Event(id=point_cluster.id, cluster=point_cluster)
            event_obj.info.estimated_pseudo_occupancy = 0
            event_obj.info.estimated_bdc = 1.0 - 0
            event_obj.info.global_correlation = 0
            event_obj.info.local_correlation = 0
            # ============================================================================>
            # Append to dataset handler
            # ============================================================================>
            dataset.events.append(event_obj)
            events.append(event_obj)

        return events

    def repr(self):
        repr = OrderedDict()
        return repr
