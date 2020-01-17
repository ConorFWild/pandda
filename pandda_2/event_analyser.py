from __future__ import print_function

from collections import OrderedDict

import numpy

from scitbx.array_family import flex

from giant.structure.select import protein, non_water, find_nearest_atoms
from giant.grid.masks import GridMask
from giant.xray.maps.bdc import calculate_varying_bdc_correlations, calculate_maximum_series_discrepancy, \
    calculate_bdc_subtracted_map


class EventAnalyser:
    def __init__(self,
                 max_bdc=1.0,
                 min_bdc=0.0,
                 increment=0.01,
                 output_multiplier=1.0, ):
        self.max_bdc = max_bdc
        self.min_bdc = min_bdc
        self.increment = increment
        self.output_multiplier = output_multiplier

    def __call__(self,
                 dataset,
                 dataset_map,
                 ref_map,
                 events,
                 grid,

                 ):
        # ============================================================================>
        # Extract the map data in non-sparse format
        # ============================================================================>
        dset_map_data = dataset_map.get_map_data(sparse=False)
        ref_map_data = ref_map.get_map_data(sparse=False)
        # ============================================================================>
        # Unpack cluster
        # ============================================================================>
        event_stats = OrderedDict()
        for event in events[2]:
            # ============================================================================>
            # Estimate the background correction of the detected feature
            # ============================================================================>
            # Extract sites for this cluster and estimate the background correction for the event

            # Generate custom grid mask for this dataset
            event_mask = GridMask(parent=grid,
                                  sites_cart=grid.grid2cart(event.cluster.points,
                                                            origin_shift=True),
                                  max_dist=2.0,
                                  min_dist=0.0,
                                  )

            # Select masks to define regions for bdc calculation
            exp_event_idxs = flex.size_t(event_mask.outer_mask_indices())
            reference_idxs = flex.size_t(grid.global_mask().inner_mask_indices())
            # ============================================================================>
            # Generate BDC-estimation curve and estimate BDC
            # ============================================================================>
            event_remains, event_corrs, global_corrs = calculate_varying_bdc_correlations(
                ref_map_data=ref_map_data,
                query_map_data=dset_map_data,
                feature_idxs=exp_event_idxs,
                reference_idxs=reference_idxs,
                min_remain=1.0 - self.max_bdc,
                max_remain=1.0 - self.min_bdc,
                bdc_increment=self.increment,
                verbose=True)
            event_remain_est = calculate_maximum_series_discrepancy(
                labels=event_remains,
                series_1=global_corrs,
                series_2=event_corrs)

            event_remain_est = min(event_remain_est * self.output_multiplier,
                                   1.0 - self.min_bdc)
            # ============================================================================>
            # Calculate the map correlations at the selected BDC
            # ============================================================================>
            event_map_data = calculate_bdc_subtracted_map(
                ref_map_data=ref_map_data,
                query_map_data=dset_map_data,
                bdc=1.0 - event_remain_est)
            global_corr = \
                numpy.corrcoef(event_map_data.select(reference_idxs), ref_map_data.select(reference_idxs))[
                    0, 1]
            local_corr = numpy.corrcoef(event_map_data.select(exp_event_idxs), ref_map_data.select(exp_event_idxs))[
                0, 1]
            # ============================================================================>
            # Update event parameters
            # ============================================================================>
            event.info.estimated_pseudo_occupancy = event_remain_est
            event.info.estimated_bdc = 1.0 - event_remain_est
            event.info.global_correlation = global_corr
            event.info.local_correlation = local_corr

            # ============================================================================>
            # Find the nearest atom to the event
            # ============================================================================>
            # TODO: restore this?
            atm = find_nearest_atoms(atoms=list(protein(dataset.model.hierarchy).atoms_with_labels()),
                                     query=dataset.model.alignment.ref2nat(
                                         grid.grid2cart(sites_grid=[map(int, event.cluster.centroid)],
                                                        origin_shift=True)))[0]

            event_stats[event.id] = OrderedDict()
            event_stats[event.id]["estimated_pseudo_occupancy"] = event_remain_est
            event_stats[event.id]["estimated_bdc"] = 1.0 - event_remain_est
            event_stats[event.id]["global_corr"] = global_corr
            event_stats[event.id]["local_corr"] = global_corr

        return event_stats

    def repr(self):
        repr = OrderedDict()
        repr["max_bdc"] = self.max_bdc
        repr["min_bdc"] = self.min_bdc
        repr["increment"] = self.increment
        repr["output_multiplier"] = self.output_multiplier

        return repr
