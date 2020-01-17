from __future__ import print_function

import gc

from collections import OrderedDict
import pathlib as p

import numpy
import pandas as pd

from scitbx.array_family import flex

from bamboo.common.path import rel_symlink, splice_ext

from giant.grid.masks import GridMask
from giant.structure.select import protein, non_water, find_nearest_atoms
from giant.xray.maps.bdc import calculate_varying_bdc_correlations, calculate_maximum_series_discrepancy, calculate_bdc_subtracted_map
from giant.xray.maps.local_align import create_native_map

from pandda.analyse.events import PointCluster, Event
from pandda.analyse import graphs as analyse_graphs


#
# class PanDDADefaultCriticism:
#
#     def __init__(self):
#
#         self.criticisms = OrderedDict()
#
#
# class OutputMaps:
#
#     def __init__(self):
#         self.map_maker = 1
#
#     def __call__(self, model, tree):
#
#         dataset_paths = {path.name: path
#                          for path
#                          in tree(("processed_datasets", ".*"))[0]}
#         print(dataset_paths)
#
#         for dtag, m in model.statistical_maps:
#             self.map_maker(m, dataset_paths[dtag] / "z_map.ccp4")
#
#
#     def __call__(self, dataset, dataset_map, ref_map_data, z_clusters, grid, args, verbose, log_strs, blob_finder, ref_ext):
#
#         # ============================================================================>
#         # Extract the map data in non-sparse format
#         # ============================================================================>
#         dset_map_data = dataset_map.get_map_data(sparse=False)
#         # ============================================================================>
#         # Process the identified features
#         # ============================================================================>
#         for event_idx, (event_points, event_values) in enumerate(z_clusters):
#             # Number events from 1
#             event_num = event_idx + 1
#             # Create a unique identifier for this event
#             event_key = (dataset.tag, event_num)
#             # ============================================================================>
#             # Create a point cluster object
#             # ============================================================================>
#             point_cluster = PointCluster(id=event_key, points=event_points, values=event_values)
#             # ============================================================================>
#             # Estimate the background correction of the detected feature
#             # ============================================================================>
#             # Extract sites for this cluster and estimate the background correction for the event
#
#             # Generate custom grid mask for this dataset
#             event_mask = GridMask(parent=grid, sites_cart=grid.grid2cart(point_cluster.points, origin_shift=True),
#                                   max_dist=2.0, min_dist=0.0)
#
#             # Select masks to define regions for bdc calculation
#             exp_event_idxs = flex.size_t(event_mask.outer_mask_indices())
#             reference_idxs = flex.size_t(grid.global_mask().inner_mask_indices())
#             # ============================================================================>
#             # Generate BDC-estimation curve and estimate BDC
#             # ============================================================================>
#             event_remains, event_corrs, global_corrs = calculate_varying_bdc_correlations(
#                 ref_map_data=ref_map_data,
#                 query_map_data=dset_map_data,
#                 feature_idxs=exp_event_idxs,
#                 reference_idxs=reference_idxs,
#                 min_remain=1.0 - args.params.background_correction.max_bdc,
#                 max_remain=1.0 - args.params.background_correction.min_bdc,
#                 bdc_increment=args.params.background_correction.increment,
#                 verbose=verbose)
#             event_remain_est = calculate_maximum_series_discrepancy(
#                 labels=event_remains,
#                 series_1=global_corrs,
#                 series_2=event_corrs)
#             analyse_graphs.write_occupancy_graph(
#                 f_name=dataset.file_manager.get_file('bdc_est_png').format(event_num),
#                 x_values=event_remains,
#                 global_values=global_corrs,
#                 local_values=event_corrs)
#             event_remain_est = min(event_remain_est * args.params.background_correction.output_multiplier,
#                                    1.0 - args.params.background_correction.min_bdc)
#             # ============================================================================>
#             # Calculate the map correlations at the selected BDC
#             # ============================================================================>
#             event_map_data = calculate_bdc_subtracted_map(
#                 ref_map_data=ref_map_data,
#                 query_map_data=dset_map_data,
#                 bdc=1.0 - event_remain_est)
#             global_corr = numpy.corrcoef(event_map_data.select(reference_idxs), ref_map_data.select(reference_idxs))[
#                 0, 1]
#             local_corr = numpy.corrcoef(event_map_data.select(exp_event_idxs), ref_map_data.select(exp_event_idxs))[
#                 0, 1]
#             # ============================================================================>
#             # Write out EVENT map (in the reference frame) and grid masks
#             # ============================================================================>
#             if args.output.developer.write_reference_frame_dataset_masks_and_maps:
#                 # Event map
#                 f_name = splice_ext(dataset.file_manager.get_file('event_map').format(event_num, event_remain_est),
#                                     ref_ext, position=-1)
#                 event_map = dataset_map.new_from_template(event_map_data, sparse=False)
#                 event_map.to_file(filename=f_name, space_group=grid.space_group())
#                 # Event map mask
#                 f_name = splice_ext(dataset.file_manager.get_file('grid_mask'),
#                                     '-event-mask-{}.{}'.format(event_num, ref_ext), position=-1)
#                 grid.write_indices_as_map(indices=event_mask.outer_mask_indices(), f_name=f_name)
#
#             # ============================================================================>
#             # Find the nearest atom to the event
#             # ============================================================================>
#             atm = find_nearest_atoms(atoms=list(protein(dataset.model.hierarchy).atoms_with_labels()),
#                                      query=dataset.model.alignment.ref2nat(
#                                          grid.grid2cart(sites_grid=[map(int, point_cluster.centroid)],
#                                                         origin_shift=True)))[0]
#
#
# class EventMapMaker:
#
#     def __init__(self, grid):
#         self.grid = grid
#
#     def __call__(self, dataset_map, ref_map, event, bdc, file_path):
#
#         # ============================================================================>
#         # Extract the map data in non-sparse format
#         # ============================================================================>
#         dset_map_data = dataset_map.get_map_data(sparse=False)
#         ref_map_data = ref_map.get_map_data(sparse=False)
#         # ============================================================================>
#         # Calculate the map correlations at the selected BDC
#         # ============================================================================>
#         event_map_data = calculate_bdc_subtracted_map(
#             ref_map_data=ref_map_data,
#             query_map_data=dset_map_data,
#             bdc=bdc)
#         # ============================================================================>
#         # Write out EVENT map (in the reference frame) and grid masks
#         # ============================================================================>
#         # Event map
#         event_map = dataset_map.new_from_template(event_map_data, sparse=False)
#         event_map.to_file(filename=file_path, space_group=self.grid.space_group())
#         # Event map mask
#         # TODO: restore this?
#         # f_name = splice_ext(dataset.file_manager.get_file('grid_mask'),
#         #                     '-event-mask-{}.{}'.format(event_num, ref_ext), position=-1)
#         # grid.write_indices_as_map(indices=event_mask.outer_mask_indices(), f_name=f_name)

def make_z_map(xmap, truncated_dataset, z_map_path, statistical_model, grid):
    statistical_map = statistical_model.evaluate(xmap)  # sample)
    statistical_map_maker = PanddaNativeStatisticalMapMaker(grid,
                                                            grid.get_mask("protein")._max_dist,
                                                            grid.grid_spacing(),
                                                            )
    file_path = p.Path(z_map_path)
    if not file_path.exists():
        file_path_string = str(file_path)
        print("Outputing statistical map")
        statistical_map_maker(truncated_dataset,
                              statistical_map,
                              file_path_string,
                              )


def make_event_map(xmap, truncated_dataset, ref_map, event, bdc, event_map_path, grid):

        # Calculate event maps and output
        event_map_maker = PanddaNativeEventMapMaker(grid,
                                                    grid.get_mask("protein")._max_dist,
                                                    grid.grid_spacing(),
                                                    )

        gc.collect()
        print("Outputing event map")
        file_path = p.Path(event_map_path)

        if not file_path.exists():
            file_path_string = str(file_path)
            event_map_maker(truncated_dataset,
                            xmap,  # sample
                            ref_map,
                            event,
                            bdc,
                            str(file_path_string),
                            )


def make_mean_map(dataset,
                  ref_map,
                  grid,
                  mean_map_path,
                  ):
    map_maker = NativeMapMaker(dataset=dataset,
                               map_obj=ref_map,
                               sites_mask=grid.global_mask().sites_cart,
                               filename=str(mean_map_path),
                               outer_mask=grid.get_mask("protein")._max_dist,
                               grid_spacing=grid.grid_spacing(),
                               )

    map_maker.run()


class PanDDADefaultMapMaker:

    def __init__(self, statistical_model=None):
        self.grid = None
        self.outer_mask = None
        self.grid_spacing = None
        self.statistical_model = statistical_model

    def attach_grid(self, grid):
        self.grid = grid
        self.outer_mask = self.grid.get_mask("protein")._max_dist
        self.grid_spacing = self.grid.grid_spacing()

    def process_single(self, xmap, truncated_dataset, ref_map, events, bdcs, dataset_path):
        gc.collect()

        if events[2]:
            # sample = sample_loader(truncated_dataset)

            # Calculate statistical map and output
            print("Making statistical map")
            statistical_map = self.statistical_model.evaluate(xmap  # sample
                                                              )
            statistical_map_maker = PanddaNativeStatisticalMapMaker(self.grid,
                                                                    self.outer_mask,
                                                                    self.grid_spacing)
            file_path = p.Path(dataset_path) / "z_map.ccp4"
            if not file_path.exists():
                file_path_string = str(file_path)
                print("Outputing statistical map")
                statistical_map_maker(truncated_dataset, statistical_map, file_path_string)

            # Calculate event maps and output
            event_map_maker = PanddaNativeEventMapMaker(self.grid,
                                                        self.outer_mask,
                                                        self.grid_spacing)
            for event in events[2]:
                gc.collect()
                print("Outputing event map")
                file_path = p.Path(dataset_path) / "event_map_{}_{}.ccp4".format(event.id[1], bdcs[event.id])
                if not file_path.exists():
                    file_path_string = str(file_path)
                    event_map_maker(truncated_dataset,
                                    xmap,  # sample
                                    ref_map, event, bdcs[event.id], file_path_string)

    def process_shell(self, ref_dataset, ref_map, dataset_path):
        # Calculate shell mean map
        file_path = p.Path(dataset_path) / "ground_map.ccp4"
        if not file_path.exists():
            file_path_string = str(file_path)
            mean_ground_maker = PanDDAGroundMapMaker(self.grid)
            mean_ground_maker(ref_dataset, ref_map, file_path_string)


class PanddaNativeEventMapMaker:

    def __init__(self, grid, outer_mask, grid_spacing):
        self.grid = grid
        self.outer_mask = outer_mask
        self.grid_spacing = grid_spacing

    def __call__(self, dataset, sample, ref_map, event, bdc, file_path):

        # ============================================================================================= #
        #####                                                                                       #####
        #               Generate native-aligned maps (in the crystallographic unit cell)                #
        #####                                                                                       #####
        # ============================================================================================= #
        # ============================================================================>
        # Make Event-map
        # ============================================================================>
        ref_event_map = (sample - ref_map * bdc)
        map_maker = NativeMapMaker(dataset=dataset,
                                   map_obj=ref_event_map,
                                   sites_mask=self.grid.global_mask().sites_cart,
                                   filename=file_path,
                                   outer_mask=self.outer_mask,
                                   grid_spacing=self.grid_spacing
                                   )

        map_maker.run()


class PanddaNativeStatisticalMapMaker:

    def __init__(self, grid, outer_mask, grid_spacing):
        self.grid = grid
        self.outer_mask = outer_mask
        self.grid_spacing = grid_spacing

    def __call__(self, dataset, ref_z_map, file_path):

        # ============================================================================================= #
        #####                                                                                       #####
        #               Generate native-aligned maps (in the crystallographic unit cell)                #
        #####                                                                                       #####
        # ============================================================================================= #
        # ============================================================================>
        # Make Event-map
        # ============================================================================>
        ref_z_map = ref_z_map.normalised_copy()
        map_maker = NativeMapMaker(dataset=dataset,
                                   map_obj=ref_z_map,
                                   sites_mask=self.grid.global_mask().sites_cart,
                                   filename=file_path,
                                   outer_mask=self.outer_mask,
                                   grid_spacing=self.grid_spacing,
                                   )

        map_maker.run()


class PanDDAGroundMapMaker:

    def __init__(self, grid):
        self.grid = grid
        self.outer_mask = self.grid.get_mask("protein")._max_dist
        self.grid_spacing = self.grid.grid_spacing()

    def __call__(self, dataset, ref_map, file_path):

        # ============================================================================================= #
        #####                                                                                       #####
        #               Generate native-aligned maps (in the crystallographic unit cell)                #
        #####                                                                                       #####
        # ============================================================================================= #
        # ============================================================================>
        # Make Event-map
        # ============================================================================>
        map_maker = NativeMapMaker(dataset=dataset,
                                   map_obj=ref_map,
                                   sites_mask=self.grid.global_mask().sites_cart,
                                   filename=file_path,
                                   outer_mask=self.outer_mask,
                                   grid_spacing=self.grid_spacing)

        map_maker.run()


class NativeMapMaker(object):

    def __init__(self, dataset, map_obj, sites_mask, filename, outer_mask, grid_spacing):
        """
        The main object for comparing each dataset-map to the ensemble-maps.
        Constructed so that the class can be initialised and then called within a multiprocessing function.

        Usage:
            new_proc = DatasetProcessor(...)
            output   = new_proc.run()
        or:
            output   = DatasetProcessor.process(...)
        """

        self.data = (dataset, map_obj, sites_mask, filename)
        self.outer_mask = outer_mask
        self.grid_spacing = grid_spacing

    @classmethod
    def process(cls, dataset, map_obj, sites_mask, filename):
        """Process the dataset immediately and return output"""
        return cls(dataset=dataset, map_obj=map_obj, sites_mask=sites_mask, filename=filename).run()

    def run(self):
        """Process the dataset"""

        dataset, map_obj, sites_mask, filename = self.data

        native_map_data = create_native_map(
                            native_crystal_symmetry = dataset.model.crystal_symmetry,
                            native_sites            = dataset.model.alignment.ref2nat(sites_mask),
                            alignment               = dataset.model.alignment,
                            reference_map           = map_obj.make_dense(),
                            site_mask_radius        = self.outer_mask,
                            step                    = self.grid_spacing,
                            filename                = filename
                        )

        return native_map_data


class PanDDADefaultEventTableShell:

    def __init__(self, order_by=None):
        self.order_by = None
        self.event_info = None

    def __call__(self,
                 datasets,
                 events,
                 table_path,
                 grid,
                 events_analysed,):
        event_table = self.make_event_table(datasets,
                                            events,
                                            grid,
                                            events_analysed,
                                            )
        self.write_event_table(table_path,
                               event_table,
                               )

        return event_table

    def make_event_table(self,
                         datasets,
                         events,
                         grid,
                         events_analysed,
                         ):
        # for dtag, dataset in datasets.items():
        #     print(dtag)
        #     # ============================================================================>
        #     # Mark as interesting and add events to the event table
        #     # ============================================================================>
        #     # TODO: dataset.all_masks() will need to be changed
        #     records = []
        #     if dtag in events:
        #         print(events[dtag])
        #         for e in events[dtag][2]:
        #             print(e)
        #             event_record = self.get_event_record(dataset=dataset,
        #                                                  event=e)
        #             records.append(event_record)

        records = []
        for dtag, events in events.items():
            for event in events[2]:
                event_record = self.get_event_record(dataset=datasets[dtag],
                                                     event=event,
                                                     grid=grid,
                                                     analysis=events_analysed[dtag][event.id],
                                                     )
                records.append(event_record)

        return pd.DataFrame(records)

    def get_event_record(self,
                         dataset,
                         event,
                         grid,
                         analysis,
                         ):
        """Add event entries to the event table"""

        # Default to site_idx of 0 if no site given
        if event.parent:
            site_idx = event.parent.id
        else:
            site_idx = 0

        # Get native frame loc
        xyz = list(dataset.model.alignment.ref2nat(
            coordinates=grid.grid2cart([event.cluster.peak],
                                       origin_shift=True,
                                       ))[0]
                   )

        # Generate record
        record = {"dtag": dataset.tag,
                  "event_idx": event.id[1],
                  "site_idx": site_idx,
                  "1-BDC": round(analysis["estimated_pseudo_occupancy"], 2),
                  "z_peak": round(event.cluster.max, 2),
                  "z_mean": round(event.cluster.mean, 2),
                  "cluster_size": event.cluster.size,
                  "x": xyz[0],
                  "y": xyz[1],
                  "z": xyz[2],
                  "global_correlation_to_average_map": analysis["global_corr"],
                  "local_correlation_to_average_map": analysis["local_corr"],
                  }

        return record

    def add_event_to_event_table(self, dataset, event, grid):
        """Add event entries to the event table"""

        # Check event has not been added previously
        self.event_info.loc[event.id, :] = None
        # Default to site_idx of 0 if no site given
        if event.parent:
            site_idx = event.parent.id
        else:
            site_idx = 0
        self.event_info.set_value(event.id, 'site_idx', site_idx)
        # Event and cluster information
        self.event_info.set_value(event.id, '1-BDC', round(1.0 - event.info.estimated_bdc, 2))
        self.event_info.set_value(event.id, 'z_peak', round(event.cluster.max, 2))
        self.event_info.set_value(event.id, 'z_mean', round(event.cluster.mean, 2))
        self.event_info.set_value(event.id, 'cluster_size', event.cluster.size)
        #        self.tables.event_info.set_value(event.id, ['refx','refy','refz'], list(self.grid.grid2cart([event.cluster.peak],origin_shift=False)[0]))
        self.event_info.set_value(event.id, ['x', 'y', 'z'], list(dataset.model.alignment.ref2nat(
            coordinates=grid.grid2cart([event.cluster.peak], origin_shift=True))[0]))
        self.event_info.set_value(event.id, 'global_correlation_to_average_map',
                                         event.info.global_correlation)
        self.event_info.set_value(event.id, 'local_correlation_to_average_map', event.info.local_correlation)

    def write_event_table(self, table_path, event_table):
        # Write the event data only once events have been recorded
        if len(event_table):
            # Sort the event data by z-peak and write out
            sort_eve = event_table.sort_values(by=['dtag', "event_idx"],
                                               ascending=[1, 1],
                                               )
            print(sort_eve)

            # TODO: not sure if removing comb_tab dependancy is reasonable
            # sort_eve = sort_eve.join(comb_tab, how='inner')

            sort_eve.to_csv(path_or_buf=table_path)


