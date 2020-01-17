import os, sys, configparser
from collections import OrderedDict

from multi_dataset_crystalography.dataset.dataloader import DefaultDataloader

from multi_dataset_crystalography.dataset.sample_loader import PanDDAGridSetup
from multi_dataset_crystalography.dataset.sample_loader import (DefaultSampleLoader,
                                                                PanddaReferenceMapLoader,
                                                                )
from multi_dataset_crystalography.dataset.sample_loader import MapLoaderDask

from multi_dataset_crystalography.dataset.transforms import PanddaDataChecker
from multi_dataset_crystalography.dataset.transforms import PanddaDiffractionScaler
from multi_dataset_crystalography.dataset.transforms import PanddaDatasetFilterer
from multi_dataset_crystalography.dataset.transforms import PanddaDatasetFiltererWilsonRMSD
from multi_dataset_crystalography.dataset.transforms import PanddaDefaultStructureAligner

from multi_dataset_crystalography.dataset.reference import DefaultReferenceGetter

from multi_dataset_crystalography.dataset.sample_loader import PanddaDiffractionDataTruncater


from pandda_2 import (statistical_model,
                      cluster_outliers,
                      filter_clusters,
                      event_analyser,
                      )


from multi_dataset_crystalography.dataset.partitioner import DefaultPanDDAPartitionsGetter
#
# from criticise import PanDDADefaultMapMaker, PanDDADefaultEventTableShell
#
# from classes import Criticiser, CriticiserAll


class Options:
    def __init__(self, config):
        # ============================================================================>
        # Dataset
        # ============================================================================>
        # # Dataloader
        self.dataloader = DefaultDataloader(config.input.data_dirs,
                                            config.input.pdb_style,
                                            config.input.mtz_style,
                                            config.input.regex.pdb_regex,
                                            config.input.regex.mtz_regex,
                                            config.input.regex.dir_regex,
                                            config.input.flags.only_datasets,
                                            config.input.flags.ignore_datasets,
                                            config.output.dataset_prefix,
                                            config.output.out_dir,
                                            config.input.lig_style
                                            )

        # # Reference
        self.get_reference = DefaultReferenceGetter(out_dir=config.output.out_dir,
                                                    reference_pdb_path=None,
                                                    reference_mtz_path=None,
                                                    reference_structure_factors=None,
                                                    structure_factors=config.params.diffraction_data.structure_factors
                                                    )

        # # Output
        # self.pandda_output = PanddaOutputSetup(config.output.out_dir,
        #                                        config.input.lig_style
        #                                        )

        # # Transforms
        self.transform_data_check = PanddaDataChecker(config.params.diffraction_data.structure_factors,
                                                      config.params.diffraction_data.checks.low_resolution_completeness,
                                                      config.params.diffraction_data.checks.all_data_are_valid_values
                                                      )

        self.transform_scale_diffraction = PanddaDiffractionScaler(
            config.params.diffraction_data.apply_b_factor_scaling
        )

        self.transform_filter_structure = PanddaDatasetFilterer(
            config.params.filtering.flags.similar_models_only,
            config.params.filtering.max_rfree,
            config.params.filtering.flags.same_space_group_only
        )

        self.transform_filter_wilson = PanddaDatasetFiltererWilsonRMSD(
            config.params.excluding.max_wilson_plot_z_score,
            config.params.diffraction_data.apply_b_factor_scaling,
            config.params.diffraction_data.structure_factors
        )

        self.transform_align = PanddaDefaultStructureAligner(
            method=config.params.alignment.method,
            cpus=config.settings.cpus
        )

        # # Sample loader
        self.grid_loader = PanDDAGridSetup(cpus=int(config.settings.cpus),
                                           mask_pdb=config.params.masks.pdb,
                                           align_mask_to_reference=bool(
                                               config.params.masks.align_mask_to_reference),
                                           alignment_method=config.params.alignment.method,
                                           outer_mask=float(config.params.masks.outer_mask),
                                           inner_mask=float(config.params.masks.inner_mask),
                                           inner_mask_symmetry=float(config.params.masks.inner_mask_symmetry),
                                           grid_spacing=float(config.params.maps.grid_spacing),
                                           padding=float(config.params.maps.padding),
                                           verbose=bool(config.settings.verbose),
                                           mask_selection_string=(None
                                                                  if str(
                                               config.params.z_map_anlysis.masks.selection_string) == "None"
                                                                  else str(
                                               config.params.z_map_anlysis.masks.selection_string)
                                                                  )
                                           )

        self.sample_loader = DefaultSampleLoader(config.params.maps.resolution_factor,
                                                 config.params.maps.density_scaling,
                                                 int(config.settings.cpus),
                                                 bool(config.settings.verbose),
                                                 self.grid_loader
                                                 )

        # # Refernce map loader
        self.reference_map_getter = PanddaReferenceMapLoader(config.params.maps.resolution_factor,
                                                             config.params.maps.density_scaling,
                                                             )

        # # Partititioner
        test = str(config.input.flags.test).split(",") if (str(config.input.flags.test) != "None") else None
        train = str(config.input.flags.train).split(",") if (str(config.input.flags.train) != "None") else None
        not_test = str(config.input.flags.not_test).split(",") if (
                str(config.input.flags.not_test) != "None") else None
        not_train = str(config.input.flags.not_train).split(",") if (
                str(config.input.flags.not_train) != "None") else None

        self.partitioner = DefaultPanDDAPartitionsGetter(test=test,
                                                         train=train,
                                                         not_test=not_test,
                                                         not_train=not_train,
                                                         )

        # ============================================================================>
        # Event Model
        # ============================================================================>
        self.statistical_model = statistical_model.PanDDANormalModel(method="adjusted+uncertainty",
                                                                     cpus=config.settings.cpus,
                                                                     )

        self.clusterer = cluster_outliers.PanDDADefaultCluster(
            grid_clustering_cutoff=float(
                config.params.z_map_anlysis.agglomerative_hierarchical.clustering_cutoff),
            negative_values=bool(config.params.z_map_anlysis.negative_values),
            cluster_method=config.params.z_map_anlysis.clustering_method,
            contour_level=float(config.params.z_map_anlysis.contour_level),
            outer_mask=float(config.params.masks.outer_mask),
            inner_mask_symmetry=float(config.params.masks.inner_mask_symmetry)
        )

        self.event_finder = filter_clusters.PanDDADefaultEventFinder(
            min_blob_z_peak=float(config.params.z_map_anlysis.min_blob_z_peak),
            grid_minimum_volume=float(
                config.params.z_map_anlysis.min_blob_volume),
            grid_spacing=float(config.params.maps.grid_spacing)
        )

        self.diffraction_data_truncator = PanddaDiffractionDataTruncater()

        #
        # self.bdc_calculator = PanDDADefaultBDCCalculator(config.params.background_correction.max_bdc,
        #                                                  config.params.background_correction.min_bdc,
        #                                                  config.params.background_correction.increment,
        #                                                  config.params.background_correction.output_multiplier
        #                                                  )
        #
        self.event_analyser = event_analyser.EventAnalyser(max_bdc=config.params.background_correction.max_bdc,
                                                           min_bdc=config.params.background_correction.min_bdc,
                                                           increment=config.params.background_correction.increment,
                                                           output_multiplier=config.params.background_correction.output_multiplier,
                                                           )
        #
        # self.map_maker = PanDDADefaultMapMaker()
        #
        # self.event_table_maker = PanDDADefaultEventTableShell(order_by=config.results.events.order_by)

        # self.event_model = PanDDAEventModel(self.statistical_model,
        #                                     self.clusterer,
        #                                     self.event_finder,
        #                                     statistics=[])

        # ############################
        # Criticism
        # ############################
        self.map_loader = MapLoaderDask(config.settings.verbose,
                                        config.params.maps.resolution_factor,
                                        config.params.maps.density_scaling,
                                        )
