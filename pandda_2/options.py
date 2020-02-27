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

# from pandda_2 import (statistical_model,
#                       cluster_outliers,
#                       filter_clusters,
#                       event_analyser,
#                       )

from pandda_2 import (config,
                      pandda_phil,
                      load_dataset,
                      transform_dataset,
                      get_dataset,
                      get_reference,
                      get_grid,
                      define_tree,
                      make_tree,
                      copy_dataset_files,
                      output,
                      create_shells,
                      map_loader,
                      get_reference_map,
                      statistical_model,
                      fit_model,
                      evaluate_model,
                      cluster_outliers,
                      filter_clusters,
                      event_analyser,
                      make_event_map,
                      make_z_map,
                      make_mean_map,
                      make_event_table,
                      process_shell,
                      process_shells,
                      processor,
                      create_sites_table,
                      output_sites_table,
                      create_event_table,
                      output_event_table,
                      standard_pandda,
                      autobuild,
                      checks,
                      )

from pandda_2 import statistical_model

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

        # Get dataset loader
        self.load_dataset = load_dataset.LoadDataset(dataloader=self.dataloader,
                                                     sample_loader=self.sample_loader,
                                                     )

        # Get reference loader
        self.get_reference = self.get_reference

        # Get dataset transformer
        self.transform_dataset = transform_dataset.TransformDataset(transform_data_check=self.transform_data_check,
                                                                    transform_scale_diffraction=self.transform_scale_diffraction,
                                                                    transform_filter_structure=self.transform_filter_structure,
                                                                    transform_filter_wilson=self.transform_filter_wilson,
                                                                    transform_align=self.transform_align,
                                                                    )

        # Get grid loader
        self.get_grid = self.grid_loader

        # Get partitioner
        self.partitioner = self.partitioner

        # Get output handler
        self.define_tree = define_tree.DefineTree(output_dir=config.output.out_dir)
        self.make_tree = make_tree.MakeTree(overwrite=config.output.overwrite)
        self.copy_dataset_files = copy_dataset_files.DatasetFileCopier()

        self.output = output.Output(define_tree=self.define_tree,
                                    make_tree=self.make_tree,
                                    copy_dataset_files=self.copy_dataset_files,
                                    )

        # Get resolution shell scheme
        self.create_shells = create_shells.CreateShells(min_train_datasets=60,
                                                        max_test_datasets=60,
                                                        cutoff=0.1,
                                                        )

        # Get Resolution shell processor
        diffraction_data_truncator_obj = self.diffraction_data_truncator
        reference_map_getter_obj = self.reference_map_getter
        get_reference_map_obj = get_reference_map.GetReferenceMap()

        map_loader_obj = self.map_loader
        statistical_model_obj = self.statistical_model
        fit_model_obj = fit_model.FitModel()
        evaluate_model_obj = evaluate_model.EvaluateModel()
        cluster_outliers_obj = self.clusterer
        filter_clusters_obj = self.event_finder
        event_analyser_obj = self.event_analyser

        make_event_map_obj = make_event_map.MakeEventMap()
        make_z_map_obj = make_z_map.MakeZMap()
        make_mean_map_obj = make_mean_map.MakeMeanMap()
        make_event_table_obj = make_event_table.MakeEventTable()

        process_in_shell = processor.ProcessorDictJoblib()

        self.process_shell = process_shell.ProcessShell(diffraction_data_truncator=diffraction_data_truncator_obj,
                                                        reference_map_getter=reference_map_getter_obj,
                                                        get_reference_map=get_reference_map_obj,
                                                        map_loader=map_loader_obj,
                                                        statistical_model=statistical_model_obj,
                                                        fit_model=fit_model_obj,
                                                        evaluate_model=evaluate_model_obj,
                                                        cluster_outliers=cluster_outliers_obj,
                                                        filter_clusters=filter_clusters_obj,
                                                        event_analyser=event_analyser_obj,
                                                        make_event_map=make_event_map_obj,
                                                        make_z_map=make_z_map_obj,
                                                        make_mean_map=make_mean_map_obj,
                                                        make_event_table=make_event_table_obj,
                                                        process=process_in_shell,
                                                        )

        # self.processer = processor.ProcessorLuigi(jobs=10,
        #                                           parallel_env="smp",
        #                                           n_cpu=12,
        #                                           run_locally=False,
        #                                           )
        self.processer = processor.Processor()

        # Get site table creator
        self.create_sites_table = create_sites_table.CreateSitesTable()

        # Get event table outputter
        self.output_sites_table = output_sites_table.OutputSitesTable()

        # Get event table processor
        self.create_event_table = create_event_table.CreateEventTable()

        # Get event table outputter
        self.output_event_table = output_event_table.OutputEventTable()

        # Autobuilders
        self.autobuilder = autobuild.AutobuildQFit()
