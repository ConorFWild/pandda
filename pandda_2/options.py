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
                      partitioner,
                      )

from pandda_2 import statistical_model

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

        # # Refernce map loader
        self.reference_map_getter = PanddaReferenceMapLoader(config.params.maps.resolution_factor,
                                                             config.params.maps.density_scaling,
                                                             )

        # # Partititioner
        # Order:
        # ignore datsets (remove from test and train)
        # only datasets (remove all others from test and train)
        # exclude datasets (remove from test)
        # ground state (remove all others from test)
        ignore = str(config.input.flags.ignore_datasets).split(",") if (
                    str(config.input.flags.ignore_datasets) != "None") else []
        only = str(config.input.flags.only_datasets).split(",") if (
                    str(config.input.flags.only_datasets) != "None") else []
        exclude = str(config.input.flags.exclude_from_characterisation).split(",") if (
                    str(config.input.flags.exclude_from_characterisation) != "None") else []
        ground_state_datasets = str(config.input.flags.ground_state_datasets).split(",") if (
                    str(config.input.flags.ground_state_datasets) != "None") else []

        not_test = set(ignore)
        not_train = set(ignore).union(set(exclude))

        test = set(only).difference(not_train)
        train = set(ground_state_datasets).intersection(set(only)).difference(not_test)

        self.partitioner = partitioner.DefaultPanDDAPartitionsGetter(test=test,
                                                                     train=train,
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

        self.event_analyser = event_analyser.EventAnalyser(max_bdc=config.params.background_correction.max_bdc,
                                                           min_bdc=config.params.background_correction.min_bdc,
                                                           increment=config.params.background_correction.increment,
                                                           output_multiplier=config.params.background_correction.output_multiplier,
                                                           )

        # ############################
        # Criticism
        # ############################
        self.map_loader = MapLoaderDask(config.settings.verbose,
                                        config.params.maps.resolution_factor,
                                        config.params.maps.density_scaling,
                                        )

        # Get dataset loader
        self.load_dataset = load_dataset.LoadDataset(dataloader=self.dataloader,
                                                     )

        # Get dataset transformer
        self.transform_dataset = transform_dataset.TransformDataset(transform_data_check=self.transform_data_check,
                                                                    transform_scale_diffraction=self.transform_scale_diffraction,
                                                                    transform_filter_structure=self.transform_filter_structure,
                                                                    transform_filter_wilson=self.transform_filter_wilson,
                                                                    transform_align=self.transform_align,
                                                                    )

        # Get grid loader
        self.get_grid = self.grid_loader

        # Get output handler
        self.define_tree = define_tree.DefineTree(output_dir=config.output.out_dir)
        self.make_tree = make_tree.MakeTree(overwrite=config.output.overwrite)
        self.copy_dataset_files = copy_dataset_files.DatasetFileCopier()

        self.output = output.Output(define_tree=self.define_tree,
                                    make_tree=self.make_tree,
                                    copy_dataset_files=self.copy_dataset_files,
                                    )

        # Get resolution shell scheme
        self.create_shells = create_shells.CreateShells(
            min_train_datasets=config.params.statistical_maps.min_build_datasets,
            max_test_datasets=config.params.statistical_maps.max_build_datasets,
            cutoff=0.1,
            )

        # Get Resolution shell processor
        self.n_cpus_shells = config.processing.process_dict_n_cpus
        self.h_vmem = config.processing.h_vmem
        self.m_mem_free = config.processing.m_mem_free

        if config.processing.process_dict == "joblib":
            process_in_shell = processor.ProcessorDictJoblib(config.processing.process_dict_n_cpus)

        elif config.processing.process_dict == "seriel":
            process_in_shell = processor.ProcessorDict()

        self.process_shell = process_shell.ProcessShell(diffraction_data_truncator=self.diffraction_data_truncator,
                                                        reference_map_getter=self.reference_map_getter,
                                                        get_reference_map=get_reference_map.GetReferenceMap(),
                                                        map_loader=self.map_loader,
                                                        statistical_model=self.statistical_model,
                                                        fit_model=fit_model.FitModel(),
                                                        evaluate_model=evaluate_model.EvaluateModel(),
                                                        cluster_outliers=self.clusterer,
                                                        filter_clusters=self.event_finder,
                                                        event_analyser=self.event_analyser,
                                                        make_event_map=make_event_map.MakeEventMap(),
                                                        make_z_map=make_z_map.MakeZMap(),
                                                        make_mean_map=make_mean_map.MakeMeanMap(),
                                                        make_event_table=make_event_table.MakeEventTable(),
                                                        process=process_in_shell,
                                                        )

        if config.processing.process_shells == "luigi":
            self.processer = processor.ProcessorLuigi(jobs=10,
                                                      parallel_env="smp",
                                                      run_locally=False,
                                                      n_cpu=self.n_cpus_shells,
                                                      h_vmem=self.h_vmem,
                                                      m_mem_free=self.m_mem_free,
                                                      )
        elif config.processing.process_shells == "seriel":
            self.processer = processor.Processor()

        # Get site table creator
        self.create_sites_table = create_sites_table.CreateSitesTable()

        # Get event table outputter
        self.output_sites_table = output_sites_table.OutputSitesTable()

        # Get event table processor
        self.create_event_table = create_event_table.CreateEventTable()

        # Get event table outputter
        self.output_event_table = output_event_table.OutputEventTable()
