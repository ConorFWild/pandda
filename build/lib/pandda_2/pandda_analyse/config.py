import os, sys, configparser
from collections import OrderedDict

module_path = os.path.abspath(os.path.join('..'))
if module_path not in sys.path:
    sys.path.append(module_path)
b = sys.path
sys.path = [module_path] + b

from multi_dataset_crystalography.dataset.dataloader import DefaultDataloader

from multi_dataset_crystalography.dataset.sample_loader import PanDDAGridSetup
from multi_dataset_crystalography.dataset.sample_loader import (DefaultSampleLoader,
                                                                PanddaReferenceMapLoader,
                                                                )
from multi_dataset_crystalography.dataset.sample_loader import MapLoaderDask

from pandda_output import PanddaOutputSetup

from event_model import (PanDDAEventModel,
                         PanDDANormalModel,
                         PanDDADefaultCluster,
                         PanDDADefaultEventFinder,
                         PanDDADefaultBDCCalculator,
                         analyse_dataset_events,
                         )

from multi_dataset_crystalography.dataset.transforms import PanddaDataChecker
from multi_dataset_crystalography.dataset.transforms import PanddaDiffractionScaler
from multi_dataset_crystalography.dataset.transforms import PanddaDatasetFilterer
from multi_dataset_crystalography.dataset.transforms import PanddaDatasetFiltererWilsonRMSD
from multi_dataset_crystalography.dataset.transforms import PanddaDefaultStructureAligner

from multi_dataset_crystalography.dataset.reference import DefaultReferenceGetter

from multi_dataset_crystalography.dataset.partitioner import DefaultPanDDAPartitionsGetter

from criticise import PanDDADefaultMapMaker, PanDDADefaultEventTableShell

from classes import Criticiser, CriticiserAll


# from event_model import ClusterHierarchical


########################################################################################################################
class Input:

    def __init__(self, config_obj):
        self.data_dirs = config_obj.data_dirs
        self.pdb_style = config_obj.pdb_style
        self.mtz_style = config_obj.mtz_style
        self.lig_style = config_obj.lig_style
        self.regex = Regex(config_obj.regex)
        self.flags = InputFlags(config_obj.flags)


class Regex:

    def __init__(self, config_obj):
        self.pdb_regex = config_obj.pdb_regex
        self.mtz_regex = config_obj.mtz_regex
        self.dir_regex = config_obj.dir_regex


class InputFlags:
    def __init__(self, config_obj):
        self.only_datasets = config_obj.only_datasets
        self.ignore_datasets = config_obj.ignore_datasets
        self.test = config_obj.test
        self.train = config_obj.train
        self.not_test = config_obj.not_test
        self.not_train = config_obj.not_train


########################################################################################################################
class Output:

    def __init__(self, config_obj):
        self.out_dir = config_obj.out_dir
        self.dataset_prefix = config_obj.dataset_prefix


########################################################################################################################
class Params:

    def __init__(self, config_obj):
        self.diffraction_data = DiffractionData(config_obj.diffraction_data)
        self.filtering = Filtering(config_obj.filtering)
        self.maps = Maps(config_obj.maps)
        self.alignment = Alignment(config_obj.alignment)
        self.masks = Masks(config_obj.masks)
        self.excluding = Excluding(config_obj.excluding)
        self.z_map_anlysis = ZMapAnalysis(config_obj.z_map_analysis)
        self.background_correction = BackgroundCorrection(config_obj.background_correction)


class DiffractionData:

    def __init__(self, config_obj):
        self.structure_factors = config_obj.structure_factors
        self.apply_b_factor_scaling = config_obj.apply_b_factor_scaling

        self.checks = DiffractionDataChecks(config_obj.checks)


class DiffractionDataChecks:
    def __init__(self, config_obj):
        self.low_resolution_completeness = config_obj.low_resolution_completeness
        self.all_data_are_valid_values = config_obj.all_data_are_valid_values


class Filtering:

    def __init__(self, config_obj):
        self.max_rfree = config_obj.max_rfree
        self.flags = FilteringFlags(config_obj.flags)


class FilteringFlags:

    def __init__(self, config_obj):
        self.same_space_group_only = config_obj.same_space_group_only
        self.similar_models_only = config_obj.similar_models_only


class Maps:

    def __init__(self, config_obj):
        self.resolution_factor = config_obj.resolution_factor
        self.grid_spacing = config_obj.grid_spacing
        self.density_scaling = config_obj.density_scaling
        self.padding = config_obj.padding


class Alignment:

    def __init__(self, config_obj):
        self.method = config_obj.method


class Masks:

    def __init__(self, config_obj):
        self.pdb = config_obj.pdb
        self.align_mask_to_reference = config_obj.align_mask_to_reference
        self.outer_mask = config_obj.outer_mask
        self.inner_mask = config_obj.inner_mask
        self.inner_mask_symmetry = config_obj.inner_mask_symmetry


class Excluding:
    def __init__(self, config_obj):
        self.max_wilson_plot_z_score = config_obj.max_wilson_plot_z_score


class ZMapAnalysis:
    def __init__(self, config_obj):
        self.clustering_method = config_obj.clustering_method
        self.contour_level = config_obj.contour_level
        self.negative_values = config_obj.negative_values

        self.min_blob_volume = config_obj.min_blob_volume
        self.min_blob_z_peak = config_obj.min_blob_z_peak
        self.masks = ZMasks(config_obj.masks)
        self.agglomerative_hierarchical = AgglomerativeHierarchical(config_obj.agglomerative_hierarchical)


class ZMasks:
    def __init__(self, config_obj):
        self.selection_string = config_obj.selection_string
        self.outer_mask = config_obj.outer_mask
        self.inner_mask = config_obj.inner_mask


class BackgroundCorrection:
    def __init__(self, config_obj):
        self.max_bdc = config_obj.max_bdc
        self.min_bdc = config_obj.min_bdc
        self.increment = config_obj.increment
        self.output_multiplier = config_obj.output_multiplier


class AgglomerativeHierarchical:
    def __init__(self, config_obj):
        self.clustering_cutoff = config_obj.clustering_cutoff


########################################################################################################################
class Results:
    def __init__(self, config_obj):
        self.events = Events(config_obj.events)


class Events:
    def __init__(self, config_obj):
        self.order_by = config_obj.order_by


########################################################################################################################
class Processing:
    def __init__(self, config_obj):
        self.backend = config_obj.backend


########################################################################################################################
class Settings:
    def __init__(self, config_obj):
        self.cpus = config_obj.cpus
        self.verbose = config_obj.verbose


########################################################################################################################
class Config:

    def __init__(self, config_obj):
        self.settings = Settings(config_obj.settings)

        self.input = Input(config_obj.pandda.input)
        self.output = Output(config_obj.pandda.output)
        self.params = Params(config_obj.pandda.params)
        self.processing = Processing(config_obj.pandda.processing)
        self.results = Results(config_obj.pandda.results)

    # def replace_param_phil(self, param_obj, phil_obj):
    #     # loop over param obj's atributes, try to match to phil attributes
    #     for param_name, param_val in param_obj.__dict__.items():
    #
    #         if param_name[0] == "_":
    #             continue
    #
    #         # Check whether a further level of nesting, replace if not
    #         if param_val is None:
    #             param_obj.__dict__[param_name] = phil_obj.__dict__[param_name]
    #
    #         # Recurse if there is a further level of nesting
    #         else:
    #             param_obj.__dict__[param_name] = self.replace_param_phil(param_val,
    #                                                                      phil_obj.__dict__[param_name]
    #                                                                      )
    #
    #     return param_obj
    #
    #




class PanDDAConfig:

    def __init__(self, phil):
        # ============================================================================>
        # Transform Variables
        # ============================================================================>
        config = Config(phil)
        self.config = config

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
        self.pandda_output = PanddaOutputSetup(config.output.out_dir,
                                               config.input.lig_style
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
                                           align_mask_to_reference=bool(config.params.masks.align_mask_to_reference),
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
        not_test = str(config.input.flags.not_test).split(",") if (str(config.input.flags.not_test) != "None") else None
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
        self.statistical_model = PanDDANormalModel(method="adjusted+uncertainty", cpus=1)

        self.clusterer = PanDDADefaultCluster(
            grid_clustering_cutoff=float(config.params.z_map_anlysis.agglomerative_hierarchical.clustering_cutoff),
            negative_values=bool(config.params.z_map_anlysis.negative_values),
            cluster_method=config.params.z_map_anlysis.clustering_method,
            contour_level=float(config.params.z_map_anlysis.contour_level),
            outer_mask=float(config.params.masks.outer_mask),
            inner_mask_symmetry=float(config.params.masks.inner_mask_symmetry)
            )

        self.event_finder = PanDDADefaultEventFinder(min_blob_z_peak=float(config.params.z_map_anlysis.min_blob_z_peak),
                                                     grid_minimum_volume=float(
                                                         config.params.z_map_anlysis.min_blob_volume),
                                                     grid_spacing=float(config.params.maps.grid_spacing)
                                                     )

        self.bdc_calculator = PanDDADefaultBDCCalculator(config.params.background_correction.max_bdc,
                                                         config.params.background_correction.min_bdc,
                                                         config.params.background_correction.increment,
                                                         config.params.background_correction.output_multiplier
                                                         )

        self.event_analyser = lambda dataset, dataset_map, ref_map, events, grid: analyse_dataset_events(dataset,
                                                                                                         dataset_map,
                                                                                                         ref_map,
                                                                                                         events,
                                                                                                         grid,
                                                                                                         max_bdc=config.params.background_correction.max_bdc,
                                                                                                         min_bdc=config.params.background_correction.min_bdc,
                                                                                                         increment=config.params.background_correction.increment,
                                                                                                         output_multiplier=config.params.background_correction.output_multiplier,
                                                                                                         )

        self.map_maker = PanDDADefaultMapMaker()

        self.event_table_maker = PanDDADefaultEventTableShell(order_by=config.results.events.order_by)

        self.event_model = PanDDAEventModel(self.statistical_model,
                                            self.clusterer,
                                            self.event_finder,
                                            statistics=[])

        # ############################
        # Criticism
        # ############################
        self.map_loader = MapLoaderDask(config.settings.verbose,
                                        config.params.maps.resolution_factor,
                                        config.params.maps.density_scaling,
                                        )
