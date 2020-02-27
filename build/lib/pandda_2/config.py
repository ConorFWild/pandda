import os, sys, configparser, copy
from collections import OrderedDict

import libtbx.phil
from libtbx.utils import Sorry

try:
    VERSION = pkg_resources.get_distribution("panddas").version
except:
    VERSION = '(developer -- see setup.py file)'

HEADER_TEXT = """
------------------------------------------------------------------>
-         .__               __
-    ____ |__|____    _____/  |_
-   / ___\\|  \\__  \\  /    \\   __\\
-  / /_/  >  |/ __ \\|   |  \\  |
-  \\___  /|__(____  /___|  /__|     A crystallographic toolbox for
-  /_____/         \\/     \\/        protein structure determination
-
------------------------------------------------------------------>

> {program}
{description}
------------------------------------------------------------------>
"""

class module_info:
    name        = 'giant'
    version     = VERSION
    header_text = HEADER_TEXT


def extract_params_default(master_phil, args, blank_arg_prepend=None, home_scope=None, module_info=None):
    """Extract the parameters by a default script"""
    show_version_and_exit_maybe(module_info, args)
    show_defaults_and_exit_maybe(master_phil, args)
    working_phil = parse_phil_args(master_phil, args, blank_arg_prepend=blank_arg_prepend, home_scope=home_scope)
    return working_phil


def show_version_and_exit_maybe(module_info, args):
    if '--version' not in args: return
    if module_info is None:
        print 'no version information available'
    else:
        print '{} version: {}'.format(module_info.name, module_info.version)
    sys.exit()


def show_defaults_and_exit_maybe(master_phil, args):
    """Show master_phil and exit if requested"""

    attributes_level = expert_level = 0

    if '-h' in args: attributes_level = 1
    if '-hh' in args: attributes_level = 2
    if '-hhh' in args: attributes_level = 3
    if '-a' in args: expert_level = 1
    if '-aa' in args: expert_level = 2
    if '-aaa' in args: expert_level = 3

    if ('?' in args) or ('--show-defaults' in args) or (not attributes_level==expert_level==0):
        print('\n====================== Showing Default Parameters =====================\n')
        master_phil.show(expert_level=expert_level, attributes_level=attributes_level)
    else:
        return

    raise SystemExit('\n============================= Now Exiting =============================\n')


def parse_phil_args(master_phil, args, blank_arg_prepend=None, home_scope=None):

    if blank_arg_prepend is None:
        pass
    elif isinstance(blank_arg_prepend, dict):
        for item in blank_arg_prepend.values():
            assert '=' in item
    elif isinstance(blank_arg_prepend, str):
        assert '=' in blank_arg_prepend
    else:
        raise Exception('blank_arg_prepend must be str or dict')

    # Copy the args so that we can remove items from the list without affecting args etc
    args = copy.copy(args)
    # Construct interpreter
    cmd_interpr = master_phil.command_line_argument_interpreter(home_scope=home_scope)

    # Process any args that are eff files
    eff_files = [f for f in args if ((f.endswith('.eff') or f.endswith('.def')) and (not f.count('=')) and os.path.isfile(f))]
    # Remove them from the original lists
    [args.remove(f) for f in eff_files]
    # Parse the 'eff' files - these should contain phils
    #eff_sources = [libtbx.phil.parse(open(f, 'r').read()) for f in eff_files]
    eff_sources = [cmd_interpr.process(open(f, 'r').read()) for f in eff_files]

    # Process input arguments
    arg_sources = []
    for arg in args:
        try:
            # Prepend if blank
            if '=' not in arg:
                if isinstance(blank_arg_prepend, dict):
                    found_key = False
                    for key in blank_arg_prepend.keys():
                        if key is None:
                            continue
                        if arg.endswith(key):
                            arg = blank_arg_prepend[key]+arg
                            found_key = True
                            break
                    if (found_key == False) and (None in blank_arg_prepend.keys()):
                        arg = blank_arg_prepend[None]+arg
                elif isinstance(blank_arg_prepend, str):
                    arg = blank_arg_prepend+arg
            # Attempt to process arg
            cmd_line_args = cmd_interpr.process(arg=arg)
        except KeyboardInterrupt:
            raise
        except Exception:
            raise Sorry("Unknown file or keyword: %s" % arg)
        else:
            arg_sources.append(cmd_line_args)
    # Extract Scope object (putting eff sources first so that they're overridden if double-defined)
    working_phil = master_phil.fetch(sources=eff_sources+arg_sources)

    return working_phil



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
        self.overwrite = config_obj.overwrite


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

class Autobuilding:
    def __init__(self, config_obj):
        self.autobuild = config_obj.autobuild

########################################################################################################################
class Config:

    def __init__(self, config_obj):
        self.settings = Settings(config_obj.settings)

        self.input = Input(config_obj.pandda.input)
        self.output = Output(config_obj.pandda.output)
        self.params = Params(config_obj.pandda.params)
        self.processing = Processing(config_obj.pandda.processing)
        self.results = Results(config_obj.pandda.results)
        self.autobuilding = Autobuilding(config_obj.pandda.autobuilding)

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
#
#
# class PanDDAConfig:
#
#     def __init__(self, phil):
#         # ============================================================================>
#         # Transform Variables
#         # ============================================================================>
#         config = Config(phil)
#         self.config = config
#
#         # ============================================================================>
#         # Dataset
#         # ============================================================================>
#         # # Dataloader
#         self.dataloader = DefaultDataloader(config.input.data_dirs,
#                                             config.input.pdb_style,
#                                             config.input.mtz_style,
#                                             config.input.regex.pdb_regex,
#                                             config.input.regex.mtz_regex,
#                                             config.input.regex.dir_regex,
#                                             config.input.flags.only_datasets,
#                                             config.input.flags.ignore_datasets,
#                                             config.output.dataset_prefix,
#                                             config.output.out_dir,
#                                             config.input.lig_style
#                                             )
#
#         # # Reference
#         self.get_reference = DefaultReferenceGetter(out_dir=config.output.out_dir,
#                                                     reference_pdb_path=None,
#                                                     reference_mtz_path=None,
#                                                     reference_structure_factors=None,
#                                                     structure_factors=config.params.diffraction_data.structure_factors
#                                                     )
#
#         # # Output
#         self.pandda_output = PanddaOutputSetup(config.output.out_dir,
#                                                config.input.lig_style
#                                                )
#
#         # # Transforms
#         self.transform_data_check = PanddaDataChecker(config.params.diffraction_data.structure_factors,
#                                                       config.params.diffraction_data.checks.low_resolution_completeness,
#                                                       config.params.diffraction_data.checks.all_data_are_valid_values
#                                                       )
#
#         self.transform_scale_diffraction = PanddaDiffractionScaler(
#             config.params.diffraction_data.apply_b_factor_scaling
#         )
#
#         self.transform_filter_structure = PanddaDatasetFilterer(
#             config.params.filtering.flags.similar_models_only,
#             config.params.filtering.max_rfree,
#             config.params.filtering.flags.same_space_group_only
#         )
#
#         self.transform_filter_wilson = PanddaDatasetFiltererWilsonRMSD(
#             config.params.excluding.max_wilson_plot_z_score,
#             config.params.diffraction_data.apply_b_factor_scaling,
#             config.params.diffraction_data.structure_factors
#         )
#
#         self.transform_align = PanddaDefaultStructureAligner(
#             method=config.params.alignment.method,
#             cpus=config.settings.cpus
#         )
#
#         # # Sample loader
#         self.grid_loader = PanDDAGridSetup(cpus=int(config.settings.cpus),
#                                            mask_pdb=config.params.masks.pdb,
#                                            align_mask_to_reference=bool(config.params.masks.align_mask_to_reference),
#                                            alignment_method=config.params.alignment.method,
#                                            outer_mask=float(config.params.masks.outer_mask),
#                                            inner_mask=float(config.params.masks.inner_mask),
#                                            inner_mask_symmetry=float(config.params.masks.inner_mask_symmetry),
#                                            grid_spacing=float(config.params.maps.grid_spacing),
#                                            padding=float(config.params.maps.padding),
#                                            verbose=bool(config.settings.verbose),
#                                            mask_selection_string=(None
#                                                                   if str(
#                                                config.params.z_map_anlysis.masks.selection_string) == "None"
#                                                                   else str(
#                                                config.params.z_map_anlysis.masks.selection_string)
#                                                                   )
#                                            )
#
#         self.sample_loader = DefaultSampleLoader(config.params.maps.resolution_factor,
#                                                  config.params.maps.density_scaling,
#                                                  int(config.settings.cpus),
#                                                  bool(config.settings.verbose),
#                                                  self.grid_loader
#                                                  )
#
#         # # Refernce map loader
#         self.reference_map_getter = PanddaReferenceMapLoader(config.params.maps.resolution_factor,
#                                                              config.params.maps.density_scaling,
#                                                              )
#
#         # # Partititioner
#         test = str(config.input.flags.test).split(",") if (str(config.input.flags.test) != "None") else None
#         train = str(config.input.flags.train).split(",") if (str(config.input.flags.train) != "None") else None
#         not_test = str(config.input.flags.not_test).split(",") if (str(config.input.flags.not_test) != "None") else None
#         not_train = str(config.input.flags.not_train).split(",") if (
#                     str(config.input.flags.not_train) != "None") else None
#
#         self.partitioner = DefaultPanDDAPartitionsGetter(test=test,
#                                                          train=train,
#                                                          not_test=not_test,
#                                                          not_train=not_train,
#                                                          )
#
#         # ============================================================================>
#         # Event Model
#         # ============================================================================>
#         self.statistical_model = PanDDANormalModel(method="adjusted+uncertainty", cpus=1)
#
#         self.clusterer = PanDDADefaultCluster(
#             grid_clustering_cutoff=float(config.params.z_map_anlysis.agglomerative_hierarchical.clustering_cutoff),
#             negative_values=bool(config.params.z_map_anlysis.negative_values),
#             cluster_method=config.params.z_map_anlysis.clustering_method,
#             contour_level=float(config.params.z_map_anlysis.contour_level),
#             outer_mask=float(config.params.masks.outer_mask),
#             inner_mask_symmetry=float(config.params.masks.inner_mask_symmetry)
#             )
#
#         self.event_finder = PanDDADefaultEventFinder(min_blob_z_peak=float(config.params.z_map_anlysis.min_blob_z_peak),
#                                                      grid_minimum_volume=float(
#                                                          config.params.z_map_anlysis.min_blob_volume),
#                                                      grid_spacing=float(config.params.maps.grid_spacing)
#                                                      )
#
#         self.bdc_calculator = PanDDADefaultBDCCalculator(config.params.background_correction.max_bdc,
#                                                          config.params.background_correction.min_bdc,
#                                                          config.params.background_correction.increment,
#                                                          config.params.background_correction.output_multiplier
#                                                          )
#
#         self.event_analyser = lambda dataset, dataset_map, ref_map, events, grid: analyse_dataset_events(dataset,
#                                                                                                          dataset_map,
#                                                                                                          ref_map,
#                                                                                                          events,
#                                                                                                          grid,
#                                                                                                          max_bdc=config.params.background_correction.max_bdc,
#                                                                                                          min_bdc=config.params.background_correction.min_bdc,
#                                                                                                          increment=config.params.background_correction.increment,
#                                                                                                          output_multiplier=config.params.background_correction.output_multiplier,
#                                                                                                          )
#
#         self.map_maker = PanDDADefaultMapMaker()
#
#         self.event_table_maker = PanDDADefaultEventTableShell(order_by=config.results.events.order_by)
#
#         self.event_model = PanDDAEventModel(self.statistical_model,
#                                             self.clusterer,
#                                             self.event_finder,
#                                             statistics=[])
#
#         # ############################
#         # Criticism
#         # ############################
#         self.map_loader = MapLoaderDask(config.settings.verbose,
#                                         config.params.maps.resolution_factor,
#                                         config.params.maps.density_scaling,
#                                         )
#
#         # self.criticiser = Criticiser()
#         # self.criticiser_all = CriticiserAll(self.map_loader,
#         #                                     event_table_maker,
#         #                                     map_maker)
#
#
# class PanDDAConfig_dep:
#
#     def __init__(self, phil):
#
#         # ============================================================================>
#         # Transform Variables
#         # ============================================================================>
#         # input
#         # input = phil.input
#         # data_dirs = str(phil.input.data_dirs)
#         # pdb_style = phil.input.pdb_style
#         # mtz_style = phil.input.mtz_style
#         #
#         # # # input.regex
#         # regex = input.regex
#         # pdb_regex = regex.pdb_regex
#         # mtz_regex = regex.mtz_regex
#         # dir_regex = regex.dir_regex
#         #
#         # # # input.flags
#         # flags = input.flags
#         # only_datasets = flags.only_datasets
#         config = Config(phil)
#         #
#         # cpus = int(config["args"]["cpus"])
#         #
#         # resolution_factor = config["maps"]["resolution_factor"]
#         # density_scaling = config["maps"]["density_scaling"]
#
#         # ============================================================================>
#         # Dataset
#         # ============================================================================>
#         # # Dataloader
#         if config["dataset"]["dataloader"] == "default":
#             self.dataloader = DefaultDataloader(config.input.data_dirs,
#                                                 config.input.pdb_style,
#                                                 config.input.mtz_style,
#                                                 config.input.regex.pdb_regex,
#                                                 config.input.regex.mtz_regex,
#                                                 config.input.regex.dir_regex,
#                                                 config.input.flags.only_datasets,
#                                                 config.input.flags.ignore_datasets,
#                                                 config.output.dataset_prefix,
#                                                 config.output.out_dir,
#                                                 config.input.lig_style
#                                                 )
#
#         # # Reference
#         if config["dataset"]["get_reference"] == "default":
#             self.get_reference = DefaultReferenceGetter(out_dir=config.output.out_dir,
#                                                         reference_pdb_path=None,
#                                                         reference_mtz_path=None,
#                                                         reference_structure_factors=config.params.diffraction_data.structure_factors,
#                                                         structure_factors=config.params.diffraction_data.structure_factors
#                                                         )
#
#         # # Output
#         if config["output"]["version"] == "pandda":
#             self.pandda_output = PanddaOutputSetup(config.output.out_dir,
#                                                    config.input.lig_style
#                                                    )
#
#         # # Transforms
#         self.dataset_transforms = OrderedDict()
#         if config["transforms"]["PanddaDataChecker"] == "True":
#             self.dataset_transforms["data_check"] = PanddaDataChecker(config.params.diffraction_data.structure_factors,
#                                                                       config.params.diffraction_data.low_resolution_completeness,
#                                                                       config.params.diffraction_data.all_data_are_valid_values
#                                                                       )
#
#         if config["transforms"]["PanddaDiffractionScaler"] == "True":
#             self.dataset_transforms["scale_diffraction"] = PanddaDiffractionScaler(
#                 config.params.diffraction_data.apply_b_factor_scaling
#             )
#
#         if config["transforms"]["PanddaDatasetFilterer"] == "True":
#             self.dataset_transforms["filter_structure"] = PanddaDatasetFilterer(
#                 config.params.filtering.flags.similar_models_only,
#                 config.params.filtering.max_rfree,
#                 config.params.filtering.flags.same_space_group_only
#             )
#
#         if config["transforms"]["PanddaDatasetFiltererWilsonRMSD"] == "True":
#             self.dataset_transforms["filter_wilson"] = PanddaDatasetFiltererWilsonRMSD(
#                 config.params.excluding.max_wilson_plot_z_score,
#                 config.params.diffraction_data.apply_b_factor_scaling,
#                 config.params.diffraction_data.structure_factors
#             )
#
#         if config["transforms"]["PanddaDefaultStructureAligner"] == "True":
#             method = str(config["transforms"]["align_method"])
#             self.dataset_transforms["align"] = PanddaDefaultStructureAligner(
#                 method=config.params.alignment.method,
#                 cpus=config.settings.cpus
#             )
#
#             # # Sample loader
#         if config["dataset"]["sample_loader"] == "default":
#             grid_loader = PanDDAGridSetup(cpus=int(config.settings.cpus),
#                                           mask_pdb=config.params.masks.pdb,
#                                           align_mask_to_reference=bool(config.params.masks.align_mask_to_reference),
#                                           alignment_method=config.params.alignment.method,
#                                           outer_mask=float(config.params.masks.outer_mask),
#                                           inner_mask=float(config.params.masks.inner_mask),
#                                           inner_mask_symmetry=float(config.params.masks.inner_mask_symmetry),
#                                           grid_spacing=float(config.params.maps.grid_spacing),
#                                           padding=float(config.params.maps.padding),
#                                           verbose=bool(config.settings.verbose),
#                                           mask_selection_string=str(config.params.z_map_anlysis.masks.selection_string)
#                                           )
#
#             self.sample_loader = DefaultSampleLoader(config.params.maps.resolution_factor,
#                                                      config.params.maps.density_scaling,
#                                                      int(config.settings.cpus),
#                                                      bool(config.settings.verbose),
#                                                      grid_loader
#                                                      )
#
#         # # Refernce map loader
#         self.reference_map_getter = PanddaReferenceMapLoader(config.params.maps.resolution_factor,
#                                                              config.params.maps.density_scaling)
#
#         # # Partititioner
#         test = str(config["partitioner"]["test"]).split(",") if (str(config["partitioner"]["test"]) != "None") else None
#         train = str(config["partitioner"]["train"]).split(",") if (
#                     str(config["partitioner"]["train"]) != None) else None
#         not_test = str(config["partitioner"]["not_test"]).split(",") if bool(
#             str(config["partitioner"]["not_test"])) else None
#         not_train = str(config["partitioner"]["not_train"]).split(",") if bool(
#             str(config["partitioner"]["not_train"])) else None
#
#         if config["dataset"]["partitioner"] == "default":
#             self.partitioner = DefaultPanDDAPartitionsGetter(test=test,
#                                                              train=train,
#                                                              not_test=not_test,
#                                                              not_train=not_train
#                                                              )
#
#         # ============================================================================>
#         # Event Model
#         # ============================================================================>
#         bdc_calculator = config["event_model"]["bdc_calculator"]
#         max_bdc = float(config["bdc_calculator"]["max_bdc"])
#         min_bdc = float(config["bdc_calculator"]["min_bdc"])
#         increment = float(config["bdc_calculator"]["increment"])
#         output_multiplier = float(config["bdc_calculator"]["output_multiplier"])
#
#         map_maker = config["event_model"]["map_maker"]
#         order_by = config["event_table"]["order_by"]
#
#         event_table_maker = config["event_model"]["event_table_maker"]
#
#         if config["event_model"]["model"] == "PanDDAEventModel":
#
#             self.statistical_model = None
#             if config["event_model"]["statistical_model"] == "PanDDANormalModel":
#                 self.statistical_model = PanDDANormalModel(method="adjusted+uncertainty", cpus=1)
#
#             self.clusterer = None
#             if config["event_model"]["clusterer"] == "PanDDADefaultCluster":
#                 self.clusterer = PanDDADefaultCluster(grid_clustering_cutoff=float(
#                     config.params.z_map_anlysis.agglomerative_hierarchical.clustering_cutoff),
#                                                       negative_values=bool(config.params.z_map_anlysis.negative_values),
#                                                       cluster_method=config.params.z_map_anlysis.clustering_method,
#                                                       contour_level=float(config.params.z_map_anlysis.contour_level),
#                                                       outer_mask=float(config.params.masks.outer_mask),
#                                                       inner_mask_symmetry=float(config.params.masks.inner_mask_symmetry)
#                                                       )
#
#             self.event_finder = None
#             if config["event_model"]["event_finder"] == "PanDDADefaultEventFinder":
#                 self.event_finder = PanDDADefaultEventFinder(
#                     min_blob_z_peak=float(config.params.z_map_anlysis.min_blob_z_peak),
#                     grid_minimum_volume=float(config.params.z_map_anlysis.min_blob_volume),
#                     grid_spacing=float(config.params.maps.grid_spacing)
#                     )
#
#             if bdc_calculator == "PanDDADefaultBDCCalculator":
#                 self.bdc_calculator = PanDDADefaultBDCCalculator(config.params.background_correction.max_bdc,
#                                                                  config.params.background_correction.min_bdc,
#                                                                  config.params.background_correction.increment,
#                                                                  config.params.background_correction.output_multiplier
#                                                                  )
#
#             if map_maker == "PanDDADefaultMapMaker":
#                 self.map_maker = PanDDADefaultMapMaker()
#
#             if event_table_maker == "PanDDADefaultEventTableShell":
#                 self.event_table_maker = PanDDADefaultEventTableShell(order_by=order_by)
#
#             self.event_model = PanDDAEventModel(self.statistical_model,
#                                                 self.clusterer,
#                                                 self.event_finder,
#                                                 statistics=[])
#
#         # ############################
#         # Criticism
#         # ############################
#         verbose = bool(config["args"]["verbose"])
#         resolution_factor = float(config["maps"]["resolution_factor"])
#         density_scaling = str(config["maps"]["density_scaling"])
#
#         self.map_loader = MapLoaderDask(verbose,
#                                         resolution_factor,
#                                         density_scaling)
#
#         self.criticiser = Criticiser()
#         self.criticiser_all = CriticiserAll(self.map_loader,
#                                             event_table_maker,
#                                             map_maker)

# class PanDDAConfig:
#
#     def __init__(self, config):
#
#         # ============================================================================>
#         # Transform Variables
#         # ============================================================================>
#         cpus = int(config["args"]["cpus"])
#
#         resolution_factor = config["maps"]["resolution_factor"]
#         density_scaling = config["maps"]["density_scaling"]
#
#         # ============================================================================>
#         # Dataset
#         # ============================================================================>
#         # # Dataloader
#         if config["dataset"]["dataloader"] == "default":
#             self.dataloader = DefaultDataloader(config["dataloader"]["data_dirs"],
#                                                 config["dataloader"]["pdb_style"],
#                                                 config["dataloader"]["mtz_style"],
#                                                 config["dataloader"]["pdb_regex"],
#                                                 config["dataloader"]["mtz_regex"],
#                                                 config["dataloader"]["dir_regex"],
#                                                 config["dataloader"]["only_datasets"],
#                                                 config["dataloader"]["ignore_datasets"],
#                                                 config["dataloader"]["dataset_prefix"],
#                                                 config["dataloader"]["out_dir"],
#                                                 config["dataloader"]["lig_style"])
#
#         # # Reference
#         if config["dataset"]["get_reference"] == "default":
#             self.get_reference = DefaultReferenceGetter(out_dir=config["dataloader"]["out_dir"],
#                                                         reference_pdb_path=None,
#                                                         reference_mtz_path=None,
#                                                         reference_structure_factors=config["diffraction_data"][
#                                                             "structure_factors"],
#                                                         structure_factors=config["diffraction_data"][
#                                                             "structure_factors"])
#
#         # # Output
#         if config["output"]["version"] == "pandda":
#             self.pandda_output = PanddaOutputSetup(config["dataloader"]["out_dir"],
#                                                    config["dataloader"]["lig_style"])
#
#         # # Transforms
#         self.dataset_transforms = OrderedDict()
#         if config["transforms"]["PanddaDataChecker"] == "True":
#             self.dataset_transforms["data_check"] = PanddaDataChecker(config["diffraction_data"]["structure_factors"],
#                                                                       config["diffraction_data"][
#                                                                           "low_resolution_completeness"],
#                                                                       config["diffraction_data"][
#                                                                           "all_data_are_valid_values"]
#                                                                       )
#
#         if config["transforms"]["PanddaDiffractionScaler"] == "True":
#             self.dataset_transforms["scale_diffraction"] = PanddaDiffractionScaler(
#                 config["diffraction_data"]["apply_b_factor_scaling"]
#             )
#
#         if config["transforms"]["PanddaDatasetFilterer"] == "True":
#             self.dataset_transforms["filter_structure"] = PanddaDatasetFilterer(
#                 config["dataset_filtering"]["similar_models_only"],
#                 config["dataset_filtering"]["max_rfree"],
#                 config["dataset_filtering"]["same_space_group_only"]
#             )
#
#         if config["transforms"]["PanddaDatasetFiltererWilsonRMSD"] == "True":
#             self.dataset_transforms["filter_wilson"] = PanddaDatasetFiltererWilsonRMSD(
#                 config["diffraction_data"]["max_wilson_plot_z_score"],
#                 config["diffraction_data"]["apply_b_factor_scaling"],
#                 config["diffraction_data"]["structure_factors"]
#             )
#
#         if config["transforms"]["PanddaDefaultStructureAligner"] == "True":
#             method = str(config["transforms"]["align_method"])
#             self.dataset_transforms["align"] = PanddaDefaultStructureAligner(
#                 method=method,
#                 cpus=cpus
#             )
#
#             # # Sample loader
#         if config["dataset"]["sample_loader"] == "default":
#             grid_loader = PanDDAGridSetup(cpus=int(config["args"]["cpus"]),
#                                           mask_pdb=config["grid"]["mask_pdb"],
#                                           align_mask_to_reference=bool(config["grid"]["align_mask_to_reference"]),
#                                           alignment_method=config["grid"]["alignment_method"],
#                                           outer_mask=float(config["grid"]["outer_mask"]),
#                                           inner_mask=float(config["grid"]["inner_mask"]),
#                                           inner_mask_symmetry=float(config["grid"]["inner_mask_symmetry"]),
#                                           grid_spacing=float(config["grid"]["grid_spacing"]),
#                                           padding=float(config["grid"]["padding"]),
#                                           verbose=bool(config["args"]["verbose"]),
#                                           mask_selection_string=str(config["grid"]["mask_selection_string"])
#                                           )
#
#             self.sample_loader = DefaultSampleLoader(config["maps"]["resolution_factor"],
#                                                      config["maps"]["density_scaling"],
#                                                      int(config["args"]["cpus"]),
#                                                      bool(config["args"]["verbose"]),
#                                                      grid_loader
#                                                      )
#
#         # # Refernce map loader
#         self.reference_map_getter = PanddaReferenceMapLoader(resolution_factor,
#                                                              density_scaling)
#
#         # # Partititioner
#         test = str(config["partitioner"]["test"]).split(",") if bool(str(config["partitioner"]["test"])) else None
#         train = str(config["partitioner"]["train"]).split(",") if bool(str(config["partitioner"]["train"])) else None
#         not_test = str(config["partitioner"]["not_test"]).split(",") if bool(str(config["partitioner"]["not_test"])) else None
#         not_train = str(config["partitioner"]["not_train"]).split(",") if bool(str(config["partitioner"]["not_train"])) else None
#
#         if config["dataset"]["partitioner"] == "default":
#             self.partitioner = DefaultPanDDAPartitionsGetter(test=test,
#                                                              train=train,
#                                                              not_test=not_test,
#                                                              not_train=not_train
#                                                              )
#
#
#         # ============================================================================>
#         # Event Model
#         # ============================================================================>
#         bdc_calculator = config["event_model"]["bdc_calculator"]
#         max_bdc = float(config["bdc_calculator"]["max_bdc"])
#         min_bdc = float(config["bdc_calculator"]["min_bdc"])
#         increment = float(config["bdc_calculator"]["increment"])
#         output_multiplier = float(config["bdc_calculator"]["output_multiplier"])
#
#         map_maker = config["event_model"]["map_maker"]
#         order_by = config["event_table"]["order_by"]
#
#         event_table_maker = config["event_model"]["event_table_maker"]
#
#         if config["event_model"]["model"] == "PanDDAEventModel":
#
#             self.statistical_model = None
#             if config["event_model"]["statistical_model"] == "PanDDANormalModel":
#                 self.statistical_model = PanDDANormalModel(method="adjusted+uncertainty", cpus=1)
#
#             self.clusterer = None
#             if config["event_model"]["clusterer"] == "PanDDADefaultCluster":
#                 self.clusterer = PanDDADefaultCluster(grid_clustering_cutoff=float(config["cluster"]["grid_clustering_cutoff"]),
#                                                       negative_values=bool(config["cluster"]["negative_values"]),
#                                                       cluster_method=config["cluster"]["method"],
#                                                       contour_level=float(config["cluster"]["contour_level"]),
#                                                       outer_mask=float(config["grid"]["outer_mask"]),
#                                                       inner_mask_symmetry=float(config["grid"]["inner_mask_symmetry"])
#                                                       )
#
#             self.event_finder = None
#             if config["event_model"]["event_finder"] == "PanDDADefaultEventFinder":
#                 self.event_finder = PanDDADefaultEventFinder(min_blob_z_peak=float(config["event"]["min_blob_z_peak"]),
#                                                              grid_minimum_volume=float(config["event"]["min_blob_volume"]),
#                                                              grid_spacing=float(config["grid"]["grid_spacing"])
#                                                              )
#
#             if bdc_calculator == "PanDDADefaultBDCCalculator":
#                 self.bdc_calculator = PanDDADefaultBDCCalculator(max_bdc, min_bdc, increment, output_multiplier)
#
#             if map_maker == "PanDDADefaultMapMaker":
#                 self.map_maker = PanDDADefaultMapMaker()
#
#             if event_table_maker == "PanDDADefaultEventTableShell":
#                 self.event_table_maker = PanDDADefaultEventTableShell(order_by=order_by)
#
#             self.event_model = PanDDAEventModel(self.statistical_model,
#                                                 self.clusterer,
#                                                 self.event_finder,
#                                                 statistics=[])
#
#
#
#         # ############################
#         # Criticism
#         # ############################
#         verbose = bool(config["args"]["verbose"])
#         resolution_factor = float(config["maps"]["resolution_factor"])
#         density_scaling = str(config["maps"]["density_scaling"])
#
#         self.map_loader = MapLoaderDask(verbose,
#                                         resolution_factor,
#                                         density_scaling)
#
#         self.criticiser = Criticiser()
#         self.criticiser_all = CriticiserAll(self.map_loader,
#                                             event_table_maker,
#                                             map_maker)
