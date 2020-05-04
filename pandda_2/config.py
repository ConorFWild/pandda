import os, sys, copy

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
        print('no version information available')
    else:
        print('{} version: {}'.format(module_info.name, module_info.version))
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
        self.ground_state_datasets = config_obj.ground_state_datasets
        self.exclude_from_characterisation = config_obj.exclude_from_characterisation


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
        self.statistical_maps = StatisticalMaps(config_obj.statistical_maps)


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


class StatisticalMaps:
    def __init__(self, config_obj):
        self.min_build_datasets = config_obj.min_build_datasets
        self.max_build_datasets = config_obj.max_build_datasets


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
        self.process_shells = config_obj.process_shells
        self.process_dict_n_cpus = config_obj.process_dict_n_cpus
        self.process_dict = config_obj.process_dict
        self.h_vmem = config_obj.h_vmem
        self.m_mem_free = config_obj.m_mem_free


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
