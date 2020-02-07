import os, sys, copy

import libtbx.phil
from libtbx.utils import Sorry

from giant import module_info

class run_default(object):

    _module_info = module_info

    def __init__(self, run, master_phil, args, blank_arg_prepend=None, program='', description=''):
        """Run a program via a standard setup of functions and objects"""
        print(self._module_info.header_text.format(program=program, description=description))
        working_phil = extract_params_default(master_phil=master_phil, args=args, blank_arg_prepend=blank_arg_prepend, module_info=self._module_info)
        out = run(params=working_phil.extract())


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
    if '-a' in args: expert_level = 3

    if ('?' in args) or ('--show-defaults' in args) or (not args) or (not attributes_level==expert_level==0):
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

