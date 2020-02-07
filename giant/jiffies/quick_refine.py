import os, sys, glob, shutil

import libtbx.phil

from bamboo.common.command import CommandManager
from bamboo.common.logs import Log
from bamboo.common.path import rel_symlink

from giant.jiffies import split_conformations

############################################################################

PROGRAM = 'giant.quick_refine'
DESCRIPTION = """
    A tool to simplify the launching of standard refinement jobs in REFMAC or PHENIX.

    1) Simple usage:
        > giant.quick_refine input.pdb input.mtz ligand.cif

    2) Defining custom names for the output folder and the output files
        > giant.quick_refine ... dir_prefix='XXX' out_prefix='XXX'

    3) Specify additonal parameter file (see giant.make_restraints)
        > giant.quick_refine ... params=restraints.params
"""

blank_arg_prepend = {   '.mtz':'input.mtz=',
                        '.pdb':'input.pdb=',
                        '.cif':'input.cif=',
                        '.params':'input.params=' }

############################################################################

master_phil = libtbx.phil.parse("""
input {
    pdb = None
        .type = str
    mtz = 'refine.mtz'
        .type = str
    cif = None
        .type = str
        .multiple = True
    params = None
        .help = "File that contains additional parameters for the refinement program."
        .type = str
    args = None
        .help = "Pass any additional arguments to the program? (Command Line arguments). Pass as quoted strings."
        .type = str
        .multiple = True
}
output {
    dir_prefix = 'refine_'
        .type = str
    out_prefix = 'output'
        .type = str
    link_prefix = 'refine'
        .type = str
}
options {
    program = *phenix refmac
        .type = choice
    split_conformations = False
        .help = "Split the output structure into different conformations for modelling."
        .type = bool
}
split_conformations {
    include scope giant.jiffies.split_conformations.master_phil
}
settings {
    verbose = True
        .type = bool
}
""", process_includes=True)

############################################################################

def run(params):

    # Identify any existing output directories
    current_dirs = sorted(glob.glob(params.output.dir_prefix+'*'))
    if not current_dirs:
        next_int = 1
    else:
        current_nums = [s.replace(params.output.dir_prefix, '') for s in current_dirs]
        next_int = sorted(map(int, current_nums))[-1]+1

    # Create output directory name from int
    out_dir = params.output.dir_prefix + '{:04}'.format(next_int)
    # Create output directory
    os.mkdir(out_dir)

    # Create log object
    log = Log(log_file=os.path.join(out_dir, params.output.out_prefix+'.quick-refine.log'), verbose=True)

    # Report
    if current_dirs:
        log('Found existing refinement directories: \n\t{}'.format('\n\t'.join(current_dirs)))
        log('')
    log('Creating new output directory: {}'.format(out_dir))

    # Validate input parameters
    log.subheading('Validating input parameters')
    assert params.input.pdb is not None, 'No PDB given for refinement'
    assert params.input.mtz is not None, 'No MTZ given for refinement'

    if os.path.islink(params.input.mtz):
        log('Converting mtz path to real path:')
        log('{} -> {}'.format(params.input.mtz, os.path.realpath(params.input.mtz)))
        params.input.mtz = os.path.realpath(params.input.mtz)

    # Link input
    log('Copying/linking files to refinement folder')
    shutil.copy(params.input.pdb, os.path.abspath(os.path.join(out_dir, 'input.pdb')))
    rel_symlink(params.input.mtz, os.path.abspath(os.path.join(out_dir, 'input.mtz')))
    # Copy parameter file to output folder
    if params.input.params:
        shutil.copy(params.input.params, os.path.abspath(os.path.join(out_dir, 'input.params')))

    # Create output prefixes
    output_prefix = os.path.join(out_dir, params.output.out_prefix)
    log('Real output file path prefixes: {}'.format(output_prefix))
    log('Link output file path prefixes: {}'.format(params.output.link_prefix))

    # Create command objects
    log.subheading('Preparing command line input for refinement program')

    # PHENIX
    if params.options.program == 'phenix':
        cm = CommandManager('phenix.refine')
        # Command line args
        cm.add_command_line_arguments([ params.input.pdb, params.input.mtz ])
        cm.add_command_line_arguments([ 'output.prefix={}'.format(output_prefix) ])
        if params.input.cif:
            cm.add_command_line_arguments( params.input.cif )
        if params.input.params and os.path.exists(params.input.params):
            cm.add_command_line_arguments([ params.input.params ])

    # REFMAC
    elif params.options.program == 'refmac':
        cm = CommandManager('refmac5')
        # Command line args
        cm.add_command_line_arguments( ['xyzin', params.input.pdb, 'hklin', params.input.mtz] )
        cm.add_command_line_arguments( ['xyzout', output_prefix+'.pdb', 'hklout', output_prefix+'.mtz'] )
        if params.input.cif:
            for cif in params.input.cif:
                cm.add_command_line_arguments( ['libin', cif] )
        # Standard input
        if params.input.params:
            cm.add_standard_input( open(params.input.params).read().split('\n') )

        cm.add_standard_input( ['END'] )

    # Pass additional command line arguments?
    if params.input.args:
        cm.add_command_line_arguments( params.input.args )

    # Report
    log(str(cm))

    log.subheading('running refinement - {}'.format(cm.program[0]))
    out = cm.run()

    log.subheading('refinement output')
    log(cm.output)

    if out != 0:
        log.subheading('errors')
        log(cm.error)

    # Find output files
    try:
        real_pdb = glob.glob(output_prefix+'*.pdb')[0]
        real_mtz = glob.glob(output_prefix+'*.mtz')[0]
    except:
        print('Refinement has failed - outptu files do not exist')
        print(real_pdb)
        print(real_mtz)
        raise

    # List of links to make at the end of the run
    link_file_pairs = [     (real_pdb, params.output.link_prefix+'.pdb'),
                            (real_mtz, params.output.link_prefix+'.mtz')    ]

    # Split conformations
    if params.options.split_conformations:
        params.split_conformations.settings.verbose = params.settings.verbose
        log.subheading('splitting refined structure conformations')
        # Running split conformations
        out_files = split_conformations.split_conformations(filename=real_pdb, params=params.split_conformations, log=log)
        # Link output files to top
        for real_file in out_files:
            link_file = params.output.link_prefix+os.path.basename(real_file.replace(os.path.splitext(real_pdb)[0], ''))
            link_file_pairs.append([real_file, link_file])

    # Link output files
    log.subheading('linking output files')
    for real_file, link_file in link_file_pairs:
        log('Linking {} -> {}'.format(link_file,real_file))
        if not os.path.exists(real_file):
            log('file does not exist: {}'.format(real_file))
            continue
        if os.path.exists(link_file) and os.path.islink(link_file):
            log('removing existing link: {}'.format(link_file))
            os.unlink(link_file)
        if not os.path.exists(link_file):
            rel_symlink(real_file, link_file)

    log.heading('finished - refinement')

#######################################

if __name__=='__main__':
    from giant.jiffies import run_default
    run_default(
        run                 = run,
        master_phil         = master_phil,
        args                = sys.argv[1:],
        blank_arg_prepend   = blank_arg_prepend,
        program             = PROGRAM,
        description         = DESCRIPTION)
