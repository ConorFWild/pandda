import os, sys, copy

import iotbx.pdb
import libtbx.phil
from libtbx.utils import Sorry

from bamboo.common.logs import Log

from giant.io.pdb import strip_pdb_to_input
from giant.structure.altlocs import expand_alternate_conformations, prune_redundant_alternate_conformations

from giant.jiffies import make_restraints

############################################################################

PROGRAM = 'giant.standardise_model'

DESCRIPTION = """
    A tool to standardise a multi-conformer model so that no backbone discontinuities are created, etc.

    1) Simple usage:
        > giant.standardise_model my.pdb
"""

############################################################################

blank_arg_prepend = {'.pdb': 'input.pdb='}

master_phil = libtbx.phil.parse("""
input {
    pdb = None
        .help = 'Input PDB files - multiple can be provided'
        .type = str
        .multiple = True
}
options {
    pruning_rmsd = 0.1
        .help = 'rmsd at which to prune all alternate conformations'
        .type = float
}
output {
    suffix = '.standardised'
        .help = 'output pdb file suffix'
        .type = str
    log = 'standardise-model.log'
        .help = 'output log file'
        .type = str
}
settings {
    overwrite = True
        .type = bool
    verbose = False
        .type = bool
}
""")

############################################################################

def standardise_multiconformer_model(hierarchy, pruning_rmsd=0.1, in_place=False, verbose=False, log=None):
    """Standardise hierarchies by expanding alternate model conformations, and then trimming alternate conformations where possible"""

    if log is None: log = Log(verbose=True)

    # Alter the original files?
    if not in_place:
        # Copy the hierarchies
        hierarchy = hierarchy.deep_copy()

    # Sort the atoms
    hierarchy.sort_atoms_in_place()

    log.heading('Preparing to standardise structure')

    log.subheading('Explicitly expanding model to all conformations of the crystal')
    expand_alternate_conformations(
        hierarchy   = hierarchy,
        in_place    = True,
        verbose     = verbose)

    log.subheading('Pruning unneccessary multi-conformer residues in the expanded structure')
    prune_redundant_alternate_conformations(
        hierarchy           = hierarchy,
        required_altlocs    = hierarchy.altloc_indices(),
        rmsd_cutoff         = pruning_rmsd,
        in_place            = True,
        verbose             = verbose)

    return hierarchy

############################################################################

def run(params):

    # Create log file
    log = Log(log_file=params.output.log, verbose=True)

    # Report
    log.heading('Validating input parameters and input files')

    # Check one or other have been provided
    assert params.input.pdb, 'No pdb files have been provided'
    for pdb in params.input.pdb:
        if not os.path.exists(pdb): raise Sorry('pdb does not exist: {}'.format(pdb))

    for pdb in params.input.pdb:

        log.subheading('Reading pdb: {}'.format(pdb))
        obj = strip_pdb_to_input(pdb, remove_ter=True)
        try:
            obj.hierarchy.only_model()
        except:
            raise Sorry('Input structures may only have one model')

        # Merge the hierarchies
        final =  standardise_multiconformer_model(
            hierarchy    = obj.hierarchy,
            pruning_rmsd = params.options.pruning_rmsd,
            in_place     = True,
            verbose      = params.settings.verbose)

        # Update the atoms numbering
        final.sort_atoms_in_place()

        # Write output file
        filename = os.path.splitext(pdb)[0]+params.output.suffix+'.pdb'
        log('Writing output structure to {}'.format(filename))
        final.write_pdb_file(file_name=filename, crystal_symmetry=obj.crystal_symmetry())

    log.heading('FINISHED')
    log.heading('Final Parameters')
    log(master_phil.format(params).as_str().strip())

    return

############################################################################

if __name__=='__main__':
    from giant.jiffies import run_default
    run_default(
        run                 = run,
        master_phil         = master_phil,
        args                = sys.argv[1:],
        blank_arg_prepend   = blank_arg_prepend,
        program             = PROGRAM,
        description         = DESCRIPTION)
