import os, sys, copy

import iotbx.pdb
import libtbx.phil

from bamboo.common.logs import Log

from giant.io.pdb import strip_pdb_to_input
from giant.structure.utils import resolve_residue_id_clashes, transfer_residue_groups_from_other
from giant.structure.altlocs import expand_alternate_conformations, find_next_conformer_idx, increment_altlocs, prune_redundant_alternate_conformations
from giant.structure.occupancy import set_conformer_occupancy

from giant.jiffies import make_restraints

############################################################################

PROGRAM = 'giant.merge_conformations'

DESCRIPTION = """
    A tool to merge two models of the same crystal into one structure.
        The alternate conformers of the two models are generated and merged. Any common sections are then
        converted back into a single conformer model. This tool is therefore particularly useful when one
        model is an edited version of the other as many of the residues will be identical and the output
        structure will not contain too many alternate conformations.

    1) Simple usage:
        > giant.merge_conformations major=structure_1.pdb minor=structure_2.pdb
    2) Lazy usage -- first file with be used as the major conformation
        > giant.merge_conformations structure_1.pdb structure_2.pdb
"""

############################################################################

blank_arg_prepend = {'.pdb': 'input.pdb='}

master_phil = libtbx.phil.parse("""
input {
    major = None
        .help = 'The major conformation of the protein (the conformation you want to have the lowest altlocs)'
        .type = str
    minor = None
        .help = 'The minor conformation of the protein (the conformation that will have the higher altlocs)'
        .type = str
    pdb = None
        .help = 'Dummy holder to allow files to be specified without major= or minor=. The first file will be used as major.'
        .type = str
        .multiple = True
}
options {
    minor_occupancy = 0.5
        .help = 'occupancy of the minor state after merging (multiplies existing occupancies)'
        .type = float
    major_occupancy = 0.5
        .help = 'occupancy of the major state after merging (multiplies existing occupancies)'
        .type = float
    prune_duplicates_rmsd = 0.05
        .help = 'rmsd at which to remove duplicate residues'
        .type = float
    reset_all_occupancies = False
        .help = 'set all conformer occupancies to the same value (1.0/<number of conformers>). overrides minor_occupancy=X and major_occupancy=X.'
        .type = bool
}
output {
    pdb = 'multi-state-model.pdb'
        .help = 'output pdb file'
        .type = str
    log = 'multi-state-merge-conformations.log'
        .help = 'output log file'
        .type = str
    make_restraints = True
        .help = 'generate refinement restraints for the output model'
        .type = bool
}
restraints {
    output {
    include scope giant.jiffies.make_restraints.output_phil
    }
    include scope giant.jiffies.make_restraints.options_phil
}
settings {
    overwrite = True
        .type = bool
    verbose = False
        .type = bool
}
""", process_includes=True)

############################################################################

def merge_complementary_hierarchies(hierarchy_1, hierarchy_2, prune_duplicates_rmsd=0.1, in_place=False, verbose=False, log=None):
    """Merge hierarchies that are alternate models of the same crystal by expanding alternate model conformations, merging, and then trimming alternate conformations where possible"""

    if log is None: log = Log(verbose=True)

    # Alter the original files?
    if not in_place:
        # Copy the hierarchies
        hierarchy_1 = hierarchy_1.deep_copy()
        hierarchy_2 = hierarchy_2.deep_copy()

    # Sort the atoms
    hierarchy_1.sort_atoms_in_place()
    hierarchy_2.sort_atoms_in_place()

    log.heading('Preparing to merge structures')

    log.subheading('Explicitly expanding models to all conformations of the crystal')
    log('Expanding alternate conformations in structure 1')
    expand_alternate_conformations(
        hierarchy   = hierarchy_1,
        in_place    = True,
        verbose     = verbose)
    log('Expanding alternate conformations in structure 2')
    expand_alternate_conformations(
        hierarchy   = hierarchy_2,
        in_place    = True,
        verbose     = verbose)
    log.subheading('Applying conformer shift to the second structure before merging')
    log('Identifying the altloc shift required from the number of alternate conformers in structure 1')
    conf_offset = find_next_conformer_idx(
        hierarchy           = hierarchy_1,
        all_ids             = iotbx.pdb.systematic_chain_ids())
    log('Incrementing all altlocs in structure 2 by {}'.format(conf_offset))
    increment_altlocs(
        hierarchy           = hierarchy_2,
        offset              = conf_offset,
        in_place            = True,
        verbose             = verbose)
    log.subheading('Renaming residues that do not align between structures')
    resolve_residue_id_clashes(
        fixed_hierarchy     = hierarchy_1,
        moving_hierarchy    = hierarchy_2,
        in_place            = True,
        verbose             = verbose)

    log.heading('Merging structures')

    log('Transferring residues from Structure 2 to Structure 1')
    transfer_residue_groups_from_other(
        acceptor_hierarchy  = hierarchy_1,
        donor_hierarchy     = hierarchy_2,
        in_place            = True,
        verbose             = verbose)

    log.heading('Post-processing structure')

    log('Pruning unneccessary multi-conformer residues in the merged structure')
    prune_redundant_alternate_conformations(
        hierarchy           = hierarchy_1,
        required_altlocs    = hierarchy_1.altloc_indices(),
        rmsd_cutoff         = prune_duplicates_rmsd,
        in_place            = True,
        verbose             = verbose)

    return hierarchy_1

############################################################################

def run(params):

    # Create log file
    log = Log(log_file=params.output.log, verbose=True)

    # Report
    log.heading('Validating input parameters and input files')

    # Check one or other have been provided
    if (params.input.major or params.input.minor) and not (params.input.pdb==[None] or params.input.pdb==[]):
        raise Exception('Have provided input.major & input.minor, as well as files to input.pdb. Specify either input.major & input.minor, or two input.pdb.')
    # Assign files to major and minor if necessary
    if not (params.input.major and params.input.minor):
        if len(params.input.pdb) != 2:
            raise Exception('Must provide zero or two pdb files to input.pdb')
        params.input.major = params.input.pdb[0]
        params.input.minor = params.input.pdb[1]
    # Check files exist
    if not os.path.exists(params.input.major):
        raise Exception('input.major does not exist: {}'.format(params.input.major))
    if not os.path.exists(params.input.minor):
        raise Exception('input.minor does not exist: {}'.format(params.input.minor))
    # Just check again...
    assert params.input.major
    assert params.input.minor
    assert params.output.pdb
    # Check existence of output pdb and delete as necessary
    if os.path.exists(params.output.pdb):
        if params.settings.overwrite:
            os.remove(params.output.pdb)
        else:
            raise Exception('Output file already exists: {}. Run with overwrite=True to remove this file'.format(params.output.pdb))

    # Check that the input occupancies are valid
    if (params.options.minor_occupancy>1.0) or (params.options.major_occupancy>1.0):
        raise Exception('minor_occupancy and major_occupancy cannot be greater than 1.0 (currently {} and {})'.format(params.options.minor_occupancy,params.options.major_occupancy))

    # Report validated parameters
    log.subheading('Processed merging parameters')
    for obj in master_phil.format(params).objects:
        if obj.name == 'restraints': continue
        log(obj.as_str().strip())

    # Read in the ligand file and set each residue to the requested conformer
    log.subheading('Reading input files')
    maj_obj = strip_pdb_to_input(params.input.major, remove_ter=True)
    min_obj = strip_pdb_to_input(params.input.minor, remove_ter=True)

    # Check that ... something
    try:
        maj_obj.hierarchy.only_model()
        min_obj.hierarchy.only_model()
    except:
        raise Sorry('Input structures may only have one model')

    # Multiply the input hierarchies by occupancy multipliers
    log.subheading('Updating input occupancies prior to merging')
    log('Multiplying occupancies of input.major by {}'.format(params.options.major_occupancy))
    maj_obj.hierarchy.atoms().set_occ(maj_obj.hierarchy.atoms().extract_occ()*params.options.major_occupancy)
    log('Multiplying occupancies of input.minor by {}'.format(params.options.minor_occupancy))
    min_obj.hierarchy.atoms().set_occ(min_obj.hierarchy.atoms().extract_occ()*params.options.minor_occupancy)

    # Merge the hierarchies
    final_struct = merge_complementary_hierarchies(
        hierarchy_1 = maj_obj.hierarchy,
        hierarchy_2 = min_obj.hierarchy,
        prune_duplicates_rmsd = params.options.prune_duplicates_rmsd,
        in_place    = True,
        verbose     = params.settings.verbose)

    # Set output occupancies
    log.subheading('Post-processing occupancies')
    # Set all main-conf occupancies to 1.0
    log('Setting all main-conf occupancies to 1.0')
    set_conformer_occupancy(hierarchy   = final_struct,
                            altlocs     = [''],
                            occupancy   = 1.0,
                            in_place    = True,
                            verbose     = params.settings.verbose)
    # Reset occupancies if required
    if params.options.reset_all_occupancies:
        # Calculate number of altlocs and associated occupancy
        altlocs = [a for a in final_struct.altloc_indices() if a]
        if altlocs:
            new_occ = 1.0/len(altlocs)
            # Set the occupancies
            log('Setting all conformer ({}) occupancies to {}'.format(','.join(altlocs), new_occ))
            set_conformer_occupancy(hierarchy   = final_struct,
                                    altlocs     = altlocs,
                                    occupancy   = new_occ,
                                    in_place    = True,
                                    verbose     = params.settings.verbose)

    # Update the atoms numbering
    final_struct.sort_atoms_in_place()
    final_struct.atoms_reset_serial()
    # Write output file
    log('Writing output structure to {}'.format(params.output.pdb))
    final_struct.write_pdb_file(file_name=params.output.pdb, crystal_symmetry=maj_obj.crystal_symmetry())

    # Run the restraint generation for the merged structure if requested
    if params.output.make_restraints:

        # Transfer the other phil objects from the master phil
        r_params = make_restraints.master_phil.extract()
        for name, obj in r_params.__dict__.items():
            if name.startswith('_'): continue
            if name not in params.restraints.__dict__:
                params.restraints.__inject__(name, obj)

        # Apply the output of merging to input of restraints
        params.restraints.input.pdb = params.output.pdb
        # Rename output files to be in same folder as output structure
        if params.restraints.output.phenix:
            params.restraints.output.phenix = os.path.join(os.path.dirname(params.output.pdb), os.path.basename(params.restraints.output.phenix))
        if params.restraints.output.refmac:
            params.restraints.output.refmac = os.path.join(os.path.dirname(params.output.pdb), os.path.basename(params.restraints.output.refmac))
        # Set log file name to this program if one given
        if params.output.log:
            params.restraints.output.log = params.output.log
        elif params.restraints.output.log:
            params.restraints.output.log = os.path.join(os.path.dirname(params.output.pdb), os.path.basename(params.restraints.output.log))
        # Which alternate conformations to generate restraints for
        params.restraints.local_restraints.altlocs = ','.join([a for a in min_obj.hierarchy.altloc_indices() if a])
        # Update settigns
        params.restraints.settings.verbose = params.settings.verbose
        params.restraints.settings.overwrite = params.settings.overwrite

        # Report
        log.heading('Parameters for generating restraints')
        log(master_phil.format(params).get('restraints').as_str().strip())
        log.heading('Generating restraints')
        # Run make_restraints
        make_restraints.run(params.restraints)

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
