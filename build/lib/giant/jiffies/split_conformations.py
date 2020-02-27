import os, sys, copy

import numpy

import iotbx.pdb
import libtbx.phil

from bamboo.common.logs import Log

from giant.io.pdb import strip_pdb_to_input, get_pdb_header
from giant.structure import calculate_residue_group_occupancy
from giant.structure.formatting import Labeller
from giant.structure.iterators import residue_groups_with_complete_set_of_conformers
from giant.structure.altlocs import prune_redundant_alternate_conformations

############################################################################

PROGRAM = 'giant.strip_conformations'

DESCRIPTION = """
    A tool to remove unwanted conformations of a multi-conformer model.
        Conformations to keep can be kept can be declared explicity (using conf=...) or implicitly,
        through keeping the conformations associated with certain residue names (res=...).

    1) Simple usage:
        > giant.strip_conformations input.pdb
"""

############################################################################

blank_arg_prepend = {'.pdb':'pdb='}

master_phil = libtbx.phil.parse("""
input  {
    pdb = None
        .help = 'A model containing multiple states/conformations (for example unbound and bound states)'
        .type = str
        .multiple = True
}
output {
    suffix_prefix = 'split'
        .help = "prefix for the suffix to be appended to the output structure -- appended to the basename of the input structure before any conformer-specific labels"
        .type = str
    log = None
        .help = ""
        .type = path
}
options {
    mode = by_conformer by_conformer_group *by_residue_name
        .help = 'How to split the model:\n\tby_conformer : output a model for each conformer\n\tby_conformer_group : output a model for a selection of conformers\n\tby_residue_name : output a model containing all states containing a specifice residue name'
        .type = choice
    by_conformer_group {
        conformers = C,D,E,F,G
            .help = "Output these conformers in one model. Multiple can be provided."
            .multiple = True
            .type = str
    }
    by_residue_name {
        resname = DRG,FRG,LIG,UNK,UNL
            .help = "Group conformers containing any of these residue names"
            .type = str
        selected_name = 'bound-state'
            .help = "Output filename component containing selected residues"
            .type = str
        unselected_name = 'ground-state'
            .help = "Output filename component containing residues not selected through 'rename'"
            .type = str
    }
    pruning {
        prune_duplicates = True
            .help = 'Remove duplicated conformers in the output structure (convert to blank altloc)'
            .type = bool
        rmsd_cutoff = 0.1
            .type = float
    }
    reset_altlocs = True
        .help = 'Relabel conformers of kept residues to begin with "A" (i.e. C,D,E -> A,B,C)'
        .type = bool
    reset_occupancies = True
        .help = 'Normalise the occupancies so that the maximum occupancy of the output structure is 1.0 (all relative occupancies are maintained)'
        .type = bool
}
settings {
    overwrite = False
        .type = bool
    verbose = False
        .type = bool
}

""")

############################################################################

def split_conformations(filename, params, log=None):

    if log is None: log = Log(verbose=True)

    # Read the pdb header - for writing later...
    header_contents = get_pdb_header(filename)

    # Read in and validate the input file
    ens_obj = strip_pdb_to_input(filename, remove_ter=True)
    ens_obj.hierarchy.only_model()

    # Create a new copy of the structures
    new_ens = ens_obj.hierarchy.deep_copy()

    # Extract conformers from the structure as set
    all_confs = set(ens_obj.hierarchy.altloc_indices())
    all_confs.discard('')

    if params.options.mode == 'by_residue_name':
        sel_resnames = params.options.by_residue_name.resname.split(',')
        sel_confs = [ag.altloc for ag in new_ens.atom_groups() if (ag.resname in sel_resnames)]
        # List of conformers to output for each structure, and suffixes
        out_confs = map(sorted, [all_confs.intersection(sel_confs), all_confs.difference(sel_confs)])
        out_suffs = [params.options.by_residue_name.selected_name, params.options.by_residue_name.unselected_name]
    elif params.options.mode == 'by_conformer':
        sel_resnames = None
        sel_confs = None
        # One structure for each conformer
        out_confs = [[c] for c in sorted(all_confs)]
        out_suffs = [''.join(c) for c in out_confs]
    elif params.options.mode == 'by_conformer_group':
        sel_resnames = None
        sel_confs = None
        # One structure for each set of supplied conformer sets
        out_confs = [s.split(',') for s in params.options.by_conformer_group.conformers]
        out_suffs = [''.join(c) for c in out_confs]
    else:
        raise Exception('Invalid selection for options.mode: {}'.format(params.options.mode))

    assert len(out_confs) == len(out_suffs), '{} not same length as {}'.format(str(out_confs), str(out_suffs))

    for confs, suffix in zip(out_confs, out_suffs):
        print 'Conformers {} -> {}'.format(str(confs), suffix)

    # Create paths from the suffixes
    out_paths = ['.'.join([os.path.splitext(filename)[0],params.output.suffix_prefix,suff,'pdb']) for suff in out_suffs]

    log.subheading('Processing {}'.format(filename[-70:]))

    for this_confs, this_path in zip(out_confs, out_paths):

        if not this_confs: continue

        # Select atoms to keep - no altloc, or altloc in selection
        sel_string = ' or '.join(['altid " "']+['altid "{}"'.format(alt) for alt in this_confs])
        # Extract selection from the hierarchy
        sel_hiery = new_ens.select(new_ens.atom_selection_cache().selection(sel_string), copy_atoms=True)

        log.bar(True, False)
        log('Outputting conformer(s) {} to {}'.format(''.join(this_confs), this_path))
        log.bar()
        log('Keeping ANY atom with conformer id: {}'.format(' or '.join(['" "']+this_confs)))
        log('Selection: \n\t'+sel_string)

        if params.options.pruning.prune_duplicates:
            log.bar()
            log('Pruning redundant conformers')
            # Remove an alternate conformers than are duplicated after selection
            prune_redundant_alternate_conformations(
                hierarchy           = sel_hiery,
                required_altlocs    = [a for a in sel_hiery.altloc_indices() if a],
                rmsd_cutoff         = params.options.pruning.rmsd_cutoff,
                in_place            = True,
                verbose             = params.settings.verbose)

        if params.options.reset_altlocs:
            log.bar()
            # Change the altlocs so that they start from "A"
            if len(this_confs) == 1:
                conf_hash = {this_confs[0]: ' '}
            else:
                conf_hash = dict(zip(this_confs, iotbx.pdb.systematic_chain_ids()))
            log('Resetting structure altlocs:')
            for k in sorted(conf_hash.keys()):
                log('\t{} -> "{}"'.format(k, conf_hash[k]))
            if params.settings.verbose: log.bar()
            for ag in sel_hiery.atom_groups():
                if ag.altloc in this_confs:
                    if params.settings.verbose:
                        log('{} -> alt {}'.format(Labeller.format(ag), conf_hash[ag.altloc]))
                    ag.altloc = conf_hash[ag.altloc]

        if params.options.reset_occupancies:
            log.bar()
            log('Resetting output occupancies (maximum occupancy of 1.0, etc.)')
            # Divide through by the smallest occupancy of any complete residues groups with occupancies of less than one
            rg_occs = [calculate_residue_group_occupancy(rg) for rg in residue_groups_with_complete_set_of_conformers(sel_hiery)]
            non_uni = [v for v in numpy.unique(rg_occs) if v < 1.0]
            if non_uni:
                div_occ = min(non_uni)
                log('Dividing all occupancies by {}'.format(div_occ))
                sel_hiery.atoms().set_occ(sel_hiery.atoms().extract_occ() / div_occ)
            # Normalise the occupancies of any residue groups with more than unitary occupancy
            log('Fixing any residues that have greater than unitary occupancy')
            for rg in sel_hiery.residue_groups():
                occ = calculate_residue_group_occupancy(rg)
                if occ > 1.0:
                    log('Normalising residue {} (occupancy {})'.format(Labeller.format(rg), occ))
                    rg.atoms().set_occ(rg.atoms().extract_occ() / occ)
            # Perform checks
            max_occ = max([calculate_residue_group_occupancy(rg) for rg in sel_hiery.residue_groups()])
            log('Maximum occupancy of output structue: {}'.format(max_occ))
            assert max_occ >= 0.0, 'maximum occupancy is less than 0.0?!?!'
            assert max_occ <= 1.0, 'maximum occupancy is greater than 1.0?!?!'

        log.bar()
        log('Writing structure: {}'.format(this_path))
        log.bar(False, True)

        # Write header contents
        with open(this_path, 'w') as fh: fh.write(header_contents)
        # Write output file
        sel_hiery.write_pdb_file(this_path, open_append=True)

    return out_paths

############################################################################

def run(params):

    # Create log file
    log = Log(log_file=params.output.log, verbose=True)

    log.heading('Validating input parameters')

    assert params.input.pdb, 'No PDB files given'

    log.heading('Splitting multi-state structures')

    # Iterate through the input structures and extract the conformation
    for pdb in params.input.pdb:
        split_conformations(filename=pdb, params=params, log=log)

    log.heading('FINISHED')

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
