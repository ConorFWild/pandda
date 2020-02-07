import os, sys, copy

import libtbx.phil

import numpy

import iotbx.pdb

from scitbx.array_family import flex

from giant.structure.occupancy import normalise_occupancies, set_conformer_occupancy

#######################################

PROGRAM = 'giant.set_occupancies'
DESCRIPTION = """
    A tool to allow the occupancies of particular conformers in a model, either by
        conformer id   - set the occupancy of conformer explicitly
        residue altloc - set the occupancy of conformer associated with a residue

    1) Simple usage:
        > giant.set_occupancies input.pdb
"""

blank_arg_prepend = {'.pdb':'input='}

master_phil = libtbx.phil.parse("""

conf_ids = None
    .type = str
res_ids = LIG,UNL,FRG
    .type = str

set_occ = 0.1
    .type = float

input = None
    .type = path
output = None
    .type = path

minimum_occ = 0.0
    .type = float
maximum_occ = 1.0
    .type = float

overwrite = False
    .type = bool
verbose = False
    .type = bool
""")

#######################################

def run(params):

    ######################################################################
    if params.verbose: print '===========================================>>>'
    if params.verbose: print 'SETTING OCCUPANCIES'
    ######################################################################

    assert params.input, 'No Input file given!'

    if not params.output: params.output='./occ-edited.pdb'
    if os.path.exists(params.output):
        if params.overwrite: os.remove(params.output)
        else: raise Exception('File already exists: {}'.format(params.output))

    assert params.conf_ids or params.res_ids, 'Must provide conformer IDs or residue IDs'

    ######################################################################
    if params.verbose: print '===========================================>>>'
    if params.verbose: print 'READING INPUT FILES'
    ######################################################################
    # Read in the ligand file and set each residue to the requested conformer
    prot_h = iotbx.pdb.hierarchy.input(params.input).hierarchy

    ######################################################################
    if params.verbose: print '===========================================>>>'
    if params.verbose: print 'VALIDATING INPUT MODELS'
    ######################################################################
    # Check that ... something
    prot_h.only_model()

    ######################################################################
    if params.verbose: print '===========================================>>>'
    if params.verbose: print 'SETTING OCCUPANCIES'
    ######################################################################

    if params.res_ids:  sel_res_ids = params.res_ids.split(',')
    else:               sel_res_ids = []
    if params.conf_ids: sel_conf_ids = params.conf_ids.split(',')
    else:               sel_conf_ids = []
    all_conf_ids = sel_conf_ids + list(set([ag.altloc for ag in prot_h.atom_groups() if ag.altloc and (ag.resname in sel_res_ids)]))

    print 'Setting Occupancy for Conformers: {}'.format(sorted(all_conf_ids))

    for conf in all_conf_ids:
        if params.verbose: print 'Setting Conformer {} to {}'.format(conf, params.set_occ)
        set_conformer_occupancy(hierarchy=prot_h, conf_id=conf, conf_occ=params.set_occ, in_place=True)

    prot_h = normalise_occupancies( hierarchy = prot_h,
                                    exclude_conformers = all_conf_ids,
                                    max_occ = params.maximum_occ,
                                    min_occ = params.minimum_occ  )

    ######################################################################
    if params.verbose: print '===========================================>>>'
    if params.verbose: print 'WRITING OUTPUT STRUCTURE'
    ######################################################################

    # Write output file
    prot_h.write_pdb_file(params.output)

    return

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
