import os, sys, copy, re
import time

import libtbx.phil

import numpy

from scitbx.array_family import flex
from libtbx.utils import Sorry, Failure

#######################################

blank_arg_prepend = {'.pdb:','query_pdb=', '.mtz:','query_mtz='}

master_phil = libtbx.phil.parse("""
input {
    query_pdb = None
        .type = path
        .multiple = True
    query_mtz = None
        .type = path
        .multiple = True
    ref_pdb = None
        .type = path
        .multiple = True
    ref_mtz = None
        .type = path
        .multiple = True
}
output {
    prefix = ''
        .type = str
    suffix = '.aligned'
        .type = str
}
verbose = True
""")

#######################################

def run(params):

    assert params.input.reference_pdb is not None, 'No reference structure provided: input.reference_pdb='
    assert params.input.pdb, 'No structures provided for alignment: input.pdb='
    for p in params.input.pdb:
        if not os.path.exists(p): raise Sorry('Input file "{}" does not exists'.format(p))

    if params.alignment.conformer:
        assert isinstance(params.alignment.conformer, str), 'params.alignment.conformer must be one letter or None'
        assert len(params.alignment.conformer) < 2, 'Can only supply one conformer for alignment, not {}'.format(params.alignment.conformer)
        confs_to_align = ['', params.alignment.conformer]
    else:
        confs_to_align = ['']

    print pretty_string('===========================>').bold()
    print pretty_string('giant.local_align').bold()
    print pretty_string('===========================>').bold()
    print pretty_string('Aligning {} structure(s) to {}'.format(len(params.input.pdb),params.input.reference_pdb))
    print pretty_string('Aligning conformers: {}'.format('"'+'" and "'.join(confs_to_align)+'"'))
    print pretty_string('===========================>').bold()

    t_start = time.time()

    # =================================================>
    # Load structure
    # =================================================>
    h_ref = iotbx.pdb.hierarchy.input(params.input.reference_pdb).hierarchy

    for p_mov in params.input.pdb:

        # =================================================>
        # Load structure
        # =================================================>
        h_mov = iotbx.pdb.hierarchy.input(p_mov).hierarchy

        print pretty_string('\n===========================>').blue()
        print pretty_string('Aligning {} to {}'.format(p_mov,params.input.reference_pdb)).blue()
        print pretty_string('===========================>').blue()

        combined_alignment = align_structures_flexible(mov_hierarchy=h_mov, ref_hierarchy=h_ref,
                                conformers=confs_to_align, cutoff_radius=params.alignment.cutoff_radius,
                                sequence_identity_threshold=params.alignment.sequence_identity_threshold,
                                one_to_one_mapping=params.alignment.one_to_one_chain_alignment,
                                require_hierarchies_identical=params.alignment.similar_structures_only,
                                verbose=params.verbose)

        # =================================================>
        # Transform coordinates
        # =================================================>
        aligned_coords = combined_alignment.nat2ref(coordinates=h_mov.atoms().extract_xyz())
        if params.verbose:
            rmsd = (h_mov.atoms().extract_xyz() - aligned_coords).rms_length()
            print 'All-atom RMSD (before v after alignment): {}'.format(round(rmsd,3))

        # =================================================>
        # Output structure
        # =================================================>
        h_mov_new = h_mov.deep_copy()
        h_mov_new.atoms().set_xyz(aligned_coords)
        h_mov_new_fname = os.path.join(os.path.dirname(p_mov), params.output.prefix+params.output.suffix.join(os.path.splitext(os.path.basename(p_mov))))
        assert os.path.abspath(h_mov_new_fname) != os.path.abspath(p_mov), 'output filename is same as input filename'
        print 'Saving output structure to {}'.format(h_mov_new_fname)
        h_mov_new.write_pdb_file(file_name=h_mov_new_fname)

    t_end = time.time()
    print pretty_string('\n===========================>').bold()
    print pretty_string('Finished!').bold()
    print pretty_string('Total Runtime: {!s}'.format(time.strftime("%H hours:%M minutes:%S seconds", time.gmtime(t_end-t_start)))).bold()
    print pretty_string('Avg/Structure: {!s}'.format(time.strftime("%H hours:%M minutes:%S seconds", time.gmtime((t_end-t_start)/len(params.input.pdb))))).bold()
    print pretty_string('===========================>').bold()

#######################################

if __name__ == '__main__':
    from giant.jiffies import run_default
    run_default(run=run, master_phil=master_phil, args=sys.argv[1:], blank_arg_prepend=blank_arg_prepend)
