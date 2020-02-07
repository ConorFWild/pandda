import os, sys, copy

import libtbx.phil

import iotbx.pdb
import iotbx.ccp4_map

import cctbx.maptbx
import cctbx.sgtbx
import iotbx.map_tools

from scitbx.array_family import flex

from giant.structure.align import align_chains

from mmtbx.maps.superpose import mask_grid, generate_p1_box

#######################################

blank_arg_prepend = None

master_phil = libtbx.phil.parse("""
align {
    reference = None
        .type = path
    chain = None
        .type = str
        .help = 'The chain to be aligned to'
    align_all_chains = False
        .type = bool
        .help = 'Align all chains in pdb to reference'
}
input {
    align_pdb = None
        .type = path
        .help = 'The pdb to use for alignment'
    pdb = None
        .type = path
        .multiple = True
    map = None
        .type = path
        .multiple = True
}
output {
    out_dir = ./
        .type = path
    pdb_suffix = '.aligned.pdb'
        .type = str
    pdb_prefix = ''
        .type = str
    map_suffix = '.aligned.map'
        .type = str
    map_prefix = ''
        .type = str
}
""")

#######################################

def run(params):

    assert params.align.reference, 'A REFERENCE PDB MUST BE PROVIDED'
    assert params.input.pdb or params.input.align_pdb, 'A PDB FILE TO BE ALIGNED MUST BE PROVIDED'
    if not params.input.align_pdb:
        params.input.align_pdb = params.input.pdb[0]
        print 'SELECTING FIRST PDB FOR ALIGNING: {}'.format(params.input.align_pdb)

    # Read in the reference structure
    ref_inp = iotbx.pdb.hierarchy.input(params.align.reference)
    if params.align.chain:
        ref_chain = [c for c in ref_inp.hierarchy.chains() if c.id == params.align.chain][0]
    else:
        ref_chain = [c for c in ref_inp.hierarchy.chains()][0]

    print 'ALIGNING TO CHAIN {!s} OF THE REFERENCE STRUCTURE'.format(ref_chain.id)

    # Read in the structure to be aligned
    align_inp = iotbx.pdb.hierarchy.input(params.input.align_pdb)

    # Build up a list of rotation matrices and write out transformed maps
    for align_chain in align_inp.hierarchy.chains():
        if not align_chain.is_protein():
            continue
        if (not params.align.align_all_chains) and (not align_chain.id==ref_chain.id):
            continue

        print 'GENERATING ALIGNMENT'

        # Align to the reference chain
        lsq = align_chains(mov_chain=align_chain, ref_chain=ref_chain)
        # Extract the rotation matrix
        rt = lsq.rt()

        # Loop through the pdbs and transform each one
        for pdb_file in params.input.pdb:

            print 'TRANSFORMING:', pdb_file
            # Read in the structure to be aligned
            mov_inp = iotbx.pdb.hierarchy.input(pdb_file)
            # Generate output filename
            output_pdb = os.path.join(params.output.out_dir, params.output.pdb_prefix + os.path.basename(pdb_file) + '.chain{!s}'.format(align_chain.id) + params.output.pdb_suffix)
            # Transform coordinates and write
            mov_copy = mov_inp.hierarchy.deep_copy()
            mov_copy.atoms().set_xyz(rt*mov_copy.atoms().extract_xyz())
            mov_copy.write_pdb_file(output_pdb)

        # Loop through the maps and transform each one
        for map_file in params.input.map:

            print 'TRANSFORMING:', map_file
            # Read in the map to be aligned
            mov_map = iotbx.ccp4_map.map_reader(map_file)
            # Unpad the input map
            gridding_first = mov_map.data.origin()
            gridding_last = tuple(flex.int(mov_map.data.all())-2)
            mov_data = cctbx.maptbx.copy(mov_map.data.as_double(), first=gridding_first, last=gridding_last)
            # Generate output filename
            output_map = os.path.join(params.output.out_dir, params.output.map_prefix + os.path.basename(map_file) + '.chain{!s}'.format(align_chain.id) + params.output.map_suffix)
            assert not os.path.exists(output_map), 'OUTPUT MAP ALREADY EXISTS: {!s}'.format(output_map)
            # Rotate the map
            rot_map_data = cctbx.maptbx.rotate_translate_map(
                unit_cell=mov_map.unit_cell(),
    #            map_data=mov_map.data.as_double(),
                map_data=mov_data,
                rotation_matrix=rt.inverse().r.elems,
                translation_vector=rt.inverse().t.elems)
            # Write the map
            iotbx.ccp4_map.write_ccp4_map(
                file_name=output_map,
                unit_cell=mov_map.unit_cell(),
                space_group=cctbx.sgtbx.space_group_info(number=mov_map.space_group_number).group(),
                map_data=rot_map_data,
                labels = flex.std_string(['rotated map. aligned chain {!s} to chain {!s}.'.format(align_chain.id, ref_chain.id)])  )

#        # New change map
#        fake_symm = generate_p1_box(ref_inp.hierarchy, buffer=10.0)
#        xray_structure = ref_inp.hierarchy.extract_xray_structure(crystal_symmetry=fake_symm)
#        f_calc = xray_structure.structure_factors(d_min=2).f_calc()
#        fake_map = f_calc.fft_map(resolution_factor=0.33)
#
#        map_data_superposed = cctbx.maptbx.superpose_maps(
#            unit_cell_1        = mov_map.unit_cell(),
#            unit_cell_2        = fake_symm.unit_cell(),
#            map_data_1         = mov_map.data.as_double(),
#            n_real_2           = fake_map.n_real(),
#            rotation_matrix    = rt.inverse().r.elems,
#            translation_vector = rt.inverse().t.elems)
#
#        map_data_superposed = mask_grid(
#            xrs      = xray_structure,
#            buffer   = 10,
#            map_data = map_data_superposed,
#            n_real   = fake_map.n_real())
#
#        iotbx.map_tools.write_ccp4_map(
#            sites_cart=xray_structure.sites_cart(),
#            unit_cell=fake_symm.unit_cell(),
#            map_data=map_data_superposed,
#            n_real=fake_map.n_real(),
#            file_name=output_map,
#            buffer=10)

#######################################

if __name__=='__main__':
    from giant.jiffies import run_default
    run_default(run=run, master_phil=master_phil, args=sys.argv[1:], blank_arg_prepend=blank_arg_prepend)
