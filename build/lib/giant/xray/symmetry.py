import os

import scipy, numpy

from scitbx.array_family import flex

import iotbx.pdb
from cctbx import crystal, sgtbx

from bamboo.stats.cluster import find_connected_groups

def match_sites_by_symmetry(ref_sites, query_sites, unit_cell, space_group, cutoff=5):
    """Pair sites in query_sites to sites in ref_sites, allowing for symmetry"""

    pairings = numpy.zeros((len(ref_sites), len(query_sites)), dtype=numpy.int)

    for i_ref, ref in enumerate(ref_sites):
        ref_frac = unit_cell.fractionalize(ref)
        sym_sites_ref = sgtbx.sym_equiv_sites(space_group   = space_group,
                                              unit_cell     = unit_cell,
                                              original_site = ref_frac  )

        for i_query, query in enumerate(query_sites):
            query_frac = unit_cell.fractionalize(query)
            min_dist = sgtbx.min_sym_equiv_distance_info(sym_sites_ref, query_frac).dist()
            if min_dist < cutoff:
                pairings[i_ref, i_query] = 1

    return pairings

def combine_hierarchies(list_of_hierarchies):
    """Combine a list of hierarchies into one hierarchy -- Requires all of the chain identifiers to be unique"""
    top_h = list_of_hierarchies[0].deep_copy()
    for next_h in list_of_hierarchies[1:]: top_h.transfer_chains_from_other(next_h.deep_copy())
    return top_h

def find_symmetry_equivalent_groups(points_frac, sym_ops, unit_cell, cutoff_cart):
    """Group sets of points by minimum distance between points, allowing for symmetry"""

    # Pre-calculate the square distance
    dist_cut_sq = cutoff_cart**2
    # Matrix for whether they are near to each other
    equiv_sites = numpy.zeros([len(sym_ops), len(points_frac), len(points_frac)], dtype=int)

    # Apply the symmetry operations to each group to see if it is near to other groups
    for i_sym_op, sym_op in enumerate(sym_ops):
        # Transformed groups under this symmetry operation
        trans_points_frac = [sym_op.as_rational().as_float() * pf for pf in points_frac]
        # Loop through original points
        for i_group_1, group_1 in enumerate(points_frac):
            # Loop through transformed points again
            for i_group_2, group_2 in enumerate(trans_points_frac):
                # Comparing group to itself - skip
                if i_group_1 == i_group_2:
                    equivalent = True
                else:
                    # Start by assuming the clusters are not related
                    equivalent = False
                    # Loop through and see if the clusters overlap
                    for qp in group_2:
                        # Calculate difference and convert back to cartesian
                        diffs_cart = unit_cell.orthogonalize(group_1 - qp)
                        # Check if any are closer than the minimum required
                        if min(diffs_cart.dot()) < dist_cut_sq:
                            equivalent = True
                            break
                # If the clusters overlap, they are equivalent
                if equivalent:
                    equiv_sites[(i_sym_op, i_group_1, i_group_2)] = 1

    # Condense the cluster equivalence - take max over the symmetry operations and group by connected paths
    sym_groups = find_connected_groups(connection_matrix=equiv_sites.max(axis=0))
    return sym_groups

def get_crystal_contact_operators(hierarchy, crystal_symmetry, distance_cutoff):
    """Use an alternate method to identify the symmetry operations required to generate crystal contacts"""

    # Extract the xray structure from the reference hierarchy
    struc = hierarchy.extract_xray_structure(crystal_symmetry=crystal_symmetry)
    atoms = hierarchy.atoms()

    # Extract the mappings that will tell us the adjacent symmetry copies
    asu_mappings = struc.asu_mappings(buffer_thickness=distance_cutoff)
    uc = asu_mappings.unit_cell()
    # Symmetry operations for each atom
    mappings = asu_mappings.mappings()

    # There should be one mappings list per atom
    assert len(struc.scatterers()) == len(mappings)

    # Get all atom pairs within distance_cutoff distance
    pair_generator = crystal.neighbors_fast_pair_generator(asu_mappings, distance_cutoff=distance_cutoff)
    sym_operations = []
    for pair in pair_generator:
      # obtain rt_mx_ji - symmetry operator that should be applied to j-th atom
      # to transfer it to i-th atom
      rt_mx_i = asu_mappings.get_rt_mx_i(pair)
      rt_mx_j = asu_mappings.get_rt_mx_j(pair)
      rt_mx_ji = rt_mx_i.inverse().multiply(rt_mx_j)
      # if it is not a unit matrix, that is symmetry related pair of atoms
      if not rt_mx_ji.is_unit_mx():
        if rt_mx_ji not in sym_operations:
          sym_operations.append(rt_mx_ji)

    return sym_operations

def apply_symmetry_operators(hierarchy, crystal_symmetry, sym_ops_mat):
    """Take a list of symmetry operations and apply them to hierarchy"""

    # Extract the xray structure from the reference hierarchy
    struc = hierarchy.extract_xray_structure(crystal_symmetry=crystal_symmetry)

    # Create a list of chain ids for the new symmetry copies
    # Already existing chain ids
    chains = [c.id for c in hierarchy.chains()]
    # Number of new chains needed - one for each reference chain for each symmetry operation
    num_new_chains = len(chains)*len(sym_ops_mat)
    # Select the new chains, ignoring the ones already in the reference
    new_chain_ids = [c for c in iotbx.pdb.systematic_chain_ids()[0:(num_new_chains+len(chains))] if c not in chains][0:num_new_chains]

    assert not [c for c in new_chain_ids if c in chains], 'GENERATED CHAIN IDS ARE NOT UNIQUE'

    # Create combinations of the different symmetry operations and the reference chains and map them to new chain ids
    sym_op_chain_combinations = []
    for sym_op in sym_ops_mat:
        for r_chain in chains:
            sym_op_chain_combinations.append((sym_op.as_xyz(), r_chain))
    # Dictionary to map symmetry operations and reference chain ids to new chain ids
    new_chain_ids_hash = dict(zip(sym_op_chain_combinations, new_chain_ids))

    # Hierarchies to be returned
    sym_hierarchies = []
    chain_mappings = dict([(c, []) for c in chains])

    for sym_op in sym_ops_mat:

        sym_op_rt = sym_op.as_rational().as_float()

        # Transform all of the coordinates in the reference structure
        transformed_coords = sym_op_rt * struc.sites_frac()

        # Create copy of the xray structure to play with, and set the transformed coordinates
        new_struc = struc.customized_copy()
        new_struc.set_sites_frac(transformed_coords)

        # Transfer the sites to a new hierarchy
        new_hierarchy = hierarchy.deep_copy()
        new_hierarchy.adopt_xray_structure(new_struc)

        # Update the chain ids in the new structure
        for chain in new_hierarchy.chains():
            old_id = chain.id
            chain.id = new_chain_ids_hash[(sym_op.as_xyz(), chain.id)]
            chain_mappings[old_id].append(chain.id)

        # Add the hierarchy to the output dict, referenced by the symmetry operation
        sym_hierarchies.append(new_hierarchy)

    return sym_hierarchies, chain_mappings

def generate_crystal_contacts(hierarchy, crystal_symmetry, distance_cutoff=10):
    """Find symmetry copies of the protein in contact with the asu and generate these copies"""

    sym_ops_mat = get_crystal_contact_operators(hierarchy=hierarchy,crystal_symmetry=crystal_symmetry, distance_cutoff=distance_cutoff)
    sym_hierarchies, chain_mappings = apply_symmetry_operators(hierarchy=hierarchy,crystal_symmetry=crystal_symmetry,sym_ops_mat=sym_ops_mat)
    return sym_ops_mat, sym_hierarchies, chain_mappings

if __name__=='__main__':

    input_file = './reference.pdb'

    inp = iotbx.pdb.input(input_file)
    hie = inp.construct_hierarchy()

    for method in [1,2]:
        sym_ops, contact_mappings, sym_hierarchies, chain_mappings = generate_adjacent_symmetry_copies(    ref_hierarchy=hie,
                                                                                                           crystal_symmetry=inp.crystal_symmetry(),
                                                                                                           buffer_thickness=50,
                                                                                                           method=method)

        print 'CHAIN MAPPINGS:'
        for ch in chain_mappings.keys():
            print '\tCHAIN {!s} maps to {!s}'.format(ch, chain_mappings[ch])
        print 'SYMMETRY HEIRARCHIES:'
        for x in sym_hierarchies:
            print '\t',x
        print 'SYMMETRY OPERATIONS:'
        for x in sorted(sym_ops, key=lambda m: str(m)):
            print '\t',x
        print '{!s} SYMMETRY COPIES GENERATED'.format(len(sym_hierarchies))

        combined_sym_hierarchy = combine_hierarchies(sym_hierarchies)

        output_file = input_file.replace('.pdb', '-contacts-method{!s}.pdb'.format(method))

        assert not os.path.exists(output_file)
        combined_sym_hierarchy.write_pdb_file(output_file)

