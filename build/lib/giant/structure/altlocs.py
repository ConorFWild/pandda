from itertools import cycle

import numpy
import iotbx.pdb

from scitbx.array_family import flex

from giant.structure import get_atom_pairs, calculate_paired_atom_rmsd, normalise_occupancies
from giant.structure.formatting import Labeller

protein_amino_acid_set = set(iotbx.pdb.common_residue_names_amino_acid + iotbx.pdb.common_residue_names_modified_amino_acid)

#####################################################
###             HIERARCHY FUNCTIONS               ###
#####################################################

def find_next_conformer_idx(hierarchy, all_ids=iotbx.pdb.systematic_chain_ids()):
    current_conf_ids = sorted(hierarchy.altloc_indices())
    for i_id, new_conf_id in enumerate(all_ids):
        if new_conf_id not in current_conf_ids:
            return i_id

def increment_altlocs(hierarchy, offset=1, in_place=False, verbose=False):
    """Increment all altlocs in the structure by a certain number of letters"""
    if not in_place: hierarchy=hierarchy.deep_copy()
    cur_altlocs = [a for a in hierarchy.altloc_indices() if a]
    all_altlocs = iotbx.pdb.systematic_chain_ids()
    new_altlocs = dict([(a, all_altlocs[all_altlocs.index(a)+offset]) for a in cur_altlocs])
    if verbose:
        print '------------------>'
        print 'Updating altlocs:'
        for a in cur_altlocs:
            print '{} -> {}'.format(a, new_altlocs[a])
        print '------------------>'
    for atom_group in hierarchy.atom_groups():
        if verbose: print '{} - updating altloc: {} -> {}'.format(Labeller.format(atom_group), atom_group.altloc, new_altlocs[atom_group.altloc])
        assert atom_group.altloc != ''
        atom_group.altloc = new_altlocs[atom_group.altloc]
    return hierarchy

def expand_alternate_conformations(hierarchy, in_place=False, verbose=False):
    """Convert all atoms to multiple conformers - full multi-conformer representation of the model"""
    if not in_place: hierarchy = hierarchy.deep_copy()
    # Get all of the altlocs that should be present for each atom
    full_altloc_set = sorted([a for a in hierarchy.altloc_indices() if a])
    # If not altlocs found, expand all to "A"
    if full_altloc_set == []:
        if verbose:
            print 'No altlocs in structure: expanding all residues to conformer "A"'
        full_altloc_set = ['A']
    if verbose:
        print 'Expanding all (appropriate) residues to have altlocs {}'.format(full_altloc_set)
        print '------------------>'
    # Iterate through and expand each residue group to have all conformers
    for chain in hierarchy.chains():
        for residue_group in chain.residue_groups():
            # If has conformers but has blank altloc atoms (add blank ag to all other ags)
            if residue_group.have_conformers() and residue_group.move_blank_altloc_atom_groups_to_front():
                if verbose: print '{} - expanding to pure conformer (altlocs {})'.format(Labeller.format(residue_group), [a.altloc for a in residue_group.atom_groups()])
                # Convert all residue_groups to pure alt-conf
                create_pure_alt_conf_from_proper_alt_conf(residue_group=residue_group, in_place=True)
            # Can go to next if all conformers are present for this residue group
            current_set = {a.altloc for a in residue_group.atom_groups()}
            if not current_set.symmetric_difference(full_altloc_set): continue
            # Only want to expand conformers for protein atoms (which should be present in all conformers)
            # or where the residue group is only present in one conformation (single conformer water)
            # but DO NOT want to expand waters in conformer A to A,B,C etc...
            if protein_amino_acid_set.intersection(residue_group.unique_resnames()) or (not residue_group.have_conformers()):
                if verbose: print '{} - populating missing conformers (current altlocs {}, target set {})'.format(Labeller.format(residue_group), current_set, full_altloc_set)
                # Populate missing conformers (from the other conformers)
                populate_missing_conformers(residue_group=residue_group, full_altloc_set=full_altloc_set, in_place=True, verbose=verbose)
                assert sorted([a.altloc for a in residue_group.atom_groups()]) == full_altloc_set
                if verbose: print '{} - updated conformer list: (current altlocs {}, target set {})'.format(Labeller.format(residue_group), [a.altloc for a in residue_group.atom_groups()], full_altloc_set)
    if verbose: print '------------------>'
    return hierarchy

def prune_redundant_alternate_conformations(hierarchy, required_altlocs=[], rmsd_cutoff=0.1, in_place=False, verbose=False):
    """Remove alternate conformers of residues if residues has conformers of required_altlocs and all conformers are within rmsd_cutoff"""
    if not in_place: hierarchy = hierarchy.deep_copy()
    required_altlocs = set(required_altlocs)
    for residue_group in hierarchy.residue_groups():
        # Skip if no conformers
        if not residue_group.have_conformers():
            continue
        # Get the blank and non-blank altloc atom_groups
        if residue_group.move_blank_altloc_atom_groups_to_front() != 0:
            main_ag = residue_group.atom_groups()[0]
            alt_ags = residue_group.atom_groups()[1:]
            assert main_ag.altloc == ''
            assert alt_ags != []
        else:
            main_ag = None
            alt_ags = residue_group.atom_groups()
        # Check no misplaced main conf
        assert '' not in [ag.altloc for ag in alt_ags]
        # Check if required altlocs are present (return if not)
        if required_altlocs.difference([ag.altloc for ag in alt_ags]):
            continue
        # Check if all pair of conformers are within rmsd cutoff
        prune = True
        for i,ag_1 in enumerate(alt_ags):
            for j,ag_2 in enumerate(alt_ags):
                if j<=i: continue
                d = calculate_paired_atom_rmsd(atoms_1=ag_1.atoms(), atoms_2=ag_2.atoms(), sort=True, truncate_to_common_set=False)
                if verbose: print 'Residue {}, alt {} - alt {}: rmsd {}'.format(Labeller.format(residue_group), i,j,d)
                if (d is None) or (d > rmsd_cutoff):
                    prune = False
                    break
            if prune is False: break
        if prune is False: continue
        # All rmsds below cutoff - prune!
        if verbose: print 'Pruning {}: altlocs {} -> [""]'.format(Labeller.format(residue_group), [ag.altloc for ag in alt_ags])
        if main_ag:
            # Merge one alt group with the main atom_group
            new_main_ag = alt_ags[0].detached_copy()
            new_main_ag.altloc = ''
            normalise_occupancies(atoms=new_main_ag.atoms(), max_occ=max(main_ag.atoms().extract_occ()))
            residue_group.merge_atom_groups(main_ag, new_main_ag)
        else:
            # Remove one atom_group and set altloc to ''
            new_main_ag = alt_ags.pop(0)
            new_main_ag.altloc = ''
            normalise_occupancies(atoms=new_main_ag.atoms(), max_occ=sum([max(ag.atoms().extract_occ()) for ag in [new_main_ag]+alt_ags]))
        # Remove all remaining alternate groups
        [residue_group.remove_atom_group(ag) for ag in alt_ags]
        assert len(residue_group.atom_groups())==1

    return hierarchy

#####################################################
###            RESIDUE GROUP FUNCTIONS            ###
#####################################################

def find_duplicated_conformers_and_generate_atom_pairs(residue_group, rmsd_cutoff):
    """Return pairs of equivalent atoms for conformers with an RMSD less than rmsd_cutoff"""

    resi_pairs = find_duplicate_conformers(residue_group=residue_group, rmsd_cutoff=rmsd_cutoff)
    atom_pairs = [get_atom_pairs(residue_1=r1, residue_2=r2, fetch_labels=True) for r1,r2 in resi_pairs]
    return atom_pairs

def find_duplicate_conformers(residue_group, rmsd_cutoff=0.1, include_truncated_atom_groups=True):
    """Return pairs of residue objects for conformers with an RMSD less than rmsd_cutoff"""

    duplicate_conformers = []
    for i, c1 in enumerate(residue_group.conformers()):
        r1 = c1.only_residue()
        a1 = r1.atoms()
        for j, c2 in enumerate(residue_group.conformers()):
            if j<=i:
                continue
            r2 = c2.only_residue()
            if r1.resname != r2.resname:
                continue
            a2 = r2.atoms()
            d = calculate_paired_atom_rmsd(atoms_1=a1, atoms_2=a2, sort=True, truncate_to_common_set=include_truncated_atom_groups)
            if d is None: continue
            if d < rmsd_cutoff:
                duplicate_conformers.append((r1.standalone_copy(),r2.standalone_copy()))

    return duplicate_conformers

def create_pure_alt_conf_from_proper_alt_conf(residue_group, in_place=False):
    """Take 'main conf' atoms and add to 'pure conformer' atom_groups"""
    if not in_place: residue_group = residue_group.detached_copy()
    main_ags = [ag for ag in residue_group.atom_groups() if ag.altloc=='']
    conf_ags = [ag for ag in residue_group.atom_groups() if ag.altloc!='']
    assert len(main_ags)==1, "Must be one atom_group in residue_group with altloc==''"
    assert len(conf_ags)!=0, "Must be at least one alternate conformer present"
    # Atom group to be merged into the other ags
    main_ag = main_ags[0]
    # Remove from the residue group
    residue_group.remove_atom_group(main_ag)
    # Iterate through and add atoms from main_conf to each of the other atom_groups
    for conf_ag in conf_ags:
        new_main_ag = main_ag.detached_copy()
        # Set the occupancy of the main_conf atoms to be transferred
        max_occ = max(conf_ag.atoms().extract_occ())
        new_main_ag.atoms().set_occ(flex.double(new_main_ag.atoms().size(), max_occ))
        # Change the altloc - very important
        new_main_ag.altloc = conf_ag.altloc
        # Merge the atom groups
        residue_group.merge_atom_groups(primary=conf_ag, secondary=new_main_ag)
    return residue_group

def populate_missing_conformers(residue_group, full_altloc_set, in_place=False, verbose=False):
    """
    Expand the alternate conformations of the residue_group to full_set, using the existing conformations
    - no conformers are deleted if existing conformers are not in full_altloc_set
    """
    if not in_place: residue_group = residue_group.detached_copy()
    # Create pool of atom groups to use for the new conformers
    if residue_group.have_conformers():
        # Use any atom_group with an altloc - sort in decreasing occupancy
        pool_ags = sorted([ag for ag in residue_group.atom_groups() if ag.altloc], key=lambda x: max(x.atoms().extract_occ()), reverse=True)
        pool_alts = [ag.altloc for ag in pool_ags]
        # Occupancy multipliers for each atom group dependant on how much they'll be duplicated
        num_cur = len(pool_alts)
        num_tar = num_cur + len(set(full_altloc_set).difference(pool_alts))
        occ_cor = [1.0/((num_tar//num_cur)+((num_tar%num_cur)>k)) for k in range(num_cur)]
    else:
        # Only one conformation
        pool_ags = [residue_group.only_atom_group()]
        pool_alts = []
        # Occupancy multipliers for each atom group dependant on how much they'll be duplicated
        occ_cor = [1.0/len(full_altloc_set)]
        # Remove the blank conformer (to add it again with different altlocs later)
        assert pool_ags[0].altloc == ''
        residue_group.remove_atom_group(pool_ags[0])
    if verbose: print 'Occupancy Corrections:', occ_cor
    # Apply occupancy multipliers
    for mult, ag in zip(occ_cor,pool_ags):
        ag.atoms().set_occ(ag.atoms().extract_occ()*mult)
    # Create cycle to iterate through ags as needed
    ag_cycle = cycle(pool_ags)
    # Iterate through and create alternate conformers as required
    for altloc in full_altloc_set:
        # Check if altloc already present
        if altloc in pool_alts: continue
        # Create copy and add to residue_group
        new_ag = ag_cycle.next().detached_copy()
        new_ag.altloc = altloc
        residue_group.append_atom_group(new_ag)
    return residue_group

