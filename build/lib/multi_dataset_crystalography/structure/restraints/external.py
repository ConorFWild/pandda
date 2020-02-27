
from scitbx.array_family import flex
from giant.structure.select import non_h

def find_atoms_around_alternate_conformers(hierarchy, altlocs=None, dist_cutoff=4.2):
    """For all alternate conformers in (or subset altlocs, if given) return atom pairs to surrounding atoms"""

    # Remove hydrograns and extract atoms
    hierarchy = non_h(hierarchy)
    h_atoms = hierarchy.atoms()

    # Get all the altlocs in the structure
    all_altlocs = list(hierarchy.altloc_indices())
    if not altlocs: altlocs=all_altlocs
    # Get the indices of each conformer in the structure
    conf_indices = hierarchy.get_conformer_indices()
    # Get selection for blank altloc atoms
    i_alt_blank = all_altlocs.index('')
    alt_blank_sel = (conf_indices == i_alt_blank).iselection()

    # Output list and squared distance cutoff
    atom_pairs = []
    dist_cut_sq = dist_cutoff**2

    # Iterate through altlocs
    for alt in altlocs:
        if alt == '':
            continue
        elif alt not in all_altlocs:
            continue
        # Get a selection for atoms with this altloc
        i_alt = all_altlocs.index(alt)
        alt_sel = (conf_indices == i_alt).iselection()
        # Combine with the blank altloc selection
        comb_sel = flex.size_t(sorted(alt_sel.concatenate(alt_blank_sel)))
        # These should be mutually exclusive sets...
        assert len(comb_sel) == len(alt_sel) + len(alt_blank_sel)
        # Extract all atoms of this conformer
        alt_ats = h_atoms.select(alt_sel)
        comb_ats = h_atoms.select(comb_sel)

        # Iterate through the atoms in this conformation
        for atom in alt_ats:
            # Find all atoms within dist_cutoff
            at_dists_sq = (comb_ats.extract_xyz() - atom.xyz).dot()
            at_dists_sel = (at_dists_sq < dist_cut_sq).iselection()
            # Iterate through nearby atoms and append
            for atom_2 in comb_ats.select(at_dists_sel):
                atom_pairs.append((atom.fetch_labels(), atom_2.fetch_labels(), round(atom.distance(atom_2),3)))

    return atom_pairs
