import numpy

from scitbx.array_family import flex

from giant.structure import get_atom_pairs

def find_duplicated_conformers_and_generate_atom_pairs(residue_group, rmsd_cutoff):
    """Return pairs of equivalent atoms for conformers with an RMSD less than rmsd_cutoff"""

    resi_pairs = find_duplicate_conformers(residue_group=residue_group, rmsd_cutoff=rmsd_cutoff)
    atom_pairs = [get_atom_pairs(residue_1=r1, residue_2=r2, fetch_labels=True) for r1,r2 in resi_pairs]
    return atom_pairs

def find_duplicate_conformers(residue_group, rmsd_cutoff=0.1):
    """Return pairs of residue objects for conformers with an RMSD less than rmsd_cutoff"""

    rmsd_cutoff_sq = rmsd_cutoff**2
    duplicate_conformers = []
    for i, c1 in enumerate(residue_group.conformers()):
        r1 = c1.only_residue()
        a1 = r1.atoms()
        a1_nam = list(a1.extract_name())
        for j, c2 in enumerate(residue_group.conformers()):
            if j<=i:
                continue
            # Extract residue and skip if not comparable
            r2 = c2.only_residue()
            if r1.resname != r2.resname:
                continue
            # Extract atoms
            a2 = r2.atoms()
            a2_nam = list(a2.extract_name())
            # Get atom overlap between conformers
            common_atoms = list(set(a1_nam).intersection(a2_nam))
            # Sort the atoms so can use unsorted residues
            a1_sel = flex.size_t([a1_nam.index(an) for an in common_atoms])
            a2_sel = flex.size_t([a2_nam.index(an) for an in common_atoms])
            # Check selection working as it should
            assert a1.extract_name().select(a1_sel) == a2.extract_name().select(a2_sel)
            # Extract ordered coordinates
            a1_xyz = a1.extract_xyz().select(a1_sel)
            a2_xyz = a2.extract_xyz().select(a2_sel)
            # Claculate RMSD and check below threshold
            d = flex.mean((a1_xyz - a2_xyz).dot())
            if d < rmsd_cutoff_sq: duplicate_conformers.append((r1.standalone_copy(),r2.standalone_copy()))

    return duplicate_conformers

