from rdkit.Chem import BRICS

def identify_rotatable_bond_atom_pairs(mol):
    """find the atom quadruplets around rotatable bonds of a molecule"""

    # List of tuples of 4 atoms
    atom_sets = []
    # Get the atoms on the ends of rotatable bonds
    atom_pairs = [p[0] for p in BRICS.FindBRICSBonds(mol) if p]
    # Go through and get one of the neighbours for each
    for a1,a2 in atom_pairs:
        # Get the neighbours for a1 (removing a2)
        a1_neighbours = [n.GetIdx() for n in mol.GetAtomWithIdx(a1).GetNeighbors()]
        a1_neighbours.remove(a2)
        # Get the neighbours for a2 (removing a1)
        a2_neighbours = [n.GetIdx() for n in mol.GetAtomWithIdx(a2).GetNeighbors()]
        a2_neighbours.remove(a1)
        # Add one from either side of the double bond
        atom_sets.append((a1_neighbours[0], a1, a2, a2_neighbours[0]))
    # Now have 4 atoms from which we can calculate a dihedral angle
    return atom_sets

