import math
from rdkit import Chem

from bamboo.rdkit_utils.bonds.rotate import identify_rotatable_bond_atom_pairs

def calculate_dihedral_angle_differences(mol1, mol2):
    """Calculates the dihedral angle differences between rotatable bonds of mol1 and mol2"""

    # Get the dihedrals to calculate for both of the molecules (possibly multiple ways of overlaying the mols if symmetry exists)
    atom_sets = calculate_dihedral_atom_equivalences(mol1, mol2)
    # list of possible rmsds for the molecule
    differences = []
    # Iterate through and calculate the rmsd for each set of atom equivalences
    for mol1_atom_set, mol2_atom_set in atom_sets:
        # Calculate the dihedrals of both
        mol1_dihedrals = calculate_dihedral_angles(mol1, mol1_atom_set)
        mol2_dihedrals = calculate_dihedral_angles(mol2, mol2_atom_set)
        # Calculate the differences squared for each angle difference
        diffs = [an1-an2 for an1, an2 in zip(mol1_dihedrals,mol2_dihedrals)]
        # Append list of angle differences
        differences.append(diffs)

    return atom_sets, differences

def calculate_dihedral_atom_equivalences(mol1, mol2):
    """Gets the list that are paired in mol1 and mol2"""

    # Check that the mols are identical-ish
    if mol1.GetNumHeavyAtoms() != mol2.GetNumHeavyAtoms():
        raise EqualityError('Molecules are not identical (Num Atoms) {!s} != {!s}.\n{!s}\n{!s}'.format(mol1.GetNumHeavyAtoms(),mol2.GetNumHeavyAtoms(),Chem.MolToSmiles(mol1),Chem.MolToSmiles(mol2)))
    if mol1.GetNumBonds() != mol2.GetNumBonds():
        raise EqualityError('Molecules are not identical (Num Bonds) {!s} != {!s}:\n{!s}\n{!s}'.format(mol1.GetNumBonds(),mol2.GetNumBonds(),Chem.MolToSmiles(mol1), Chem.MolToSmiles(mol2)))

    # Gets a list of lists of atoms in mol1 (12,16,3, ...) that match the atoms in mol2 (1,2,3, ...)
    match_patterns = mol1.GetSubstructMatches(mol2, uniquify=False)
    # Get the quadruplets to calculate the dihedrals from for mol1
    mol1_atom_sets = identify_rotatable_bond_atom_pairs(mol1)
    num_atms = mol1.GetNumHeavyAtoms()
    # List for returning
    paired_atom_sets = []
    # Iterate through the different ways of overlaying the molecule (ensures we get the minimum rmsd)
    for match_pattern in match_patterns:
        # Translate from the atoms in mol1 to the atoms in mol2 (for this match_pattern)
        trans_dict = dict(zip(match_pattern, range(0,num_atms)))
        # Translate the atoms in mol1 to the atoms in mol2
        mol2_atom_sets = [ tuple([trans_dict[atm] for atm in bond_set]) for bond_set in mol1_atom_sets]
        # Add to list
        paired_atom_sets.append((mol1_atom_sets, mol2_atom_sets))
        # Check that the atom types are identical (test)
        mol1_atom_types = [ tuple([mol1.GetAtomWithIdx(atm).GetAtomicNum() for atm in bond_set]) for bond_set in mol1_atom_sets]
        mol2_atom_types = [ tuple([mol2.GetAtomWithIdx(atm).GetAtomicNum() for atm in bond_set]) for bond_set in mol2_atom_sets]
        assert mol1_atom_types == mol2_atom_types, "ATOM TYPES ARE NOT THE SAME ON THE DIHEDRAL ANGLE TO BE CALCULATED - THERE'S BEEN A MATCHING ERROR"
    # Return the list of lists of paired atoms between the structures
    return paired_atom_sets

def calculate_dihedral_angles(mol, dihedral_atom_sets):
    """find the dihedral angles of rotatable bonds in a molecule"""

    # Create list for the dihedrals (to be ordered in the same order as the input dihedral sets)
    dihedral_angles = []
    # Now calculate the dihedral angles between the sets identified previously
    conf = mol.GetConformer()
    # Loop through the angles => 2-3 is the rotatable bonds, 1,4 are the neighbours of 2,3 respectively
    for at1, at2, at3, at4 in dihedral_atom_sets:
        # Get the coordinates of the positions
        pos1 = conf.GetAtomPosition(at1)
        pos2 = conf.GetAtomPosition(at2)
        pos3 = conf.GetAtomPosition(at3)
        pos4 = conf.GetAtomPosition(at4)
        # Need to calculate three vectors 1->2, 2->3, 3->4
        vec1 = pos2 - pos1
        vec2 = pos3 - pos2
        vec3 = pos4 - pos3
        # Get the normals to the two planes (vec1-vec2 plane and vec2-vec3 plane))
        cross12 = vec1.CrossProduct(vec2)
        cross23 = vec2.CrossProduct(vec3)
        # Normalise the normals
        cross12.Normalize()
        cross23.Normalize()
        # Calculate dot-product and then inverse cosine to get the angle
        dot_prod = cross12.DotProduct(cross23)
        dihedral_rad = math.acos(dot_prod)
        dihedral_deg = 180*dihedral_rad/math.pi
        dihedral_angles.append(dihedral_deg)
    return dihedral_angles

