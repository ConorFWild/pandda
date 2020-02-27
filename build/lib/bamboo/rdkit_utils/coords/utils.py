import os, sys

from bamboo.rdkit_utils.mol import check_pdb_readable

from rdkit import Chem

def get_centroid_from_file(model):
    """Get the centroid of an isolated ligand"""
    return get_centroid_from_mol(Chem.MolFromPDBFile(model))

def get_centroid_from_mol(mol):
    """Get the centroid of a mol"""

    num_atm = mol.GetNumHeavyAtoms()
    # Get the conformer (actual coords)
    conf = mol.GetConformer()
    # Get the atom objects
    atoms = mol.GetAtoms()
    # Zero the sums
    sumx, sumy, sumz, count = (0, 0, 0, 0)

    for a in atoms:
        # Skip Hydrogens (should have been removed anyway, but...)
        if a.GetAtomicNum() == 1:
            continue
        # Get the 3D coords of the atom
        pos = conf.GetAtomPosition(a.GetIdx())
        # Number of atoms used
        count += 1
        # Add the coords to the sum
        sumx += pos.x
        sumy += pos.y
        sumz += pos.z

    # Check we've got as many as we were expecting
    assert count == num_atm, 'Number of atoms used to calculate centroid != number of heavy atoms'

    return (sumx/count,sumy/count,sumz/count)

def calculate_coordinate_differences(model1, model2):
    """Calculate the differences between the atom coordinates of two identical structures"""

    # Read in mols and check for validity
    mol1 = check_pdb_readable(model1)
    mol2 = check_pdb_readable(model2)
    if (not mol1) or (not mol2):
        return None

    # Check that the mols are identical-ish
    if mol1.GetNumHeavyAtoms() != mol2.GetNumHeavyAtoms():
        raise EqualityError('Molecules are not identical (Num Atoms) {!s} != {!s}.\n{!s}\n{!s}'.format(mol1.GetNumHeavyAtoms(),mol2.GetNumHeavyAtoms(),Chem.MolToSmiles(mol1),Chem.MolToSmiles(mol2)))
    if mol1.GetNumBonds() != mol2.GetNumBonds():
        raise EqualityError('Molecules are not identical (Num Bonds) {!s} != {!s}:\n{!s}\n{!s}'.format(mol1.GetNumBonds(),mol2.GetNumBonds(),Chem.MolToSmiles(mol1), Chem.MolToSmiles(mol2)))

    # Gets atoms in mol1 (e.g. 14,5,3...) that match mol2 (1,2,3...)
    matchpatterns = mol1.GetSubstructMatches(mol2, uniquify=False)
    # Check to see if the molecules actually DO contain common substructures
    if not matchpatterns: return None

    differences = []
    # Get the conformers to access the coords
    conf1 = mol1.GetConformer(0)
    conf2 = mol2.GetConformer(0)
    # May be more than one matching pattern. Calculate all of them.
    for matchlist in matchpatterns:
        # reset the vector of coord difference for each match pattern
        match_diffs = []
        # idx2 = 0,1,2... idx1 = 14,5,3...
        for idx2, idx1 in enumerate(matchlist):
            # Get the atom coords
            atm1 = conf1.GetAtomPosition(idx1)
            atm2 = conf2.GetAtomPosition(idx2)
            # Append tuple of differences
            match_diffs.append((atm1.x - atm2.x, atm1.y - atm2.y, atm1.z - atm2.z))
        differences.append(match_diffs)
    # Return the differences corresponding to all of the ways of matching the molecules
    return differences

def get_atomic_equivalences(mol1, mol2):
    """Returns the list of atoms in mol1 that match with the atoms in mol2 => [pattern1, pattern2] ... pattern1 = [(mol1atm, mol2atm), ...]"""

    # Get the match patterns between mol1, mol2
    match_patterns = mol1.GetSubstructMatches(mol2, uniquify=False)

    # List of paired atoms in mol1 -> mol2
    atom_pairings = []

    if mol1.GetNumHeavyAtoms() != mol2.GetNumHeavyAtoms():
        raise EqualityError('Heavy Atom Numbers are not equal!')

    for match_pattern in match_patterns:
        atom_pairings.append(zip(match_pattern,range(0,mol1.GetNumHeavyAtoms())))

    return atom_pairings

