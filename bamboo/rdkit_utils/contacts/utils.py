from rdkit import Chem

from bamboo.constants import WATER_NAMES
from bamboo.rdkit_utils.mol import check_pdb_readable

def order_structures_by_minimum_distance_to_reference(refpdb, pdbs, mincutoff=1.5, maxcutoff=6):
    """Return the molecules have at least one non-H atom within a certain distance of the reference structure"""

    # Deal with cases where mincutoff, maxcutoff are not specified
    if not mincutoff:
        mincutoff = 0
    if not maxcutoff:
        maxcutoff = 9999999

    pdbs_to_return = []

    refmol = check_pdb_readable(refpdb)

    for pdbfile in pdbs:
        pdbmol = check_pdb_readable(pdbfile)

        min_dist = calculate_minimum_distance_between_mols(refmol, pdbmol)

        # Reject if there is a clash
        if min_dist >= mincutoff and min_dist <= maxcutoff:
            pdbs_to_return.append((min_dist,pdbfile))

    # Order by minimum distance
    sorted_files = sorted(pdbs_to_return, key=lambda tup: tup[0])

    return [t[1] for t in sorted_files]

def order_structures_by_number_of_contacts_to_reference(refpdb, pdbs, cutoff=3.5):
    """Calculate the number of contacts between `pdb` and refpdb, and return a list of decreasing contacts"""

    refmol = check_pdb_readable(refpdb)

    contacts_list = []

    for pdbfile in pdbs:
        pdbmol = check_pdb_readable(pdbfile)

        # Add the number of contacts and the pdbfile to the list
        contacts_list.append((calculate_number_of_pairwise_distances_within_cutoff(refmol, pdbmol, cutoff=cutoff), pdbfile))

    return sorted(contacts_list, key=lambda tup: tup[0], reverse=True)

def calculate_minimum_distance_between_mols(mol1, mol2):
    """Takes two molecules and returns the minimum distance between the two structures"""

    distances = calculate_pairwise_distances_between_mols(mol1, mol2)

    if not distances:
        return None
    else:
        return min(distances)

def calculate_number_of_pairwise_distances_within_cutoff(mol1, mol2, cutoff=3.5):
    """Count number of times a distance less than `cutoff` is observed between an atom of mol1 and an atom of mol2"""

    distances = calculate_pairwise_distances_between_mols(mol1, mol2)

    if not distances:
        return None

    within_cutoff = [(d<cutoff) for d in distances]
    counts = within_cutoff.count(True)

    return counts

def calculate_pairwise_distances_between_mols(mol1, mol2):
    """Calculates pairwise distances between atoms in mol1 and mol2"""

    # TODO Change this so that it runs through all of the conformations...
    conf1 = mol1.GetConformer()
    conf2 = mol2.GetConformer()

    # Get the indexes of the atoms that aren't water
    mol1_idxs = [a.GetIdx() for a in mol1.GetAtoms() if a.GetMonomerInfo().GetResidueName() not in WATER_NAMES]
    mol2_idxs = [a.GetIdx() for a in mol2.GetAtoms() if a.GetMonomerInfo().GetResidueName() not in WATER_NAMES]

    # Iterate through all pairs of atoms and calculate pairwise distances
    distances = [conf1.GetAtomPosition(i1).Distance(conf2.GetAtomPosition(i2)) for i2 in mol2_idxs for i1 in mol1_idxs]

    return distances
