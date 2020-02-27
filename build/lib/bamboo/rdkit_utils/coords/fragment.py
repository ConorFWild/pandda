import os, sys
import math

from bamboo.rdkit_utils.mol import check_pdb_readable
from bamboo.rdkit_utils.coords.utils import get_atomic_equivalences
from bamboo.ccp4_utils import map_to_reference_using_symmetry

def calculate_rmsd_between_fragged_files_using_symmetry(file1, file2, delete_structures=True):
    """Take a pair of homologous fragmented molecule files and return a list of the rmsds between the paired fragments in each molecule, minimising over equivalent models"""

    # Map file2 to file1 using symmetry
    tempfile = map_to_reference_using_symmetry(refpdb=file1, movpdb=file2, pdbout=file2.replace('.pdb','.temp'))
    # Calculate the centroid as normal
    rmsds = calculate_rmsd_between_fragged_files(file1, tempfile)
    # Remove tempfiles and return
    if delete_structures: os.remove(tempfile)
    return rmsds

def calculate_rmsd_between_fragged_files(file1, file2):
    """Take a pair of homologous fragmented molecule files and return a list of the rmsds between the paired fragments in each molecule"""

    mol1 = check_pdb_readable(file1)
    mol2 = check_pdb_readable(file2)

    return calculate_rmsd_between_fragged_mols(mol1, mol2)

def calculate_rmsd_between_fragged_mols(mol1, mol2):
    """Take a pair of homologous fragmented molecules mols and return a list of the rmsds between the paired fragments in each molecule"""

    # TODO Extend this so it loops over conformers
    conf1 = mol1.GetConformer()
    conf2 = mol2.GetConformer()

    # Get the atom pairings between mol1, mol2 [[(mol1atm,mol2atm),...],...]
    try:
        match_patterns = get_atomic_equivalences(mol1, mol2)
    except EqualityError:
        raise

    # Assert that there is only one residue present (albeit with different insertion codes)
    assert len(list(set([at.GetMonomerInfo().GetChainId() for at in mol1.GetAtoms()])))==1, 'MORE THAN ONE CHAIN PRESENT IN MOL!'
    assert len(list(set([at.GetMonomerInfo().GetResidueName() for at in mol1.GetAtoms()])))==1, 'MORE THAN ONE RESIDUE NAME PRESENT IN MOL!'
    assert len(list(set([at.GetMonomerInfo().GetResidueNumber() for at in mol1.GetAtoms()])))==1, 'MORE THAN ONE RESIDUE NUMBER PRESENT IN MOL!'

    # Collect the fragments for mol1
    mol1_inscodes = sorted(list(set([at.GetMonomerInfo().GetInsertionCode() for at in mol1.GetAtoms()])))

    match_pattern_fragment_rmsds = []

    # Iterate through the different equivalent ways of matching the atoms between the molecules
    for match_pat in match_patterns:
        # Enable us to translate between mol1 and mol2
        match_dict = dict(match_pat)

        # RMSDs between fragments for this matching patterns
        fragment_rmsds = []
        total_dev_sq = []
        # List of the equivalent inscodes in mol2
        mol2_inscodes = []

        # Iterate through and form and compare the fragments for mol1 and mol2
        for mol1_inscode in mol1_inscodes:
            # Get all of the atoms in mol1 with the inscode == `inscode`
            mol1_frag_atms = [at for at in mol1.GetAtoms() if at.GetMonomerInfo().GetInsertionCode() == mol1_inscode]
            # Get the equivalent atoms in mol2
            mol2_frag_atms = [mol2.GetAtomWithIdx(match_dict[at.GetIdx()]) for at in mol1_frag_atms]

            # Get the equivalent mol2 inscode for this mol1 inscode
            mol2_inscodes.append(list(set([at.GetMonomerInfo().GetInsertionCode() for at in mol2_frag_atms]))[0])

            # Sum of length sqs of vectors joining atoms
            dists_sq = []
            for at1, at2 in zip(mol1_frag_atms, mol2_frag_atms):
                p1 = conf1.GetAtomPosition(at1.GetIdx())
                p2 = conf2.GetAtomPosition(at2.GetIdx())
                # Distance between points
                d = p1.Distance(p2)
                dists_sq.append(d**2)

            # Calculate rmsd for fragment and append
            rmsd = math.sqrt(sum(dists_sq)/len(dists_sq))
            fragment_rmsds.append(rmsd)

            # Keep a running total to enable a global rmsd to be calculated
            total_dev_sq.extend(dists_sq)

        # Calculate whole-molecule rmsd to allow for comparison between match_patterns
        total_rmsd = math.sqrt(sum(total_dev_sq)/len(total_dev_sq))

        # Record the fragment mappings, the fragment rmsds and the global rmsd
        match_pattern_fragment_rmsds.append((mol1_inscodes, mol2_inscodes, fragment_rmsds, total_rmsd))

    if not match_pattern_fragment_rmsds:
        raise Exception('No RMSDS calculated between fragments!')

    # Find the minimum of the total rmsds to find the optimum mapping
    min_rmsd = min([r[3] for r in match_pattern_fragment_rmsds])

    # Get those associated scores
    minimising_match_pattern = [r for r in match_pattern_fragment_rmsds if r[3] == min_rmsd][0]

    return zip(*minimising_match_pattern[0:3])

