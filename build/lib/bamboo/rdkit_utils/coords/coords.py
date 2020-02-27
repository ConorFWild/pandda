import os, sys
import math

from bamboo.rdkit_utils.coords.utils import get_centroid_from_file, calculate_coordinate_differences
from bamboo.ccp4_utils import map_to_reference_using_symmetry

from rdkit import Chem

def calculate_fraction_of_atom_pairs_within_tolerance(model1, model2, tolerance=1, suppress=False):
    """Calculate the percentage of atoms in model1 that are within `tolerance` of the corresponding atom in model2"""

    try: differences = calculate_coordinate_differences(model1, model2)
    except EqualityError as err:
        if suppress:
            print '{!s}: {!s}'.format(type(err),str(err))
            return None
        else: raise
    # Check for failures
    if not differences:
        if suppress: return None
        else:        raise Exception('Error calculating coordinate differences.\n\t{!s}\n\t{!s}'.format(model1,model2))
    # Calculate the lengths and whether that length is within tolerance
    within_tolerance_all = []
    for diffset in differences:
        within = sum([1.0 for d in diffset if (d[0]**2 + d[1]**2 + d[2]**2 < tolerance**2)])
        within_tolerance_all.append(within/len(diffset))
    # Return the rmsd corresponding to the minimal matching pattern
    return max(within_tolerance_all)

def calculate_fraction_of_atom_pairs_within_tolerance_using_symmetry(model1, model2, tolerance=1, suppress=False, delete_structures=True):
    """Calculate the percentage of atoms in model1 that are within `tolerance` of the corresponding atom in model2, maximising over equivalent models"""

    # Map model2 to model1 using symmetry
    tempfile = map_to_reference_using_symmetry(refpdb=model1, movpdb=model2, pdbout=model2.replace('.pdb','.temp'))
    # Calculate the RMSD as normal
    percentage = calculate_fraction_of_atom_pairs_within_tolerance(model1, tempfile, tolerance, suppress)
    # Remove tempfiles and return
    if delete_structures: os.remove(tempfile)
    return percentage

def calculate_rmsd(model1, model2, suppress=False):
    """Calculates the standard Root Mean Squared Deviation between the coordinates of the models"""

    try: differences = calculate_coordinate_differences(model1, model2)
    except EqualityError as err:
        if suppress:
            print '{!s}: {!s}'.format(type(err),str(err))
            return None
        else: raise
    # Check for failures
    if not differences:
        if suppress: return None
        else:        raise Exception('Error calculating coordinate differences.\n\t{!s}\n\t{!s}'.format(model1,model2))
    # Calculate the rms of the differences
    RMSDs = []
    for diffset in differences:
        sqrsum = sum([d[0]**2 + d[1]**2 + d[2]**2 for d in diffset])
        RMSDs.append(math.sqrt(sqrsum/len(diffset)))
    # Return the rmsd corresponding to the minimal matching pattern
    return min(RMSDs)

def calculate_rmsd_using_symmetry(model1, model2, suppress=False, delete_structures=True):
    """Calculates the standard Root Mean Squared Deviation between the coordinates of the models, minimising over equivalent models"""

    # Map model2 to model1 using symmetry
    tempfile = map_to_reference_using_symmetry(refpdb=model1, movpdb=model2, pdbout=model2.replace('.pdb','.temp'))
    # Calculate the RMSD as normal
    RMSD = calculate_rmsd(model1, tempfile, suppress)
    # Remove tempfiles and return
    if delete_structures: os.remove(tempfile)
    return RMSD

def calculate_centroid_difference(model1, model2):
    """Calculates the length of the vector connecting the centroids of the two models"""

    centroid1 = get_centroid_from_file(model1)
    centroid2 = get_centroid_from_file(model2)
    diff_vect = [a-b for a,b in zip(centroid1,centroid2)]
    len_diff_vect = math.sqrt(sum([c**2 for c in diff_vect]))
    return len_diff_vect

def calculate_centroid_difference_using_symmetry(model1, model2, delete_structures=True):
    """Calculates the length of the vector connecting the centroids of the two models, minimising over equivalent models"""

    # Map model2 to model1 using symmetry
    tempfile = map_to_reference_using_symmetry(refpdb=model1, movpdb=model2, pdbout=model2.replace('.pdb','.temp'))
    # Calculate the centroid as normal
    cent = calculate_centroid_difference(model1, tempfile)
    # Remove tempfiles and return
    if delete_structures: os.remove(tempfile)
    return cent

