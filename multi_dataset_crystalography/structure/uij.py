
import math, numpy

from scitbx.array_family import flex
from scitbx import linalg

EIGHT_PI_SQ = 8*math.pi*math.pi

def _reshape_uij(vals):
    vals = numpy.array(vals)
    if vals.shape == (6,):
        single = True
        vals = vals.reshape((1,6))
    else:
        single = False
    assert len(vals.shape) == 2
    assert vals.shape[1] == 6
    return vals, single

def _revert_output(vals, single):
    if single is True:
        return vals[0]
    return vals

def uij_eigenvalues(uij):
    assert len(uij) == 6
    return linalg.eigensystem_real_symmetric(uij).values()

def uij_positive_are_semi_definite(uij, tol=1e-6):
    """Calculate eigenvalues for each uij and check that all are greater than zero (within tolerance)"""
    # Convert to list if only one uij
    uij, single = _reshape_uij(vals=uij)
    # Check tolerance negative
    tol = -1.0 * abs(tol)
    # Extract eigenvalues for each atom
    eigenvalues = numpy.apply_along_axis(uij_eigenvalues, 1, uij)
    # Check all greater than zero
    neg = (eigenvalues < tol)
    out = neg.sum(axis=1).astype(bool)
    # Convert back to value or return list as appropriate
    return _revert_output(vals=out, single=single)

def uij_to_b(uij):
    """Convert anistropic uijs to isotropic B"""
    uij, single = _reshape_uij(vals=uij)
    out = EIGHT_PI_SQ*numpy.mean(uij[:,0:3],axis=1)
    return _revert_output(vals=out, single=single)

def anistropic_uij_to_isotropic_uij(uij):
    """Convert anistropic uijs to isotropic uijs"""
    uij, single = _reshape_uij(vals=uij)
    iso = numpy.mean(uij[:,0:3],axis=1)
    n = uij.shape[0]
    out1 = numpy.ones((n,3)) * iso.reshape((n,1))
    out2 = numpy.zeros((n,3))
    out = numpy.concatenate([out1, out2], axis=1)
    assert out.shape == uij.shape
    return _revert_output(vals=out, single=single)

def calculate_uij_anisotropy_ratio(uij, tolerance=1e-6):
    """Calculate the ratio of the maximal and the minimal eigenvalues of uij"""
    uij, single = _reshape_uij(vals=uij)
    eigenvalues = numpy.apply_along_axis(uij_eigenvalues, 1, uij)
    maxe = numpy.max(eigenvalues, axis=1)
    mine = numpy.min(eigenvalues, axis=1)
    # Find the atoms with zero eigenvalues
    zero_size = (maxe < tolerance)
    # Set min and max eigenvalues to 1.0 so that ratio gives 1.0
    mine[zero_size] = 1.0
    maxe[zero_size] = 1.0
    # Calculate ratio of minimum to maximum eigenvalues
    out = mine / maxe
    return _revert_output(vals=out, single=single)

def scale_uij_to_target_by_selection(hierarchy, selections, target=1.0, tolerance=1e-12):
    """Change the scale of groups of uij so that the maximum eigenvalue of each selection is `target` (normally for visualisation only). also scale B by the same value."""
    hierarchy = hierarchy.deep_copy()
    cache = hierarchy.atom_selection_cache()
    all_uij = hierarchy.atoms().extract_uij()
    # Pull out copies of the uijs and b for scaling
    out_u = hierarchy.atoms().extract_uij()
    out_b = hierarchy.atoms().extract_b()
    # Iterate through the selections
    for sel_str in selections:
        # Convert to boolsean selection
        sel = cache.selection(sel_str)
        # Extract uijs to calculate scale factor
        uijs = numpy.array(all_uij.select(sel))
        assert uijs.shape[1] == 6
        # Extract the maximum axis length of each uij
        eigs = numpy.apply_along_axis(uij_eigenvalues, axis=1, arr=uijs)
        maxs = numpy.max(eigs, axis=1)
        # Calculate average of the maxima
        mean_max = numpy.mean(maxs)
        # Scale to zero if max eigenvalue is approximately zero
        if mean_max < tolerance:
            mult = 0.0
        else:
            # Calculate the scaling that modifies the mean to correspond to target
            mult = float(target / mean_max)
        # Apply scaling and set szz value
        out_u.set_selected(sel, flex.sym_mat3_double(out_u.select(sel).as_double()*mult))
        out_b.set_selected(sel, out_b.select(sel)*mult)
    # Apply the scaled values to the hierarchy
    hierarchy.atoms().set_uij(out_u)
    hierarchy.atoms().set_b(out_b)
    return hierarchy

