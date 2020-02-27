
import numpy, pandas
import iotbx.pdb

from scitbx.array_family import flex
from giant.maths.geometry import is_within

import iotbx.pdb.hierarchy as iotbx_pdbh

###############################################################################
###                                                                         ###
###                      General Hierarchy Functions                        ###
###                                                                         ###
###############################################################################

def sanitise_hierarchy(hierarchy, in_place=True):
    """Remove several common structure problems"""

    if not in_place: hierarchy = hierarchy.deep_copy()

    # Fix element problems
    hierarchy.atoms().set_chemical_element_simple_if_necessary()
    # Sort atoms (allows easier comparison)
    hierarchy.sort_atoms_in_place()

    return hierarchy

###############################################################################
###                                                                         ###
###                                LABELS                                   ###
###                                                                         ###
###############################################################################

def make_label(obj):
    """Return the relevant label for a supplied hierarchy/atom object"""
    if isinstance(obj, iotbx_pdbh.residue_group):
        return label_from_residue_group(obj)
    elif isinstance(obj, iotbx_pdbh.atom_group):
        return label_from_atom_group(obj)
    elif isinstance(obj, iotbx_pdbh.conformer):
        return label_from_conformer(obj)
    elif isinstance(obj, iotbx_pdbh.residue):
        return label_from_residue(obj)
    elif isinstance(obj, iotbx_pdbh.atom):
        return label_from_atom(obj)
    elif isinstance(obj, iotbx_pdbh.atom_with_labels):
        return label_from_atom_with_labels(obj)
    else:
        raise Exception('Invalid object type provided: {}'.format(type(obj)))

def label_from_residue_group(rg):
    """Return (chain_id, resid)"""
    ch = rg.parent()
    return (ch.id, rg.resid())

def label_from_atom_group(ag):
    """Return (chain_id, resid, altloc)"""
    rg = ag.parent()
    ch = rg.parent()
    return (ch.id, rg.resid(), ag.altloc)

def label_from_conformer(cnf):
    """Return (chain_id, resid, altloc). Must only have one residue."""
    ch = cnf.parent()
    res = cnf.only_residue()
    return (ch.id, res.resid(), cnf.altloc)

def label_from_residue(res):
    """Return (chain_id, resid, altloc)."""
    cnf = res.parent()
    ch = cnf.parent()
    return (ch.id, res.resid(), cnf.altloc)

def label_from_atom(at):
    """Return (chain_id, resid, altloc)"""
    ag = at.parent()
    rg = ag.parent()
    ch = rg.parent()
    return (ch.id, rg.resid(), ag.altloc)

def label_from_atom_with_labels(at):
    """Return (chain_id, resid, altloc)"""
    return (at.chain_id, at.resid(), at.altloc)

###############################################################################
###                                                                         ###
###                         CONVENIENCE FUNCTIONS                           ###
###                                                                         ###
###############################################################################

def get_atom_pairs(residue_1, residue_2, fetch_labels=True):
    atom_pairs = []
    a1 = residue_1.atoms()
    a2 = residue_2.atoms()
    a1_d = a1.build_dict()
    a2_d = a2.build_dict()
    for a_name in a1.extract_name():
        if a_name not in a2_d: continue
        if fetch_labels:
            atom_pairs.append((a1_d[a_name].fetch_labels(),a2_d[a_name].fetch_labels()))
        else:
            atom_pairs.append((a1_d[a_name],a2_d[a_name]))
    return atom_pairs

def calculate_residue_group_occupancy(residue_group):
    """
    Extract the total occupancy of a residue, allowing for alternate conformers.
    - If conformers are present, the occupancies of each conformer are summed.
    - If multiple occupancies are present for a conformer, the maximum is taken for each.
    """
    rg = residue_group
    if [c.altloc for c in rg.conformers()] == ['']:
        res_occ = max(rg.atoms().extract_occ())
    else:
        res_occ = sum([max(c.atoms().extract_occ()) for c in rg.conformers() if c.altloc])
    return res_occ

def calculate_paired_conformer_rmsds(conformers_1, conformers_2):
    """Return a list of rmsds between two list of conformers, paired by the minimum rmsd. Each conformer may only be paired once"""

    rmsds = numpy.empty((len(conformers_1), len(conformers_2)))
    rmsds[:] = None

    for i, c1 in enumerate(conformers_1):
        for j, c2 in enumerate(conformers_2):
            rmsds[i,j] = calculate_paired_atom_rmsd(atoms_1=c1.atoms(), atoms_2=c2.atoms(), sort=True, truncate_to_common_set=True, remove_H=True)
    ret_list = []
    print rmsds
    while not numpy.isnan(rmsds).all():
        min_val = numpy.nanmin(rmsds)
        i1,i2 = zip(*numpy.where(rmsds==min_val))[0]
        # Clear these values so that a conformer cannot be used more than once
        rmsds[i1, :] = None
        rmsds[:, i2] = None
        # Return conformers and the rmsd
        ret_list.append((conformers_1[i1], conformers_2[i2], min_val))
        print rmsds
        print ret_list[-1]
    return ret_list

def calculate_paired_atom_rmsd(atoms_1, atoms_2, sort=True, truncate_to_common_set=True, remove_H=True):
    """
    Calculate the RMSD between two sets of equivalent atoms (i.e. atoms with the same names)
    - sort                   - sort atoms prior to rmsd calculation (so that equivalent atoms are compared)
    - truncate_to_common_set - only calculate rmsd over the common set of atoms, return None if different atoms are provided.
    Raises error if atom names are present more than once in each list (e.g. alternate conformers of atoms)
    """

    # Need to sort by atoms to check if atoms are the same in both sets
    if truncate_to_common_set:
        sort=True

    # Remove hydrogen atoms
    if remove_H:
        atoms_1 = atoms_1.select(atoms_1.extract_element() != ' H')
        atoms_2 = atoms_2.select(atoms_2.extract_element() != ' H')

    # Sort atoms by name if required
    if sort:
        # Extract atom names
        names_1 = list(atoms_1.extract_name())
        names_2 = list(atoms_2.extract_name())
        # Get name overlap between atom sets
        common_atoms = list(set(names_1).intersection(names_2))
        if (not truncate_to_common_set):
            if not (len(common_atoms)==len(atoms_1)==len(atoms_2)):
                return None
        # If sorting, require atom names unique
        assert len(set(names_1)) == len(names_1), 'atom names can only be present once in each atom list -- no alternate conformers of atoms'
        assert len(set(names_2)) == len(names_2), 'atom names can only be present once in each atom list -- no alternate conformers of atoms'
        # Get reorderings of atoms
        sort_1 = flex.size_t([names_1.index(an) for an in common_atoms])
        sort_2 = flex.size_t([names_2.index(an) for an in common_atoms])
        # Reorder the atoms
        atoms_1 = atoms_1.select(sort_1)
        atoms_2 = atoms_2.select(sort_2)

    # Check same number of atoms and atom names are the same, etc...
    try:
        assert atoms_1.size() == atoms_2.size()
        assert atoms_1.extract_name().all_eq(atoms_2.extract_name())
    except:
        return None
    # Calculate RMSD and return
    return flex.mean((atoms_1.extract_xyz() - atoms_2.extract_xyz()).dot())**0.5

def normalise_occupancies(atoms, max_occ=1.0):
    """Normalise the maximum occupancy of a group of atoms to max_occ"""
    occ = atoms.extract_occ()
    occ_mult = max_occ/max(occ)
    atoms.set_occ(occ*occ_mult)
    return atoms

