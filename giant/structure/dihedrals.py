import numpy

import iotbx.pdb
from scitbx.math import dihedral_angle

from bamboo.maths.angles import angle_difference
from giant.structure.select import extract_backbone_atoms
from giant.structure.iterators import generate_residue_triplets

########################################################################################

def get_all_phi_psi_for_hierarchy(pdb_h, deg=True):
    """Calculate all phi-psi angles for a structure"""
    for prev, current, next in generate_residue_triplets(pdb_h):
        phi, psi = calculate_residue_phi_psi_angles(prev=prev, current=current, next=next, deg=deg)
        yield (current, (phi,psi))

########################################################################################

def calculate_residue_phi_psi_angles(prev, current, next, deg=True):
    """Calculate phi-psi angles for a set of three residues. Returns (phi, psi)."""
    phi = dihedral_angle(sites=extract_phi_sites(prev=prev, current=current), deg=deg)
    psi = dihedral_angle(sites=extract_psi_sites(current=current, next=next), deg=deg)
    return (phi, psi)

########################################################################################

def extract_psi_sites(current, next):
    """Extract sites to calculate psi angle from the selected residues"""
    return [a.xyz for a in extract_psi_atoms(current=current, next=next)]

def extract_phi_sites(prev, current):
    """Extract sites to calculate phi angle from the selected residues"""
    return [a.xyz for a in extract_phi_atoms(prev=prev, current=current)]

########################################################################################

def extract_psi_atoms(current, next):
    """Extract atoms to calculate phi angle from the selected residues"""
    n_curr = ca_curr = c_curr = n_next = None
    n_next, ca_next, c_next = extract_backbone_atoms(residue=next)
    n_curr, ca_curr, c_curr = extract_backbone_atoms(residue=current)
    return (n_curr, ca_curr, c_curr, n_next)

def extract_phi_atoms(prev, current):
    """Extract sites to calculate phi angle from the selected residues"""
    c_prev = n_curr = ca_curr = c_curr = None
    n_prev, ca_prev, c_prev = extract_backbone_atoms(residue=prev)
    n_curr, ca_curr, c_curr = extract_backbone_atoms(residue=current)
    return (c_prev, n_curr, ca_curr, c_curr)

########################################################################################

def get_peptide_bond_type(prev_Ca, prev_C, curr_N, curr_Ca):
    "Return the type of peptide bond"

    assert prev_Ca.name.strip() == 'CA'
    assert prev_C.name.strip() == 'C'
    assert curr_N.name.strip() == 'N'
    assert curr_Ca.name.strip() == 'CA'

    ang = dihedral_angle(sites=[a.xyz for a in (prev_Ca, prev_C, curr_N, curr_Ca)], deg=True)
    ang = numpy.round(ang, 2)

    if curr_N.parent().resname == 'PRO':
        cis, trans = ('PCIS','PTRANS')
    else:
        cis, trans = ('CIS', 'TRANS')

    # Use lenient angle tolerances to detect CIS or TRANS
    if angle_difference(a1=0, a2=ang, deg=True, abs_val=True) < 45:
        bond_type = cis
    elif angle_difference(a1=180, a2=ang, deg=True, abs_val=True) < 45:
        bond_type = trans
    else:
        print 'WARNING: BOND IS NOT CIS OR TRANS (angle {:7}) for link between {} and {} - DEFAULTING TO TRANS'.format(ang, prev_C.id_str(), curr_N.id_str())
        bond_type = trans

    return bond_type
