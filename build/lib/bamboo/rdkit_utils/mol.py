import os,sys,re

from bamboo.common.errors import *
from bamboo.rdkit_utils.smile import check_smile_readable

from rdkit.Chem import AllChem as Chem

return_codes = {0 : 'Reference {!s} matches the query {!s} as far as I can tell!',
                1 : 'Reference {!s} has no matching substructure with query {!s}',
                2 : 'Reference {!s} has more atoms than query {!s}',
                3 : 'Reference {!s} has fewer atoms than query {!s}',
                4 : 'Reference {!s} has more bonds than query {!s}',
                5 : 'Reference {!s} has fewer bonds than query {!s}',
                6 : 'Reference {!s} has a different charge to query {!s}',
                7 : 'Reference {!s} does not have a full substructure match to query {!s}'}

def check_pdb_readable(file):
    """Checks a file is readable by Rdkit. Returns a mol object if successful"""

    mol = Chem.MolFromPDBFile(file)
    if not mol: raise LigandCheckError('File not rdkit-readable: {!s}'.format(file))
    else:       return mol

def has_bond_orders(mol):
    """Checks to see if the given mol has bond orders assigned"""

    if len([1 for b in mol.GetBonds() if b.GetBondTypeAsDouble() != 1.0]):
        return True
    else:
        return False

def check_and_assign_bond_orders_dynamic(mol1, mol2):
    """Checks both mol1 and mol2 to see if they have bond orders. If one does, but the other does not, then it attempts to transfer bond orders appropriately"""

    try:
        if has_bond_orders(mol1) and has_bond_orders(mol2):
            # Both have bond orders - do nothing
            pass
        elif has_bond_orders(mol1) and (not has_bond_orders(mol2)):
            # mol1 has, but mol2 does not - transfer from mol1 -> mol2
            mol2 = Chem.AssignBondOrdersFromTemplate(mol1, mol2)
        elif (not has_bond_orders(mol1)) and has_bond_orders(mol2):
            # mol2 has, but mol1 does not - transfer from mol2 -> mol1
            mol1 = Chem.AssignBondOrdersFromTemplate(mol2, mol1)
        else:
            # Neither has bond orders - pass
            pass
    except ValueError as err:
        raise BondAssignmentError('Failed to assign bond orders: {!s}'.format(err))

    if not mol1:
        raise Exception('Mol1 does not exist!')
    elif not mol2:
        raise Exception('Mol2 does not exist!')

    return mol1, mol2

def check_and_assign_bond_orders_static(mol1, mol2):
    """Transfers bond orders from mol1 to mol2, if bond orders are present"""

    try:
        # transfer from mol1 -> mol2
        mol2 = Chem.AssignBondOrdersFromTemplate(mol1, mol2)
    except ValueError as err:
        raise BondAssignmentError('Failed to assign bond orders: {!s}'.format(err))

    if not mol1:
        raise Exception('Mol1 does not exist!')
    elif not mol2:
        raise Exception('Mol2 does not exist!')

    return mol1, mol2

def check_ligands_from_file_identical(ref, query, dynamicBondAssignment=True):
    """Check two ligands (pdb files) are identical - gives similarity relative to the first file"""
    rc = check_ligand_mols_identical(check_pdb_readable(ref), check_pdb_readable(query), dynamicBondAssignment=dynamicBondAssignment)
    return rc, return_codes[rc].format('file','file')

def check_ligand_from_file_identical_to_smile(smile, query, dynamicBondAssignment=False):
    """Checks ligand (pdb file) is the same as the smile"""
    rc = check_ligand_mols_identical(check_smile_readable(smile), check_pdb_readable(query), dynamicBondAssignment=dynamicBondAssignment)
    return rc, return_codes[rc].format('smile','file')

def check_ligand_identical_to_smile(smile, mol, dynamicBondAssignment=False):
    """Checks ligand (mol) is the same as the smile"""
    rc = check_ligand_mols_identical(ref=check_smile_readable(smile), mol=mol, dynamicBondAssignment=dynamicBondAssignment)
    return rc, return_codes[rc].format('smile','mol')

def check_ligand_mols_identical(ref, mol, dynamicBondAssignment=True):
    """Check two ligands (mols) are identical - gives similarity relative to the first mol"""

    # Assign bond orders (and return if fails)
    try:
        if dynamicBondAssignment:
            # Keeps existing bond orders, or transfers if not present
            ref, mol = check_and_assign_bond_orders_dynamic(ref, mol)
        else:
            # Assigns bond orders from ref -> mol
            ref, mol = check_and_assign_bond_orders_static(ref, mol)
    except BondAssignmentError:
        return 1

    # Check Number of Atoms
    if ref.GetNumAtoms() < mol.GetNumAtoms():
        return 2
    if ref.GetNumAtoms() > mol.GetNumAtoms():
        return 3
    # Check Number of Bonds
    if ref.GetNumBonds() < mol.GetNumBonds():
        return 4
    if ref.GetNumBonds() > mol.GetNumBonds():
        return 5
    # Check Charge (This is where a lot of programs go 'wrong')
    if Chem.GetFormalCharge(ref) != Chem.GetFormalCharge(mol):
        return 6
#    try:
    # Check Substructure Matches
    if ref.GetNumAtoms() != len(mol.GetSubstructMatch(ref)):
        return 7
#    except ValueError as err:
#        raise LigandCheckError(' Value Error in GetSubstructMatch: {!s}'.format(err))

    return 0
