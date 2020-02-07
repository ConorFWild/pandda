import os, sys, math

from rdkit import Chem
from rdkit.Chem import BRICS
from bamboo.macro.molecule import MacroMol
from bamboo.rdkit_utils.mol import check_pdb_readable

def break_on_rotatable_bonds_to_mol(inmol):
    """Takes a mol and breaks it on all of the rotatable bonds (well, most of them) - returns a single mol, or None if there are no rotatable bonds"""

    # Get the indices of the atoms around the bonds to be broken
    atom_pairs = [p[0] for p in BRICS.FindBRICSBonds(inmol) if p]
    # Return if no bonds found... as it was given
    if not atom_pairs:
        return inmol
    # Get the bond indices
    bonds = [inmol.GetBondBetweenAtoms(at1,at2).GetIdx() for (at1, at2) in atom_pairs]
    # Fragment the molecule
    fragged_mol = Chem.rdmolops.FragmentOnBonds(inmol, bonds, addDummies=False)

    return fragged_mol

def break_on_rotatable_bonds_to_mols(inmol):
    """Takes a mol and breaks it on all of the rotatable bonds (well, most of them) - returns a set of mols: one for each fragment, or [] if there are no rotatable bonds"""

    fragged_mol = break_on_rotatable_bonds_to_mol(inmol)
    # No fragged mol => No bonds to fragment => return original mol
    if not fragged_mol:
        raise RDkitFragmentationError('Failed to break on rotatable bonds')
    # Get the individual mols for the fragments
    frags = Chem.GetMolFrags(fragged_mol,asMols=True)

    return frags

def break_and_rename_mol_to_file(infile, outfile):
    """Takes a ligand, breaks it up and renames the different fragments with inscodes. compiles into one file"""

    inmol = check_pdb_readable(infile)
    # PDB-rdkit disagreement...
    if not inmol:
        raise RDkitReadError('Cannot read {!s}'.format(infile))

    orig_mol = MacroMol(infile)

    assert len(orig_mol.getChains())==1, 'MORE THAN ONE CHAIN PRESENT IN FILE! {!s}'.format(infile)
    assert len(orig_mol.getResidues())==1, 'MORE THAN ONE RESIDUE PRESENT IN FILE! {!s}'.format(infile)

    # Keep the atom numbers (these are conserved in the mappings)
    atom_numbers = dict([(a.atomname,a.serialnum) for a in orig_mol.getResidues()[0].atoms])

    # Break the molecule into fragments (or get inmol back if no rotatable bonds)
    frags = break_on_rotatable_bonds_to_mols(inmol)

    # List of output pdb blocks
    out_blocks = []

    for i, frag in enumerate(frags):

        raw_block = Chem.MolToPDBBlock(frag)
        raw_frag = MacroMol(raw_block)
        # Check there's only one residue in file
        assert len(raw_frag.getChains())==1, 'FRAGMENT HAS MORE THAN ONE CHAIN? {!s}\n{!s}'.format(infile, raw_block)
        assert len(raw_frag.getResidues())==1, 'FRAGMENT HAS MORE THAN ONE RESIDUE? {!s}\n{!s}'.format(infile, raw_block)
        # Form the new inscode
        new_inscode = chr(i + ord('A'))
        # Change the inscode
        residue = raw_frag.getResidues()[0]
        residue.update_inscode(new_inscode)
        # Change the atom numbers back to those in the original file for continuity
        residue.update_atom_numbers(atom_numbers)
        # Get the pdb string and use
        proc_block = residue.getPdbString()
        # Check for consistency and append
        out_blocks.append(proc_block)

    # Add easy-header
    header_block = ''.join(orig_mol.easy_header).strip('\n')
    # Add easy-footer
    footer_block = ''.join(orig_mol.easy_footer).strip('\n')
    # Create atom block
    out_block = ''.join(out_blocks).replace('END','').strip('\n')
    # Join the blocks together
    block_to_write = '\n'.join([header_block, out_block, footer_block, 'END\n'])
    # Check for consistency
    block_to_write = block_to_write.replace('\n\n','\n')

    open(outfile,'w').write(block_to_write)

    return outfile

