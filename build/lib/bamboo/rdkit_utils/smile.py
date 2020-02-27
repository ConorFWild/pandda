import os, sys, tempfile

from rdkit.Chem import AllChem as Chem

def is_valid_smiles(string):
    """Checks if a string is a valid smile string (not smarts)"""

    string = string.strip()
    # Check for invalid characters
    invalid = [',','~',' ']
    if [1 for inv in invalid if inv in string]:
        return False
    # Check brackets
    brackets = [('(',')'),('[', ']')]
    if [1 for bra in brackets if string.count(bra[0]) != string.count(bra[1])]:
        return False
    else:
        return True

def check_smile_readable(smile):
    """Checks a smile string is interpretable by Rdkit. Returns a mol object if successful"""

    mol = Chem.MolFromSmiles(smile)
    if not mol: raise SmileCheckError('Smile not rdkit-readable: {!s}'.format(smile))
    else:       return mol

def get_canonical_smiles(smile):
    """Makes a Smile String Canonical with the RDKIT framework"""

    if is_valid_smiles(smile):
        mol = Chem.MolFromSmiles(smile)
        if mol:
            return Chem.MolToSmiles(mol, isomericSmiles=True)
        else:
            return None
    else:
        mol = Chem.MolFromSmarts(smile)
        if mol:
            return Chem.MolToSmarts(mol, isomericSmiles=True)
        else:
            return None

def create_pdb_from_smile(smile, filename):
    """Given a Smile String, creates the PDB File"""

    # NEED TO EMBED MOLECULE!!!

    # Create Mol
    mol = Chem.MolFromSmiles(smile)
    # Embed Mol
    Chem.EmbedMolecule(mol)
    # Save File
    out = Chem.MolToPDBFile(mol, filename)

    return out

def get_smiles_from_file_list(filelist):
    """Given a list of PDB Files, find the Smile Strings for the Compounds in the PDB Files"""

    smileslist = []
    # Iterate through the files
    for filename in filelist:
        # Read Smiles
        smile = get_smile_from_file(filename)
        # Append to list
        if smile:
            smileslist.append(smile)

    return smileslist

def get_smile_from_file(file):
    """Given a PDB File containing ONLY a compound, returns the Smile String of the Compound"""
    block = open(file,'r').read()
    return get_smile_from_block(block)

def get_smile_from_block(block):
    """Given a PDB Block containing ONLY a compound, returns the Smile String of the Compound"""

    try:
        # Try to get Mol from MolBlock
        mol = Chem.MolFromPDBBlock(block, removeHs=False, sanitize=True)
    except Exception as Err:
        print Err
        return None
    # Check for existence of Mol
    if not mol: return None
    # Get Smile
    smile = Chem.MolToSmiles(mol, isomericSmiles=False)

    return smile

def match_smile_to_list(search, list, assignbonds=False):
    """Find 'smile' in the list of smile strings"""

    # Filter out those that are none
    filtlist = [(smi,smi) for smi in list if (smi is not None)]
    # Handle the ones with dots
    for full in [smi for smi in list if (smi is not None) and ('.' in smi)]:
        filtlist.append((full,full))
        for part in full.split('.'):
            filtlist.append((full, part))

    for filt in ['C','N','O','F']:
        # Count the number of occurences of 'filt' in smile
        count = search.upper().count(filt)
        # Get those that match
        newfiltlist = [tup for tup in filtlist if tup[1].upper().count(filt) == count]
        # update
        filtlist = newfiltlist

    # Check cyclicity of compound (filter on whether or not cyclicity is present)
    if '1' in search:
        newfiltlist = [tup for tup in filtlist if '1' in tup[1]]
        filtlist = newfiltlist
    else:
        newfiltlist = [tup for tup in filtlist if '1' not in tup[1]]
        filtlist = newfiltlist

    if assignbonds:
        # Assign bonds from the search list
        try:
            bonds = [(tup[0],assign_bond_orders_from_template_smiles(refsmile=tup[1],rawsmile=search)) for tup in filtlist]
        except:
            return [fil[0] for fil in filtlist if fil[1]]
        # Filter out no template matches
        filtlist = [pair for pair in bonds if pair[1]]

    return [fil[0] for fil in filtlist if fil[0]]

def find_structure_matches(search, ref):
    """Finds full or part matches between 'search' smile and 'ref' smile"""

    fullmatch = None
    partmatch = None
    # Assign bonds to the search mol (from PDB file)
    search_bond = assign_bond_orders_from_template_smiles(ref, search)
    # No matches - return None
    if not search_bond:
        return fullmatch, partmatch
    # Canonicalise ref
    ref = get_canonical_smiles(ref)
    # Check for type of match
    if search_bond == ref:
        fullmatch = search_bond
    elif check_for_substructure_smile_matches(ref,search_bond):
        partmatch = search_bond

    return fullmatch, partmatch

def check_for_substructure_smile_matches(refsmile,partsmile):
    """Checks to see if CheckSmile is contained within RefSmile"""

    # Generate Full Molecule from Smile
    refmol = Chem.MolFromSmiles(refsmile)
    # Generate Part Molecule from Smart
    partmol = Chem.MolFromSmarts(partsmile)

    # Check for substructure matches
    return refmol.HasSubstructMatch(partmol)

def assign_bond_orders_from_template_smiles(refsmile, rawsmile, return_mol=False):
    """Takes a template smile (with aromaticity etc.) and applies bond orders to other smile"""

    # Create rawmol
    rawmol = check_smile_readable(rawsmile)
    # Create refmol
    if is_valid_smiles(refsmile):
        refmol = Chem.MolFromSmiles(refsmile)
    else:
        refmol = Chem.MolFromSmarts(refsmile)
    # Check for existence
    if not refmol:
        return None
    # Remove Hydrogens from reference molecule
    try:
        Chem.RemoveHs(refmol, implicitOnly=True)
    except ValueError:
        return None
    # Generated Molecule with transferred bond orders
    try:
        newmol = Chem.AssignBondOrdersFromTemplate(refmol,rawmol)
    except ValueError:
        return None

    if return_mol:
        return newmol
    else:
        # Convert Back To Smiles
        newsmiles = Chem.MolToSmiles(newmol, isomericSmiles=True)
        return newsmiles

