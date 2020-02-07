import os, sys, glob, shutil
import wget

def get_pdb_index(filename='pdb_entry_type.txt'):
    if os.path.exists(filename): return filename
    tempfile = wget.download('ftp://ftp.wwpdb.org/pub/pdb/derived_data/pdb_entry_type.txt')
    if tempfile != filename: shutil.move(tempfile, filename)
    assert os.path.exists(filename)
    return filename

def download_pdb(pdb_id, filename):
    """Download pdb file to a filename"""
    if os.path.exists(filename): return filename
    try:    tempfile = wget.download('http://www.rcsb.org/pdb/files/{}.pdb'.format(pdb_id))
    except: print('Failed to download:{}'.format(pdb_id)); raise
    if tempfile != filename: shutil.move(tempfile, filename)
    assert os.path.exists(filename)
    return filename

def download_structure_factors(pdb_id, filename, zipped=True):
    """Download pdb file to a filename"""
    if os.path.exists(filename): return filename
    try:
        if zipped:  tempfile = wget.download('https://files.rcsb.org/download/{}-sf.cif.gz'.format(pdb_id))
        else:       tempfile = wget.download('https://files.rcsb.org/download/{}-sf.cif'.format(pdb_id))
    except: print('Failed to download:{}'.format(pdb_id)); raise
    if tempfile != filename: shutil.move(tempfile, filename)
    assert os.path.exists(filename)
    return filename

def download_and_organise_pdb_as_necessary(pdb_id, out_dir='pdb-download'):
    """Download and organise the pdb file for a pdb code"""
    # Check the top directory
    if not os.path.exists(out_dir):
        try:    os.mkdir(out_dir)
        except: pass
    # Create subdirectory for the pdb file
    pdb_dir = os.path.join(out_dir, pdb_id[0])
    if not os.path.exists(pdb_dir):
        try:    os.mkdir(pdb_dir)
        except: pass
    # Get the pdb file
    pdb_file = os.path.join(pdb_dir, pdb_id+'.pdb')
    # Download file if necessary
    if not os.path.exists(pdb_file): download_pdb(pdb_id=pdb_id, filename=pdb_file)
    return pdb_file
