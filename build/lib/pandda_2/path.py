import os, sys, shutil, glob, tempfile

#:::::::::::::::::::::::::::::::::#
# ############################### #
# ### Path Utility Functions  ### #
# ############################### #
#:::::::::::::::::::::::::::::::::#

def foldername(f):
    return os.path.basename(os.path.dirname(f))

def filename(f):
    return os.path.basename(os.path.splitext(f)[0])

def easy_directory(directory):
    """Checks a directory exists and creates it if not"""
    if not os.path.exists(directory):
        os.mkdir(directory)
    return directory

def rel_symlink(orig, link):
    """Make a rela  tive symlink from link to orig"""
    assert os.path.exists(orig), 'FILE DOES NOT EXIST: {!s}'.format(orig)
    assert not os.path.exists(link), 'LINK ALREADY EXISTS: {!s}'.format(link)
    orig = os.path.abspath(orig)
    link = os.path.abspath(link)
    assert not link.endswith('/'), 'LINK CANNOT END WITH /'
    os.symlink(os.path.relpath(orig, start=os.path.dirname(link)), link)

def splice_ext(path, new, position=-1):
    dirname, basename = os.path.split(path)
    split = basename.split(os.path.extsep)
    assert abs(position) <= len(split), 'invalid position selected ({}) to insert "{}" into {}'.format(position, new, split)
    spliced = split[:position] + [new] + split[position:]
    joined = os.path.extsep.join(spliced)
    return os.path.join(dirname, joined)

def delete_with_glob(glob_str, verbose=True):
    for file in sorted(glob.glob(glob_str)):
        # I know that unlink is the same as remove...
        if os.path.islink(file):
            if verbose: print('Deleting link: {}'.format(file))
            os.unlink(file)
        elif os.path.isdir(file):
            pass
        else:
            if verbose: print('Deleting file: {}'.format(file))
            os.remove(file)

def delete_temporary_directory(tempdir):
    """Delete A Temporary Directory"""

    # Remove Temporary Files
    if tempdir.startswith(tempfile.gettempdir()):
        shutil.rmtree(tempdir)
        return 0
    else:
        FlagError("Will Not Delete This Directory - Not in '{}' - Remove Manually : {}".format(tempfile.gettempdir(), tempdir))
        return 1

def list_directory_contents(directory, templatename=None, templatestyle=None):
    """List Folders in a Directory"""

    # List the contents of the directory
    contents = os.listdir(directory)
    # Decide which are directories and which to keep
    if   templatename and templatestyle=='End':
        subdirs = [dir for dir in contents if os.isdir(os.path.join(directory,dir)) and dir.endswith(templatename)]
    elif templatename and templatestyle=='Start':
        subdirs = [dir for dir in contents if os.isdir(os.path.join(directory,dir)) and dir.startswith(templatename)]
    elif templatename and templatestyle=='Contains':
        subdirs = [dir for dir in contents if os.isdir(os.path.join(directory,dir)) and templatename in dir]
    else:
        subdirs = [dir for dir in contents if os.isdir(os.path.join(directory,dir))]
    # Construct Full Paths
    subdirs = [os.path.join(directory,dir) for dir in subdirs]

    return subdirs

