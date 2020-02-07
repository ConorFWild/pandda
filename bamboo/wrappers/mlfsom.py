import os, shutil, tempfile, tarfile

import bamboo.resources

from bamboo.common.command import CommandManager

def run_xia2(img_dir, out_dir, method='dials', verbose=True):
    """Run xia2 on the images in img_dir, and put the results in out_dir"""

    assert method in ['2d','3d','3dii','dials'], 'Invalid method provided'

    assert os.path.exists(img_dir), 'Image directory does not exist'
    if not os.path.exists(out_dir): os.mkdir(out_dir)

    # Convert to absolute path to allow folder changing
    img_dir = os.path.abspath(img_dir)
    out_dir = os.path.abspath(out_dir)

    # Store current directory
    cur_dir = os.getcwd()

    try:

        # Output logfile
        xia2_log = os.path.join(out_dir, 'xia2.log')

        # Got to directory to run MLFSOM
        os.chdir(out_dir)

        xia2 = CommandManager('xia2')
        xia2.add_command_line_arguments('-'+method, img_dir)
        if verbose: print('Running xia2 on {}').format(img_dir)
        rc = xia2.run()
        if rc != 0:
            print('STDOUT:')
            print(xia2.output)
            print('STDERR:')
            print(xia2.error)
        with open(xia2_log, 'a') as fh:
            fh.write('==================================================>\n')
            fh.write('xia2 log\n')
            fh.write('==================================================>\n')
            fh.write(xia2.output+'\n')
            fh.write('==================================================>\n')
            fh.write('xia2 errors\n')
            fh.write('==================================================>\n')
            fh.write(xia2.error+'\n')
            fh.write('==================================================>\n')

    finally:
        # Change back to original directory
        os.chdir(cur_dir)

    return xia2

def run_mlfsom(pdb_file, res_h, out_dir='mlfsom_out', n_images=180, osc=1.0, energy=12660, verbose=True):
    """Run MLFSOM on the given file."""

    ftd_pre = ['1H87.pdb', 'example.com', 'example.mat', 'mlfsom.log', 'refined.pdb', 'pristine.mtz']

    assert not os.path.exists(out_dir), 'output directory already exists!'
    assert os.path.exists(pdb_file)
    assert isinstance(n_images, int)

    # Store current directory
    cur_dir = os.getcwd()

    try:
        # Create a temporary directory for mlfsom to run in (and go to it)
        tmp_dir = tempfile.mkdtemp(prefix='mlfsom_')
        mlfsom_pdb = os.path.join(tmp_dir, 'refined.pdb')
        mlfsom_mtz = os.path.join(tmp_dir, 'pristine.mtz')
        mlfsom_log = os.path.join(tmp_dir, 'mlfsom.log')
        mlfsom_img = 'image_###.img'

        # Get the mlfsom tarball and move it to the temporary directory
        mlfsom_file = os.path.join(bamboo.resources.__path__[0], 'mlfsom.tar.gz')
        assert os.path.exists(mlfsom_file), 'MLFSOM tarball not found'
        shutil.copy(mlfsom_file, tmp_dir)
        mlfsom_file = os.path.join(tmp_dir, 'mlfsom.tar.gz')
        assert os.path.exists(mlfsom_file), 'MLFSOM tarball not in temporary directory'

        # Untar the file
        t = tarfile.open(mlfsom_file, 'r:gz')
        try: t.extractall(tmp_dir)
        finally: t.close()

        # Delete unnecessary files
        for f in ftd_pre:
            fp = os.path.join(tmp_dir, f)
            if os.path.exists(fp): os.remove(fp)

        # Move the input files in
        shutil.copy(pdb_file, mlfsom_pdb)

        # Got to directory to run MLFSOM
        os.chdir(tmp_dir)

        # Convert the pdb_file to structure factors
        ano_sfall = CommandManager('./ano_sfall.com')
        ano_sfall.add_command_line_arguments(os.path.basename(mlfsom_pdb))
        ano_sfall.add_command_line_arguments('energy={}'.format(energy))
        ano_sfall.add_command_line_arguments('{}A'.format(res_h))
        if verbose: print('Running ./ano_sfall.com on {}').format(pdb_file)
        rc = ano_sfall.run()
        if rc != 0:
            print('STDOUT:')
            print(ano_sfall.output)
            print('STDERR:')
            print(ano_sfall.error)
        shutil.move('ideal_ano.mtz', mlfsom_mtz)
        with open(mlfsom_log, 'a') as fh:
            fh.write('==================================================>\n')
            fh.write('ano_sfall log\n')
            fh.write('==================================================>\n')
            fh.write(ano_sfall.output+'\n')
            fh.write('==================================================>\n')
            fh.write('ano_sfall errors\n')
            fh.write('==================================================>\n')
            fh.write(ano_sfall.error+'\n')
            fh.write('==================================================>\n')

        # Run MLFSOM on the resulting mtzfile
        mlfsom = CommandManager('./mlfsom.com')
        mlfsom.add_command_line_arguments(mlfsom_img)
        mlfsom.add_command_line_arguments('frames={}'.format(n_images))
        mlfsom.add_command_line_arguments('osc={}'.format(osc))
        if verbose: print('Running ./mlfsom.com on {}. This may take a while...').format(pdb_file)
        rc = mlfsom.run()
        if rc != 0:
            print('STDOUT:')
            print(mlfsom.output)
            print('STDERR:')
            print(mlfsom.error)
        with open(mlfsom_log, 'a') as fh:
            fh.write('==================================================>\n')
            fh.write('mlfsom log\n')
            fh.write('==================================================>\n')
            fh.write(mlfsom.output+'\n')
            fh.write('==================================================>\n')
            fh.write('mlfsom errors\n')
            fh.write('==================================================>\n')
            fh.write(mlfsom.error+'\n')
            fh.write('==================================================>\n')

    finally:
        # Change back to original directory
        os.chdir(cur_dir)
        # Delete temporary directory
        shutil.move(tmp_dir, out_dir)

    return ano_sfall, mlfsom
