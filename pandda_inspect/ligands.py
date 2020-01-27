import os, shutil, string
import gtk

from bamboo.ccp4_utils import generate_ligand_with_acedrg

def modal_msg(msg):
    """Display an error window - model"""
    d = gtk.MessageDialog(  type    = gtk.MESSAGE_INFO,
                            buttons = gtk.BUTTONS_CLOSE,
                            message_format = msg )
    d.set_position(gtk.WIN_POS_CENTER)
    d.set_keep_above(True)
    d.run()
    d.destroy()

def post_new_ligand_window(output_directory='.'):
    """Display an error window - model"""

    if not os.path.exists(output_directory): os.mkdir(output_directory)

    dialog = gtk.Dialog("Create New Ligand",
                        None,
                        gtk.DIALOG_MODAL,
                        (gtk.STOCK_CANCEL, gtk.RESPONSE_DELETE_EVENT, gtk.STOCK_OK, gtk.RESPONSE_ACCEPT))

    # Name of the ligand
    name_hbox = gtk.HBox(homogeneous=False, spacing=5)
    dialog.vbox.pack_start(name_hbox)
    label = gtk.Label('Ligand Name:')
    label.props.width_chars = 20
    label.set_alignment(0.5, 0.0)
    name_hbox.pack_start(label, expand=True)
    name_entry = gtk.Entry(max=100)
    name_entry.set_text('new-compound')
    name_hbox.pack_start(name_entry)

    # ID for the ligand
    id_hbox = gtk.HBox(homogeneous=False, spacing=5)
    dialog.vbox.pack_start(id_hbox)
    label = gtk.Label('3-letter code:')
    label.props.width_chars = 20
    label.set_alignment(0.5, 0.0)
    id_hbox.pack_start(label)
    id_entry = gtk.Entry(max=3)
    id_entry.set_text('UNL')
    id_hbox.pack_start(id_entry)

    # SMILES for the ligand
    smiles_hbox = gtk.HBox(homogeneous=False, spacing=5)
    dialog.vbox.pack_start(smiles_hbox)
    label = gtk.Label('Smiles:')
    label.props.width_chars = 20
    label.set_alignment(0.5, 0.0)
    smiles_hbox.pack_start(label)
    smiles_entry = gtk.Entry(max=300)
    smiles_entry.set_text('')
    smiles_hbox.pack_start(smiles_entry)

    dialog.show_all()

    disallowed = set(string.punctuation + ' ')

    success = False
    while success is False:
        ligand_pdb = ligand_cif = None

        response = dialog.run()
        # Delete window/cancel?
        if response in [int(gtk.RESPONSE_REJECT), int(gtk.RESPONSE_DELETE_EVENT)]:
            print 'close "new-ligand" window'
            dialog.destroy()
            return None

        assert response is int(gtk.RESPONSE_ACCEPT), 'invalid response received ({} should be {})'.format(response, int(gtk.RESPONSE_ACCEPT))

        ligand_name     = name_entry.get_text().strip(' ')
        ligand_id       = id_entry.get_text().strip(' ')
        ligand_smiles   = smiles_entry.get_text().strip(' ')

        if disallowed.difference(['-','_']).intersection(ligand_name):
            modal_msg('ligand name cannot contain space or punctuation except for {}'.format(' or '.join(['-','_'])))
            continue
        if disallowed.difference(['_']).intersection(ligand_id):
            modal_msg('ligand name cannot contain spaces or punctuation except for {}'.format(' or '.join(['_'])))
            continue

        if len(ligand_name) == 0:
            modal_msg('No ligand name provided')
            continue
        if len(ligand_id) == 0:
            modal_msg('No ligand id provided')
            continue
        if len(ligand_smiles) == 0:
            modal_msg('No ligand smiles provided')
            continue

        ligand_prefix = os.path.join(output_directory, ligand_name)

        print 'ID: {}\nNAME: {}\nSMILES: {}'.format(ligand_id, ligand_name, ligand_smiles)

        ligand_pdb = ligand_prefix+'.pdb'
        ligand_cif = ligand_prefix+'.cif'
        ligand_dir = ligand_prefix+'_TMP'

        if os.path.exists(ligand_prefix+'.pdb') or os.path.exists(ligand_prefix+'.cif'):
            modal_msg('PDB/CIF files already exists called {}.\nPlease select another name.'.format(ligand_name))
            continue

        try:
            ligand_pdb, ligand_cif = generate_ligand_with_acedrg(smiles = ligand_smiles,
                                                                 name   = ligand_id,
                                                                 prefix = ligand_prefix,
                                                                 verbose = True)
        except Exception as e:
            msg = 'Error running acedrg'
            msg += '\n--------------->\n'
            msg += e.message
            if hasattr(e, 'command'):
                msg += '\n--------------->\n'
                msg += e.command.output
            for f in [ligand_pdb,ligand_cif,ligand_dir]:
                if os.path.isfile(f):
                    os.remove(f)
            if os.path.exists(ligand_dir):
                shutil.rmtree(ligand_dir)
            modal_msg(msg)
            continue

        for f in ['_approx.list']:
            f_full = os.path.join(output_directory, f)
            if os.path.exists(f_full):
                os.remove(f_full)

        if os.path.exists(ligand_pdb) and os.path.exists(ligand_cif):
            break

    dialog.destroy()

    return ligand_id, ligand_name, ligand_smiles, ligand_pdb, ligand_cif

