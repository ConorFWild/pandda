import os, sys, copy, re, glob

import libtbx.phil
import libtbx.easy_mp

import numpy, pandas

from bamboo.html import BAMBOO_HTML_ENV
from giant.jiffies.score_model import prepare_output_directory, score_model, make_residue_radar_plot, format_parameters_for_plot, residue_plot_phil

#######################################

bar = '=======================++>'

PROGRAM = 'giant.score_model_multiple'
DESCRIPTION = """
    A tool to quickly score multiple ligand models against crystallographic electron density

    1) Simple usage (for ligand called LIG or UNL, and files in each folder called refine.pdb, refine.mtz)
        > giant.score_model dir='/path/to/directories/*'

    See giant.score_model for more information.
"""

blank_arg_prepend = {None:'dir='}

master_phil = libtbx.phil.parse("""
input {
    dir = None
        .type = path
        .multiple = True

    pdb_1_style = 'refine.pdb'
        .type = str
    mtz_1_style = 'refine.mtz'
        .type = str

    pdb_2_style = '*-ensemble-model.pdb'
        .type = str
    mtz_2_style = '*-input.mtz'
        .type = str

    labels = basename *folder_name
        .type = choice
        .multiple = False
}
selection {
    res_names = LIG,UNL,DRG
        .type = str
        .help = "Comma-separated list of residue names to score -- if None then scores all"
}
output {
    out_dir = ./
        .type = path
    generate_summary = True
        .type = bool
}
settings{
    cpus = 1
        .type = int
}
include scope giant.jiffies.score_model.residue_plot_phil
""", process_includes=True)

#######################################

def process_folder_mp(args):
    try:
        return process_folder(*args)
    except Exception as e:
        return e

def process_folder(dir, params):
    """Score the files in a folder"""

    try:    pdb1 = glob.glob(os.path.join(dir, params.input.pdb_1_style))[0]
    except: raise Exception('FAILED TO FIND FILE: {} in {}'.format(params.input.pdb_1_style, dir))
    try:    mtz1 = glob.glob(os.path.join(dir, params.input.mtz_1_style))[0]
    except: raise Exception('FAILED TO FIND FILE: {} in {}'.format(params.input.mtz_1_style, dir))
    if params.input.pdb_2_style:
        try:    pdb2 = glob.glob(os.path.join(dir, params.input.pdb_2_style))[0]
        except: raise Exception('FAILED TO FIND FILE: {} in {}'.format(params.input.pdb_2_style, dir))
    else: pdb2 = None
    if params.input.mtz_2_style:
        try:    mtz2 = glob.glob(os.path.join(dir, params.input.mtz_2_style))[0]
        except: raise Exception('FAILED TO FIND FILE: {} in {}'.format(params.input.mtz_2_style, dir))
    else: mtz2 = None

    if   params.input.labels == 'basename':     label = os.path.split_ext(os.path.basename(pdb1))[0]
    elif params.input.labels == 'folder_name':  label = os.path.basename(dir)

    print 'Scoring Model: {}'.format(label)

    data_table = score_model(   params = params,
                                pdb1   = pdb1,
                                mtz1   = mtz1,
                                pdb2   = pdb2,
                                mtz2   = mtz2,
                                label_prefix = label
                            )

    return data_table

#######################################

def run(params):

    assert params.input.dir, 'No directories provided'
    assert params.input.pdb_1_style
    assert params.input.mtz_1_style
    # REMOVE THIS WHEN IMPLEMENTED
    assert params.selection.res_names
    # REMOVE THIS WHEN IMPLEMENTED
    if params.selection.res_names:      params.selection.__inject__("res_names_list", params.selection.res_names.split(','))
    else:                               params.selection.__inject__("res_names_list", None)

    output_dir, images_dir = prepare_output_directory(params)
    scores_file  = os.path.join(output_dir, 'residue_scores.csv')
    summary_file = os.path.join(output_dir, 'residue_scores.html')

    arg_list = []
    # Filter directories
    for dir in params.input.dir:
        if not os.path.isdir(dir): continue
        # Build arglist for mapping
        arg_list.append((dir, params))
    # Score files in parallel
    print bar
    print 'Scoring all models...'
    returned_objs = libtbx.easy_mp.pool_map(fixed_func=process_folder_mp, args=arg_list, processes=params.settings.cpus)
    # Extract results
    data_table = None
    for obj in returned_objs:
        if isinstance(obj, pandas.DataFrame):
            if data_table is None:  data_table = obj
            else:                   data_table = data_table.append(obj, verify_integrity=True)
        else:
            try:    print '{}: {}'.format(obj.__class__, obj.message)
            except: print 'Unknown Error: {}'.format(str(obj))
    print '...Done'
    print bar

    data_table.to_csv(scores_file)
    print 'Output written to {}'.format(scores_file)
    print bar

    ###################################################################
    # Image parameters
    ###################################################################
    columns = format_parameters_for_plot(params=params.plot.parameters)

    ###################################################################
    # Output Images
    ###################################################################
    all_images = []
    print 'Generating Output Images...'
    for label, row in data_table.iterrows():
        image_path = os.path.join(images_dir,'{}.png'.format(label))
        print 'Making: {}...'.format(image_path)
        make_residue_radar_plot(path = image_path,
                                data = row.to_frame().T,
                                columns = columns,
                                remove_blank_entries = params.plot.remove_blank_entries,
                                print_axis_values    = params.plot.print_axis_values  )
        all_images.append(image_path)
    print '...Done'
    print bar

    if params.output.generate_summary:
        print 'Generating output HTML...'
        # Get template to be filled in
        template = BAMBOO_HTML_ENV.get_template('summary_page.html')
        # Output directory (for relative symlinks)
        out_dir  = os.path.abspath(os.path.dirname(summary_file))

        # ===========================================================>
        # Construct the data object to populate the template
        output_data = {}
        output_data['header'] = 'Residue Score Summaries'
        output_data['title'] = 'Residue Score Summaries'
        output_data['introduction'] = 'Model Quality and Validation checks.'
        # ===========================================================>
        # Header Images
        output_data['small_images'] = []
        for img in all_images:
            output_data['small_images'].append({ 'path': './'+os.path.relpath(path=img, start=out_dir),
                                                 'title': 'Scores for {}'.format(os.path.splitext(os.path.basename(img))[0]) })
        # ===========================================================>
        # Write Output
        with open(summary_file, 'w') as out_html:
            out_html.write(template.render(output_data))
        print '...Done'
        print bar

#######################################

if __name__=='__main__':
    from giant.jiffies import run_default
    run_default(
        run                 = run,
        master_phil         = master_phil,
        args                = sys.argv[1:],
        blank_arg_prepend   = blank_arg_prepend,
        program             = PROGRAM,
        description         = DESCRIPTION)
