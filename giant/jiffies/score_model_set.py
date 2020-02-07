import os, sys, copy, re, glob

import libtbx.phil
import libtbx.easy_mp

import numpy, pandas

from bamboo.html import BAMBOO_HTML_ENV
from bamboo.stats.cluster import generate_group_idxs
from giant.jiffies.score_model import prepare_output_directory, score_model, make_residue_radar_plot, format_parameters_for_plot

from matplotlib import pyplot
pyplot.rc('axes', color_cycle=['r', 'g', 'b', 'y'])

#######################################

bar = '=======================++>'

PROGRAM = 'giant.score_model_set'
DESCRIPTION = """
    A tool to quickly score a set of related ligand models against crystallographic electron density
    The radar plots for the ligands will be output on the same graph for comparison

    1) Simple usage:
        (mtz files must be named the same as the pdb files, e.g. structure_1.mtz)
        > giant.score_model structure_1.pdb structure_2.pdb structure_2.pdb

    See giant.score_model for more information.
"""

blank_arg_prepend = {None:'pdb='}

master_phil = libtbx.phil.parse("""
input {
    pdb = None
        .type = path
        .multiple = True
    ref_pdb = None
        .type = path
        .multiple = False
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
radar_plot {
    limits = automatic *manual
        .help = 'Select how the axis limits are determined'
        .type = choice
        .multiple = False
}
include scope giant.jiffies.score_model.residue_plot_phil
""", process_includes=True)

#######################################

def run(params):

    assert params.input.pdb, 'No pdb files provided'
    # REMOVE THIS WHEN IMPLEMENTED
    assert params.selection.res_names
    # REMOVE THIS WHEN IMPLEMENTED
    if params.selection.res_names:      params.selection.__inject__("res_names_list", params.selection.res_names.split(','))
    else:                               params.selection.__inject__("res_names_list", None)

    output_dir, images_dir = prepare_output_directory(params)
    scores_file =  os.path.join(output_dir, 'residue_scores.csv')

    all_data = pandas.DataFrame()

    for pdb in params.input.pdb:

        mtz   = pdb.replace('.pdb','.mtz')

        if   params.input.labels == 'basename':    label = os.path.split_ext(os.path.basename(pdb))[0]
        elif params.input.labels == 'folder_name': label = os.path.basename(os.path.dirname(os.path.abspath(pdb)))

        print bar
        print 'Scoring model {} against {}'.format(pdb, mtz)

        data_table = score_model(   params = params,
                                    pdb1   = pdb,
                                    mtz1   = mtz,
                                    pdb2   = params.input.ref_pdb,
                                    label_prefix = label
                                )
        all_data = all_data.append(data_table, verify_integrity=True)

    print '...Done'
    print bar

    all_data.to_csv(scores_file)
    print 'Output written to {}'.format(scores_file)
    print bar

    ###################################################################
    # Image parameters
    ###################################################################
    columns = format_parameters_for_plot(params=params.plot.parameters)

    ###################################################################
    # Output Images - 1 image per residue per structure
    ###################################################################
    all_images = []
    print 'Generating Output Images...'
    for label, row in all_data.iterrows():
        image_path = os.path.join(images_dir,'{}.png'.format(label))
        print 'Making: {}...'.format(image_path)
        make_residue_radar_plot(path = image_path,
                                data = row.to_frame().T,
                                columns = columns,
                                remove_blank_entries = params.plot.remove_blank_entries,
                                print_axis_values    = params.plot.print_axis_values  )
        all_images.append(image_path)

    ###################################################################
    # Output Images - 1 image per residue (allowing comparisons)
    ###################################################################
    if params.radar_plot.limits == 'automatic': columns.pop('limits', None)
    elif params.radar_plot.limits == 'manual':  pass

    for res_label, index_idxs in generate_group_idxs([i.split('-')[-2:] for i in all_data.index]):
        res_label = '-'.join(res_label)
        image_path = os.path.join(images_dir,'compare-{}.png'.format(res_label))
        print 'Making: {}...'.format(image_path)
        make_residue_radar_plot(path = image_path,
                                data = all_data.iloc[index_idxs],
                                columns = columns,
                                remove_blank_entries = params.plot.remove_blank_entries,
                                print_axis_values    = params.plot.print_axis_values  )

    print '...Done.'
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
