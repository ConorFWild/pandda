import os, sys, copy, re, glob

#################################
import matplotlib
matplotlib.use('Agg')
matplotlib.interactive(0)
from matplotlib import pyplot
pyplot.style.use('ggplot')
#################################

import libtbx.phil

import numpy, pandas

from bamboo.plot import Radar

from giant.io.pdb import strip_pdb_to_input
from giant.xray.edstats import Edstats
from giant.structure import sanitise_hierarchy, calculate_residue_group_occupancy, calculate_paired_conformer_rmsds
from giant.structure.b_factors import calculate_residue_group_bfactor_ratio
from giant.structure.formatting import ShortLabeller
from giant.structure.select import non_h, protein, backbone, sidechains

#######################################

bar = '=======================++>'

PROGRAM = 'giant.score_model'
DESCRIPTION = """
    A tool to quickly score a ligand model against crystallographic electron density

    1) Simple usage (for ligand called LIG or UNL):
        > giant.score_model refined.pdb refined.mtz

    2) Score a residue in the file (replace XXX with ligand 3-letter ID)
        > giant.score_model ... res_names=XX1,XX2,XX3

    3) Define a "fitted" model to compare the refined model against
        > giant.score_model ... pdb1=fitted.pdb mtz1=fitted.mtz
"""

blank_arg_prepend = {'.pdb':'pdb1=', '.mtz':'mtz1='}

residue_plot_phil =  """
plot {
    remove_blank_entries = False
        .type = bool
    print_axis_values = True
        .type = bool
    parameters {
        rscc {
            title = 'Model\nQuality\n(RSCC)'
                .type = str
            axis_min = 0.60
                .type = float
            axis_max = 0.85
                .type = float
            axis_invert = True
                .type = bool
        }
        rszd {
            title = 'Model\nAccuracy\n(RSZD)'
                .type = str
            axis_min = 1.50
                .type = float
            axis_max = 4.00
                .type = float
            axis_invert = False
                .type = bool
        }
        rszo {
            title = 'Density\nPrecision\n(RSZO/OCC)'
                .type = str
            axis_min = 0.00
                .type = float
            axis_max = 2.00
                .type = float
            axis_invert = True
                .type = bool
        }
        b_factor_ratio {
            title = 'B-Factor\nStability\n(B-factor Ratio)'
                .type = str
            axis_min = 1.00
                .type = float
            axis_max = 3.00
                .type = float
            axis_invert = False
                .type = bool
        }
        rmsd {
            title = 'Coordinate\nStability\n(RMSD)'
                .type = str
            axis_min = 0.00
                .type = float
            axis_max = 1.50
                .type = float
            axis_invert = False
                .type = bool
        }
    }
}
"""

master_phil = libtbx.phil.parse("""
input {
    pdb1 = None
        .type = str
        .multiple = False
    mtz1 = None
        .type = str
        .multiple = False

    pdb2 = None
        .type = str
        .multiple = False
    mtz2 = None
        .type = str
        .multiple = False

    label = ''
        .type = str
        .multiple = False
}
selection {
    res_names = LIG,UNL,DRG
        .type = str
        .help = "Comma-separated list of residue names to score -- if None then scores all residues"
}
output {
    out_dir = ./
        .type = path
}
include scope giant.phil.settings_phil
"""+residue_plot_phil, process_includes=True)

#######################################

def prepare_table():

    columns_fix = [ 'Model RMSD' ]
    columns_num = [ 'RSCC','RSZD','RSZO','RSR',
                    'RSZO/OCC', 'Occupancy',
                    'Surroundings B-factor Ratio',
                    'Surroundings Residue Labels',
                    'Average B-factor (Residue)',
                    'Average B-factor (Surroundings)' ]
    columns_end = [ 'PDB','MTZ' ]

    columns = []
    columns.extend(columns_fix)
    columns.extend(columns_num)
    #columns.extend([c+'-2' for c in columns_num])
    columns.extend(columns_end)
    columns.extend([c+'-2' for c in columns_end])

    return pandas.DataFrame(index=[], columns=columns)

def prepare_output_directory(params):
    if not os.path.exists(params.output.out_dir): os.mkdir(params.output.out_dir)
    images_dir = os.path.join(params.output.out_dir, 'residue_plots')
    if not os.path.exists(images_dir): os.mkdir(images_dir)
    return params.output.out_dir, images_dir

def score_model(params, pdb1, mtz1, pdb2=None, mtz2=None, label_prefix='', verbose=False):
    """
    Score residues against density, and generate other model quality indicators.
    Identified residues in pdb1 are scored against mtz1 (and mtz2, if provided) using edstats.
    Identified residues in pdb1 are compared to the equivalent residues in pdb2, if provided.
    B-factors ratios of identified residues to surrounding sidechains are calculated.
    """

    if label_prefix: label_prefix = label_prefix + '-'

    # Extract the residues to look for
    res_names = params.selection.res_names_list

    print 'Reading input structure:', pdb1

    # Extract Structure
    h1_all = non_h(strip_pdb_to_input(pdb1, remove_ter=True, remove_end=True).hierarchy)
    # Normalise hierarchy (standardise atomic naming, etc...)
    sanitise_hierarchy(h1_all)
    h1_pro = protein(h1_all)
    h1_bck = backbone(h1_all)
    h1_sch = sidechains(h1_all)

    # Pull out residues to analyse
    if res_names:
        rg_for_analysis = [rg for rg in h1_all.residue_groups() if [n for n in rg.unique_resnames() if n in res_names]]
        print 'Selecting residues named {}: {} residue(s)'.format(' or '.join(res_names), len(rg_for_analysis))
    else:
        rg_for_analysis = h1_all.residue_groups()
        print 'Analysing all residues ({} residues)'.format(len(rg_for_analysis))

    # Check residues to analyse or skip
    if not rg_for_analysis:
        raise Exception('There are no residues called {} in {}'.format(' or '.join(params.selection.res_names_list), pdb1))

    # Extract PDB2
    if pdb2 is not None:
        print 'Reading input structure:', pdb2
        h2_all = non_h(strip_pdb_to_input(pdb2, remove_ter=True, remove_end=True).hierarchy)
        sanitise_hierarchy(h2_all)

    # Score MTZ1
    if mtz1 is not None:
        print 'Scoring model against mtz file'
        print 'Scoring {} >>> {}'.format(pdb1, mtz1)
        mtz1_edstats_scores = Edstats(mtz_file=mtz1, pdb_file=pdb1)
    else:
        mtz1_edstats_scores = None
    # Score MTZ2
    if mtz2 is not None:
        print 'Scoring model against mtz file'
        print 'Scoring {} >>> {}'.format(pdb1, mtz2)
        mtz2_edstats_scores = Edstats(mtz_file=mtz2, pdb_file=pdb1)
    else:
        mtz2_edstats_scores = None

    # Prepare output table
    data_table = prepare_table()

    for rg_sel in rg_for_analysis:

        # Create label for the output table
        #rg_label = (label_prefix+rg_sel.unique_resnames()[0]+'-'+rg_sel.parent().id+'-'+rg_sel.resseq+rg_sel.icode).replace(' ','')
        #rg_label = (label_prefix+rg_sel.parent().id+'-'+rg_sel.resseq+rg_sel.icode).replace(' ','')
        rg_label = ShortLabeller.format(rg_sel).replace(' ','')
        tab_label = label_prefix + rg_label

        if len(rg_sel.unique_resnames()) != 1:
            raise Exception(tab_label+': More than one residue name associated with residue group -- cannot process')

        # Append empty row to output table
        data_table.loc[tab_label] = None

        data_table.set_value(index = tab_label,
                             col   = 'PDB',
                             value = pdb1 )
        data_table.set_value(index = tab_label,
                             col   = 'Occupancy',
                             value = calculate_residue_group_occupancy(residue_group=rg_sel) )

        data_table = calculate_residue_group_bfactor_ratio(residue_group = rg_sel,
                                                           hierarchy     = h1_sch,
                                                           data_table    = data_table,
                                                           rg_label      = tab_label)

        if pdb2 is not None:
            data_table.set_value(index = tab_label,
                                 col   = 'PDB-2',
                                 value = pdb2 )

            # Extract the equivalent residue in pdb2
            rg_sel_2 = [rg for rg in h2_all.residue_groups() if ShortLabeller.format(rg).replace(' ','') == rg_label]

            try:
                assert rg_sel_2, 'Residue is not present in pdb file: {} not in {}'.format(rg_label, pdb2)
                assert len(rg_sel_2) == 1, 'More than one residue has been selected for {} in {}'.format(rg_label, pdb2)
            except:
                raise

            # Extract occupancy
            data_table.set_value(index = tab_label,
                                 col   = 'Occupancy-2',
                                 value = calculate_residue_group_occupancy(residue_group=rg_sel_2[0]) )

            # Calculate the RMSD between the models
            try:
                confs1, confs2, rmsds = zip(*calculate_paired_conformer_rmsds(conformers_1=rg_sel.conformers(), conformers_2=rg_sel_2[0].conformers()))
                data_table.set_value(index=tab_label, col='Model RMSD', value=min(rmsds))
            except:
                raise
                print 'Could not calculate RMSD between pdb_1 and pdb_2 for residue {}'.format(rg_label)
                pass

        # Extract Density Scores - MTZ 1
        if mtz1 is not None:
            data_table.set_value(index=tab_label, col='MTZ', value=mtz1)
        if mtz1_edstats_scores is not None:
            data_table = mtz1_edstats_scores.extract_residue_group_scores(  residue_group  = rg_sel,
                                                                            data_table     = data_table,
                                                                            rg_label       = tab_label )
            # Normalise the RSZO by the Occupancy of the ligand
            data_table['RSZO/OCC'] = data_table['RSZO']/data_table['Occupancy']

        # Extract Density Scores - MTZ 2
        if mtz2 is not None:
            data_table.set_value(index=tab_label, col='MTZ-2', value=mtz2)
        if mtz2_edstats_scores is not None:
            data_table = mtz2_edstats_scores.extract_residue_group_scores(  residue_group  = rg_sel,
                                                                            data_table     = data_table,
                                                                            rg_label       = tab_label,
                                                                            column_suffix  = '-2' )
            # Normalise the RSZO by the Occupancy of the ligand
            data_table['RSZO/OCC-2'] = data_table['RSZO-2']/data_table['Occupancy-2']

    return data_table

def make_residue_radar_plot(path, data, columns, linetype=None, remove_blank_entries=False, print_axis_values=True):
    "Plot radar graph. data is a list of (label, scores) tuples. label is a string. scores is a pandas Series."

#    column_limit = [(0.6,0.85),(1.5,4),(0,2),(1,3),(0,1.5)]
#    column_names = ['RSCC','RSZD','RSZO/OCC','Surroundings B-factor Ratio','Model RMSD']
#    column_title = ['Model\nQuality\n(RSCC)', 'Model\nAccuracy\n(RSZD)', 'Model\nPrecision\n(RSZO/OCC)', 'B-Factor\nRatio', 'Model\nRMSD']
#    column_invse = [1,0,1,0,0]

    # ----------------------->
    assert isinstance(columns, dict)
    # ----------------------->
    # Column Titles are compulsory
    column_title = columns['titles']
    # ----------------------->
    # Column names are compulsory (for pulling from data_frame)
    column_names = columns['names']
    assert len(column_names) == len(column_title)
    # ----------------------->
    # Limits are optional
    if 'limits' in columns:
        column_limit = columns['limits']
        assert len(column_limit) == len(column_title)
    else:
        column_limit = None
    # ----------------------->
    # Inverse is optional
    if 'invert' in columns:
        column_invse = columns['invert']
        assert len(column_invse) == len(column_title)
    else:
        column_invse = None

    # ----------------------->
    # Extract the plot data from the data_frame
    plot_data = data[column_names]

    if remove_blank_entries:
        # ----------------------->
        # Filter the entries based on whether there is at least one values for each column
        data_mask = [True if data[c].any() else False for c in column_names]
        # ----------------------->
        # Filter against the mask
        column_title = [column_title[i] for i in range(len(data_mask)) if data_mask[i]]
        column_names = [column_names[i] for i in range(len(data_mask)) if data_mask[i]]
        if column_invse: column_invse = [column_invse[i] for i in range(len(data_mask)) if data_mask[i]]
        if column_limit: column_limit = [column_limit[i] for i in range(len(data_mask)) if data_mask[i]]
        # ----------------------->
        # Reselect the plot_data
        plot_data = data[column_names]

    # ----------------------->
    # Round column values
    col_tick_vals = []
    col_tick_labs = []
    for col in column_names:
        tvs=[]; tls=[]
        for i, v in enumerate(plot_data[col]):
            try: l = round(v, 2)
            except:
                l = 'Error'; v = None;
            # Add point marker
            tvs.append(v); tls.append(l)
        # Add to column ticks
        col_tick_vals.append(tvs); col_tick_labs.append(tls)

    # ----------------------->
    # Prepare markers
    # ----------------------->
    if plot_data.index.size > 1:
        mrs = ['o','^','s','D','*']
        markers = mrs*(1+(plot_data.index.size-1)//len(mrs))
        markersize = 10
    else:
        markers = [None]
        markersize = None

    # ----------------------->
    # Prepare linetypes
    # ----------------------->
    if not linetype:
        if plot_data.index.size > 1:
            #lts = ['-']*len(mrs)+['--']*len(mrs)
            #lts = ['--','-','-.',':']
            lts = ['-','-','--']
        else:
            lts = ['-']
        linetype = lts*(1+(plot_data.index.size-1)//len(lts))
    elif not isinstance(linetype, list):
        linetype = [linetype]*plot_data.index.size
    assert len(linetype) >= plot_data.index.size

    # ----------------------->
    r = Radar(titles=column_title)
    # ----------------------->
    # Add each row of the data frame as a separate line
    for label, row in plot_data.iterrows():
        r.add(row[column_names].tolist(), linetype.pop(0), marker=markers.pop(0), markersize=markersize, label=label, lw=3)
    # ----------------------->
    # Set axis meta manually
    if column_invse: r.set_inversion(column_invse)
    if column_limit: r.set_limits(column_limit)
    # ----------------------->
    if print_axis_values: r.set_ticks(values=col_tick_vals, labels=col_tick_labs)
    else:                 r.set_ticks(values=col_tick_vals, labels=[' ']*len(col_tick_vals))
    # ----------------------->
    # Plot, modify and save
    r.plot()
    #r.ax.legend(loc='lower center', fancybox=True, bbox_to_anchor=(1.0, 1.0), fontsize=20)
    r.ax.patch.set_facecolor('lightgray')
    r.ax.legend(loc='lower left', bbox_to_anchor=(0.65, 0.95), bbox_transform=pyplot.gcf().transFigure, fontsize=20)
    r.savefig(path)
    r.close()

    return

def format_parameters_for_plot(params):
    """Convert plot scope parameters to parameter dict"""

    p = params
    columns = {}
    columns['titles'] = [p.rscc.title, p.rszd.title, p.rszo.title, p.b_factor_ratio.title, p.rmsd.title]
    columns['names']  = ['RSCC','RSZD','RSZO/OCC','Surroundings B-factor Ratio','Model RMSD']
    columns['invert'] = [p.rscc.axis_invert, p.rszd.axis_invert, p.rszo.axis_invert, p.b_factor_ratio.axis_invert, p.rmsd.axis_invert ]
    columns['limits'] = [(p.rscc.axis_min,           p.rscc.axis_max),
                         (p.rszd.axis_min,           p.rszd.axis_max),
                         (p.rszo.axis_min,           p.rszo.axis_max),
                         (p.b_factor_ratio.axis_min, p.b_factor_ratio.axis_max),
                         (p.rmsd.axis_min,           p.rmsd.axis_max)  ]

    return columns

#######################################

def run(params):

    assert params.input.pdb1, 'No pdb1 provided'
    assert params.input.mtz1, 'No mtz1 provided'
    # REMOVE THIS WHEN IMPLEMENTED
    assert params.selection.res_names
    # REMOVE THIS WHEN IMPLEMENTED
    if params.selection.res_names:      params.selection.__inject__("res_names_list", params.selection.res_names.split(','))
    else:                               params.selection.__inject__("res_names_list", None)

    output_dir, images_dir = prepare_output_directory(params)
    scores_file = os.path.join(output_dir, 'residue_scores.csv')

    print bar
    print 'Scoring model...'
    data_table = score_model(   params = params,
                                pdb1   = params.input.pdb1,
                                mtz1   = params.input.mtz1,
                                pdb2   = params.input.pdb2,
                                mtz2   = params.input.mtz2,
                                label_prefix = params.input.label,
                                verbose = params.settings.verbose
                            )
    print '...Done'
    print bar

    data_table.dropna(axis=1, how='all').to_csv(scores_file)
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
        print 'Making: {}...'.format(label)
        image_path = os.path.join(images_dir,'{}.png'.format(label))
        make_residue_radar_plot(path = image_path,
                                data = row.to_frame().T,
                                columns = columns,
                                remove_blank_entries = params.plot.remove_blank_entries,
                                print_axis_values    = params.plot.print_axis_values  )
        all_images.append(image_path)
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
