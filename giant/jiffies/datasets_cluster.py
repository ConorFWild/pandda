import os, sys, copy, re, shutil

import libtbx.phil

from libtbx.utils import Sorry, Failure

import numpy
import pandas

from bamboo.common.logs import Log

from giant.xray.crystal import CrystalSummary
from giant.xray.crystal.cluster import CrystalGroup

try:
    import matplotlib
    matplotlib.interactive(False)
    from matplotlib import pyplot
    pyplot.switch_backend('agg')
    pyplot.style.use('ggplot')
except:
    raise Sorry('Failed to import matplotlib -- needed for plotting')

############################################################################

PROGRAM = 'giant.datasets.cluster'

DESCRIPTION = """
    A tool to cluster sets of pdb and/or mtz files.

    1) Simple usage:
        > giant.datasets.cluster *.pdb pdb_label=filename
    or
        > giant.datasets.cluster */model.pdb pdb_label=foldername

"""

############################################################################

blank_arg_prepend = {'.mtz':'mtz=', '.pdb':'pdb='}

master_phil = libtbx.phil.parse("""
input {
    mtz = None
        .type = path
        .multiple = True
    pdb = None
        .type = path
        .multiple = True
    labels {
        mtz_label = None* filename foldername
            .type = choice
        mtz_regex = None
            .type = str
        pdb_label = None* filename foldername
            .type = choice
        pdb_regex = None
            .type = str
    }
}
clustering{
    lcv_cutoff = 0.2
        .type = float
    label_nodes_above = 0.2
        .help = "Label nodes above this value"
        .type = float
}
output {
    out_dir = clustered-datasets
        .type = path
    file_mode = copy *symlink
        .type = choice
    split_pdbs_and_mtzs = False
        .type = bool
    html_out = None
        .type = path
}
""")

############################################################################

def run(params):

    # Validate input files
    if not (params.input.pdb or params.input.mtz):
        raise Sorry('No pdb/mtz files have been provided: specify with input.pdb or input.mtz')
    # Check and create output directory
    if not params.output.out_dir:
        raise Sorry('No output directory has been specified: specify with output.out_dir')
    if not os.path.exists(params.output.out_dir):
        os.mkdir(params.output.out_dir)
    # Define and create image directory
    img_dir = os.path.join(params.output.out_dir, 'dendrograms')
    if not os.path.exists(img_dir):
        os.mkdir(img_dir)

    # Create log object
    log = Log(log_file=params.output.out_dir+'.clustering.log', verbose=True)

    # Define output_file_function to copy or symlink files as needed
    if params.output.file_mode == 'symlink':
        out_file_func = os.symlink
    elif params.output.file_mode == 'copy':
        out_file_func = shutil.copy

    log.heading('Processing input pdb/mtz files')
    log('Making dataset labels for {} pdb(s) and {} mtz(s)'.format(len(params.input.pdb), len(params.input.mtz)))

    try:
        if   params.input.labels.pdb_label=='filename':   p_labels = [os.path.basename(os.path.splitext(f)[0]) for f in params.input.pdb]
        elif params.input.labels.pdb_label=='foldername': p_labels = [os.path.basename(os.path.dirname(f))     for f in params.input.pdb]
        elif params.input.labels.pdb_regex:               p_labels = [re.findall(params.input.labels.pdb_regex, f)[0] for f in params.input.pdb]
        else:                                             p_labels = ['PDB-{:06d}'.format(i) for i in range(len(params.input.pdb))]
        if   params.input.labels.mtz_label=='filename':   m_labels = [os.path.basename(os.path.splitext(f)[0]) for f in params.input.mtz]
        elif params.input.labels.mtz_label=='foldername': m_labels = [os.path.basename(os.path.dirname(f))     for f in params.input.mtz]
        elif params.input.labels.mtz_regex:               m_labels = [re.findall(params.input.labels.mtz_regex, f)[0] for f in params.input.mtz]
        else:                                             m_labels = ['MTZ-{:06d}'.format(i) for i in range(len(params.input.mtz))]
    except:
        print 'Error reading file: {}'.format(f)
        raise

    # Check labels are unique
    set_m_labels = set(m_labels)
    set_p_labels = set(p_labels)
    if len(set_m_labels) != len(m_labels):
        raise Sorry('MTZ labels are not unique. Repeated labels: {}'.format(' '.join(['{}'.format(l) for l in set_m_labels if m_labels.count(l)!=1])))
    if len(set_p_labels) != len(p_labels):
        raise Sorry('PDB labels are not unique. Repeated labels: {}'.format(' '.join([l for l in set_p_labels if p_labels.count(l)!=1])))

    # Report labels
    if p_labels:
        log.subheading('PDB Labels')
        log(', '.join(p_labels))
    if m_labels:
        log.subheading('MTZ Labels')
        log(', '.join(m_labels))

    # Load crystal summaries
    log.bar(True, True)
    log('Reading data for {} pdb(s) and {} mtz(s)'.format(len(params.input.pdb), len(params.input.mtz)))

    if params.input.pdb: pdb_summaries = [CrystalSummary.from_pdb(pdb_file=f, id=lab) for f,lab in zip(params.input.pdb, p_labels)]
    else:                pdb_summaries = []
    if params.input.mtz: mtz_summaries = [CrystalSummary.from_mtz(mtz_file=f, id=lab) for f,lab in zip(params.input.mtz, m_labels)]
    else:                mtz_summaries = []

    # Group by SpaceGroup
    log.subheading('Grouping {} crystals by space group...'.format(len(pdb_summaries+mtz_summaries)))
    crystal_groups = CrystalGroup.by_space_group(crystals=pdb_summaries+mtz_summaries)
    log('Grouped crystals into {} space groups'.format(len(crystal_groups)))

    log.heading('Analysing variation of unit cells for each space group')

    for cg in crystal_groups:

        sg_name = 'sg-{}'.format(cg.space_groups[0].split(' (')[0].replace(' ','_'))

        log.subheading('Space Group {}: {} dataset(s)'.format(cg.space_groups[0], len(cg.crystals)))

        log('Unit Cell Variation:')
        log(numpy.round(cg.uc_stats.as_pandas_table().T, 2))

        log('')
        log('Making unit cell dendrogram for all crystals with this spacegroup')
        if len(cg.crystals) > 1:
            cg.dendrogram(  fname = os.path.join(img_dir,'{}-all.png'.format(sg_name)),
                            xlab  = 'Crystal',
                            ylab  = 'Linear Cell Variation',
                            annotate_y_min = params.clustering.label_nodes_above)

        log('')
        log('Clustering {} unit cells...'.format(len(cg.crystals)))
        sg_crystal_groups = cg.by_unit_cell(cg.crystals, cutoff=params.clustering.lcv_cutoff)
        log('Clustered crystals into {} groups'.format(len(sg_crystal_groups)))

        for i_cg2, cg2 in enumerate(sg_crystal_groups):

            cluster_name = '{}-cluster-{}'.format(sg_name, i_cg2+1)

            log.bar(True, False)
            log('Processing cluster: {}'.format(cluster_name))
            log.bar(False, True)

            log('Unit Cell Variation:')
            log(numpy.round(cg.uc_stats.as_pandas_table().T, 2))

            log('')
            log('Making unit cell dendrogram for this cluster of crystals')
            if len(cg2.crystals) > 1:
                cg2.dendrogram( fname = os.path.join(img_dir, '{}.png'.format(cluster_name)),
                                xlab  = 'Crystal',
                                ylab  = 'Linear Cell Variation',
                                ylim  = (0, params.clustering.lcv_cutoff),
                                annotate_y_min = params.clustering.label_nodes_above)

            log('Copying files to output directory')

            # Go through and link the datasets for each of the spacegroups into a separate folder
            sub_dir = os.path.join(params.output.out_dir, cluster_name)
            if not os.path.exists(sub_dir): os.mkdir(sub_dir)

            # Split the mtzs and pdbs into separate directories -- or not
            if params.output.split_pdbs_and_mtzs:
                mtz_dir = os.path.join(sub_dir, 'mtzs')
                if not os.path.exists(mtz_dir): os.mkdir(mtz_dir)
                pdb_dir = os.path.join(sub_dir, 'pdbs')
                if not os.path.exists(pdb_dir): os.mkdir(pdb_dir)
            else:
                mtz_dir = pdb_dir = sub_dir

            for c in cg2.crystals:
                # Set parameters based on pdb or mtz
                if c.mtz_file:
                    sub_sub_dir = os.path.join(mtz_dir, c.id)
                    def_file = os.path.abspath(c.mtz_file)
                    def_suff = '.mtz'
                    pos_suff = '.pdb'
                elif c.pdb_file:
                    sub_sub_dir = os.path.join(pdb_dir, c.id)
                    def_file = os.path.abspath(c.pdb_file)
                    def_suff = '.pdb'
                    pos_suff = '.mtz'
                # Create subdirectory
                if not os.path.exists(sub_sub_dir): os.mkdir(sub_sub_dir)
                # Output file base template
                out_base = os.path.join(sub_sub_dir, c.id)
                # Export file
                out_file = out_base+def_suff
                if not os.path.exists(out_file):
                    out_file_func(def_file, out_file)
                # output other as well if filenames are the same
                pos_file = def_file.replace(def_suff,pos_suff)
                out_file = out_base+pos_suff
                if os.path.exists(pos_file) and not os.path.exists(out_file):
                    out_file_func(pos_file, out_file)

    log.heading('finished')

############################################################################

if __name__=='__main__':
    from giant.jiffies import run_default
    run_default(
        run                 = run,
        master_phil         = master_phil,
        args                = sys.argv[1:],
        blank_arg_prepend   = blank_arg_prepend,
        program             = PROGRAM,
        description         = DESCRIPTION)
