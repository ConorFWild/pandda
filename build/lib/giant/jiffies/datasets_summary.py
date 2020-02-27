import os, sys, copy, re

import libtbx.phil

from libtbx.utils import Sorry, Failure

import numpy
import pandas

from bamboo.common.logs import Log

from giant.xray.crystal import CrystalSummary
from giant.xray.crystal.cluster import CrystalGroup

from scitbx.array_family import flex
from scitbx.python_utils import robust_statistics

############################################################################

PROGRAM = 'giant.datasets.summary'

DESCRIPTION = """
    A tool to summarise the variation in sets of pdb and/or mtz files.
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
    file_label = *filename foldername
        .type = choice
}
check_for{
    column_label = None
        .type = str
        .multiple = True
}
summary {
    column_label = None
        .type = str
        .multiple = True
}
output {
    log_file = dataset_summary.log
        .type = path
}
""")

############################################################################

def format_column(col, col_width, position, fmt):
    COL_STR = '{:'+position+str(col_width)+fmt+'}'
    return COL_STR.format(col)

def format_row(cols, row='', col_width=20, header=False):
    if header: fmt = ''
    else:      fmt = '.3f'
    s = ('| '+format_column(row,col_width,'',''))*bool(row)+'|'+'|'.join([format_column(c,col_width,'^',fmt) for c in cols])+'|'
    if header:
        s = '|-'+('-'*col_width+'|')*(len(cols)+bool(row))+'\n'+s+'\n|-'+('-'*col_width+'|')*(len(cols)+bool(row))
    return s

def crystal_statistics(label, crystals, value_func, header=False, footer=False):
    sorted_crystals = sorted(crystals, key=value_func)
    sorted_values   = flex.double([value_func(c) for c in crystals])
    min_max_mean = sorted_values.min_max_mean()
    stddev       = sorted_values.sample_standard_deviation()
    hinges       = robust_statistics.hinges(sorted_values)
    median       = robust_statistics.median(sorted_values)
    s_out = []
    if header: s_out += [format_row(cols=['Min', 'Lower Quartile', 'Median', 'Upper Quartile', 'Max', 'Mean', 'Standard Dev'], row=' ',   header=True)]
    s_out += [format_row(cols=[min_max_mean.min, hinges[0], median, hinges[1], min_max_mean.max, min_max_mean.mean, stddev],   row=label, header=False)]
    if footer: s_out += [format_row(cols=['Min', 'Lower Quartile', 'Median', 'Upper Quartile', 'Max', 'Mean', 'Standard Dev'], row=' ',   header=True)]
    return '\n'.join(s_out)

def crystal_min_max(label, crystals, value_func):
    sorted_crystals = sorted(crystals, key=value_func)
    s_out = []
    s_out += [label+' Smallest: {}, {}, {}'.format(sorted_crystals[0].id, sorted_crystals[0].pdb_file, sorted_crystals[0].mtz_file)]
    s_out += [label+' Largest:  {}, {}, {}'.format(sorted_crystals[-1].id, sorted_crystals[-1].pdb_file, sorted_crystals[-1].mtz_file)]
    return '\n'.join(s_out)

############################################################################

def run(params):

    log = Log(log_file=params.output.log_file, verbose=True)

    # Process MTZs
    if params.input.mtz:

        log.heading('Processing {} MTZ Files'.format(len(params.input.mtz)))

        if   params.input.file_label=='filename':   labels = [os.path.basename(os.path.splitext(f)[0]) for f in params.input.mtz]
        elif params.input.file_label=='foldername': labels = [os.path.basename(os.path.dirname(f)) for f in params.input.mtz]
        else: raise Exception('MTZ labelling function not supported: {}'.format(params.input.file_label))

        log.bar()
        log('Grouping {} mtz files by space group'.format(len(params.input.mtz)))
        crystal_groups = CrystalGroup.by_space_group(crystals=[CrystalSummary.from_mtz(mtz_file=f, id=lab) for f,lab in zip(params.input.mtz, labels)])
        log('> Clustered into {} space group(s)'.format(len(crystal_groups)))
        log.bar()

        for cg in crystal_groups:

            log.subheading('Space group {} - {} datasets'.format(','.join(cg.space_groups), len(cg.crystals)))

            error = False
            for c in cg.crystals:
                for label in params.check_for.column_label:
                    if label is None: continue
                    if label not in c.column_labels:
                        log('Checking: column "{}" not in diffraction data of {}. columns present are {}'.format(label, c.mtz_file, c.column_labels))
                for label in params.summary.column_label:
                    if label is None: continue
                    if label not in c.column_labels:
                        log('Required: column "{}" not in diffraction data of {}. columns present are {}'.format(label, c.mtz_file, c.column_labels))
                        error = True
            if error is True: raise Sorry('There are datasets that do not contain the right columns.')

            log(crystal_statistics('Wavelength',         cg.crystals, value_func=lambda c: c.mtz_object().crystals()[1].datasets()[0].wavelength(), header=True))
            log(crystal_statistics('Resolution (high)',  cg.crystals, value_func=lambda c: c.high_res,                                              header=False))
            log(crystal_statistics('Resolution (low)',   cg.crystals, value_func=lambda c: c.low_res,                                               header=False))
            log(crystal_statistics('Unit cell - vol',    cg.crystals, value_func=lambda c: c.unit_cell.volume(),                                    header=False))
            log(crystal_statistics('Unit cell - a',      cg.crystals, value_func=lambda c: c.unit_cell.parameters()[0],                             header=False))
            log(crystal_statistics('Unit cell - b',      cg.crystals, value_func=lambda c: c.unit_cell.parameters()[1],                             header=False))
            log(crystal_statistics('Unit cell - c',      cg.crystals, value_func=lambda c: c.unit_cell.parameters()[2],                             header=False))
            log(crystal_statistics('Unit cell - alpha',  cg.crystals, value_func=lambda c: c.unit_cell.parameters()[3],                             header=False))
            log(crystal_statistics('Unit cell - beta',   cg.crystals, value_func=lambda c: c.unit_cell.parameters()[4],                             header=False))
            log(crystal_statistics('Unit cell - gamma',  cg.crystals, value_func=lambda c: c.unit_cell.parameters()[5],                             header=False, footer=True))

            for label in params.summary.column_label:
                if label is None: continue
                log(crystal_statistics('Column: {}'.format(label), cg.crystals, value_func=lambda c: c.mtz_object().get_column(label).n_valid_values(),     header=False, footer=True))

            log.bar(True, False)
            log('Smallest + Largest Values')
            log.bar()

            log(crystal_min_max('Resolution', cg.crystals, value_func=lambda c: c.high_res))

    # Process PDBs
    if params.input.pdb:

        log.heading('Processing {} PDB Files'.format(len(params.input.pdb)))

        if   params.input.file_label=='filename':   labels = [os.path.basename(os.path.splitext(f)[0]) for f in params.input.pdb]
        elif params.input.file_label=='foldername': labels = [os.path.basename(os.path.dirname(f)) for f in params.input.pdb]
        else: raise Exception('PDB labelling function not supported: {}'.format(params.input.file_label))

        log.bar()
        log('Grouping {} pdb files by space group'.format(len(params.input.pdb)))
        crystal_groups = CrystalGroup.by_space_group(crystals=[CrystalSummary.from_pdb(pdb_file=f, id=lab) for f,lab in zip(params.input.pdb, labels)])
        log('> Clustered into {} space group(s)'.format(len(crystal_groups)))

        for cg in crystal_groups:

            log.subheading('Space group: {} - {} datasets'.format(','.join(cg.space_groups), len(cg.crystals)))

            log(crystal_statistics('R-work', cg.crystals, value_func=lambda c: c.pdb_input().get_r_rfree_sigma().r_work, header=True))
            log(crystal_statistics('R-free', cg.crystals, value_func=lambda c: c.pdb_input().get_r_rfree_sigma().r_free, header=False, footer=True))

            log.bar(True, False)
            log('Smallest + Largest Values')
            log.bar()

            log(crystal_min_max('R-free',     cg.crystals, value_func=lambda c: c.pdb_input().get_r_rfree_sigma().r_free))

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
