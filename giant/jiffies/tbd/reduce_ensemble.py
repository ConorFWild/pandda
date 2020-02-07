import os, sys, glob, re

#################################
import matplotlib
matplotlib.use('Agg')
matplotlib.interactive(0)
from matplotlib import pyplot
pyplot.style.use('ggplot')
#################################

import libtbx.phil

import numpy, pandas

import iotbx.pdb

from scitbx.array_family import flex

from giant.structure.ensembles import StructureCollection

#######################################

bar = '=======================++>'

#######################################

blank_arg_prepend = {'.pdb':'pdb='}

master_phil = libtbx.phil.parse("""
input {
    pdb = None
        .help = 'Input structure files'
        .type = str
        .multiple = True
    label = file_name *folder_name
        .help = 'How should the datasets be labelled?'
        .type = choice
        .multiple = False
}
output {
    out_dir = reduced_ensemble
        .type = path
}
""")


#######################################

def run(params):

    assert params.input.pdb, 'No PDBs provided'

    # Sort input files
    params.input.pdb = sorted(params.input.pdb)

    if params.input.label == 'folder_name': labels = [os.path.basename(os.path.dirname(p)) for p in params.input.pdb]
    elif params.input.label == 'file_name': labels = [os.path.basename(p) for p in params.input.pdb]

    e = StructureCollection.from_files(params.input.pdb, labels=labels)
    e.load_all()
    e.write_all(out_dir=params.output.out_dir)

#######################################

if __name__=='__main__':
    from giant.jiffies import run_default
    run_default(run=run, master_phil=master_phil, args=sys.argv[1:], blank_arg_prepend=blank_arg_prepend)
