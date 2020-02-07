import os, sys, copy, re, shutil

import libtbx.phil
from bamboo.common.command import CommandManager

#######################################

blank_arg_prepend = {'.mtz':'input.mtz=', None:'column.label='}

master_phil = libtbx.phil.parse("""
input {
    mtz = None
        .multiple = True
        .type = path
}
output {
    mtz = merged.mtz
        .type = path
    log = merged.log
        .type = path
}
options {
    label_suffix = *incremental filename foldername
        .help = 'suffix added to the column name in the output file'
        .type = choice
        .multiple = False
    column {
        label = None
            .help = 'columns to be extracted from the input mtz files'
            .type = str
            .multiple = True
    }
}
""")

#######################################

def run(params):

    in_mtzs = params.input.mtz
    assert len(params.input.mtz) > 1, 'Need to provide at least one mtz files'
#    assert len(params.input.mtz) < 3, 'Need to provide at most two mtz files'

    out_mtz = params.output.mtz
    out_log = params.output.log

    ###########################################
    # COMMAND LINE COMMANDS
    ###########################################
    cm = CommandManager('cad')

    ###########################################
    # ITERATE THROUGH INPUT MTZS
    ###########################################
    for i_mtz, mtz in enumerate(in_mtzs):
        # Use numbering from 1
        n_mtz = i_mtz + 1
        # Create an mtz suffix
        if   params.options.label_suffix == 'incremental': suffix = '-{}'.format(n_mtz)
        elif params.options.label_suffix == 'filename':    suffix = '-{}'.format(os.path.splitext(os.path.basename(mtz))[0])
        elif params.options.label_suffix == 'foldername':  suffix = '-{}'.format(os.path.basename(os.path.dirname(mtz)))
        else:                                              raise Exception('Not yet implemented, sorry')
        # Add to the command line
        cm.add_command_line_arguments(['hklin{}'.format(n_mtz), mtz])
        # Select column labels
        cm.add_standard_input('labin file_number {0} {1}'.format(n_mtz, ' '.join(['E{0}={1}'.format(i+1,c) for i,c in enumerate(params.options.column.label)])))
        cm.add_standard_input('labout file_number {0} {1}'.format(n_mtz, ' '.join(['E{0}={1}'.format(i+1,c+suffix) for i,c in enumerate(params.options.column.label)])))

    ###########################################
    # OUTPUT FILES
    ###########################################
    cm.add_command_line_arguments(['hklout', out_mtz])

    ###########################################
    # RUN
    ###########################################
    print 'Running CAD:'
    cm.print_settings()
    ret_code = cm.run()
    cm.write_output(log_file=out_log)
    if ret_code != 0:
        print '============================>'
        print 'Refmac returned with an error'
        print '============================>'
        print cm.output
        print '============================>'
        print cm.error
        print '============================>'

#######################################

if __name__=='__main__':
    from giant.jiffies import run_default
    run_default(run=run, master_phil=master_phil, args=sys.argv[1:], blank_arg_prepend=blank_arg_prepend)
