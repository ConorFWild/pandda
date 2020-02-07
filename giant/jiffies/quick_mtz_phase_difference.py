import os, sys, glob

import libtbx.phil

from libtbx.utils import Sorry, Failure

from bamboo.common.command import CommandManager
from bamboo.common.path import rel_symlink

############################################################################

PROGRAM = 'giant.calculate_phase_differences'

DESCRIPTION = """
Calculate the average phase difference between sets of mtz files using cphasematch from CCP4.
"""

############################################################################

blank_arg_prepend = {   '.mtz':'mtz='   }

master_phil = libtbx.phil.parse("""
input {
    mtz = 'refine.mtz'
        .multiple = True
        .type = str
}
output {
    merged_mtz = merged.mtz
        .type = str
    merging_log = phase_difference_merging.log
        .type = str
    phasematch_mtz = phasematch.mtz
        .type = str
    phasematch_log_template = phase_difference
        .type = str

}
settings {
    f_obs = F,SIGF
        .type = str
    f_calc = F-model,PHIF-model
        .type = str
}

""")

############################################################################

def run(params):

    assert params.input.mtz, 'No MTZs given for comparison'

    assert not os.path.exists(params.output.merged_mtz), 'The output file ({}) already exists. Please delete it before re-running.'.format(params.output.merged_mtz)
    assert not os.path.exists(params.output.phasematch_mtz), 'The output file ({}) already exists. Please delete it before re-running.'.format(params.output.phasematch_mtz)

    merge = CommandManager('giant.mtz.merge')
    merge.add_command_line_arguments( 'label_suffix=incremental' )
    merge.add_command_line_arguments( params.input.mtz )
    merge.add_command_line_arguments( 'output.mtz='+params.output.merged_mtz )
    merge.add_command_line_arguments( params.settings.f_obs.split(',') )
    merge.add_command_line_arguments( params.settings.f_calc.split(',') )
    merge.run()
    merge.write_output(params.output.merging_log)

    if not os.path.exists(params.output.merged_mtz):
        raise Failure('giant.mtz.merge has failed to merge the mtz files')

    fo1, fo2 = params.settings.f_obs.split(',')
    fc1, fc2 = params.settings.f_calc.split(',')

    for i_2 in range(1, len(params.input.mtz)+1):
        for i_1 in range(1, len(params.input.mtz)+1):
            if i_1 == i_2:
                break
            match = CommandManager('cphasematch')
            match.add_command_line_arguments( ['-mtzin', params.output.merged_mtz] )
            match.add_command_line_arguments( ['-mtzout', params.output.phasematch_mtz] )
            match.add_command_line_arguments( ['-colin-fo',   '{}-{},{}-{}'.format(fo1,i_1,fo1,i_1)] )
            match.add_command_line_arguments( ['-colin-fc-1', '{}-{},{}-{}'.format(fc1,i_1,fc1,i_1)] )
            match.add_command_line_arguments( ['-colin-fc-2', '{}-{},{}-{}'.format(fc2,i_2,fc2,i_2)] )
            match.run()
            march.write_output(params.output.phase_log_template+'-{}-{}.log'.format(i_1.i_2))

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
