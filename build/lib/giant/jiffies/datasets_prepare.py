#!/usr/bin/env ccp4-python

import os, sys, copy, re, shutil

import libtbx.phil
import libtbx.easy_mp

from libtbx.utils import Sorry, Failure
from bamboo.common.path import splice_ext, foldername, filename, easy_directory
from bamboo.common.command import CommandManager
from giant.xray.crystal import CrystalSummary

############################################################################

PROGRAM = 'giant.datasets.prepare'

DESCRIPTION = """
    Take input MTZs, reindex, transfer R-free flags from a reference, fill missing reflections, and optionally run refinement.
"""

############################################################################

blank_arg_prepend = {'.mtz':'mtz=', '.pdb':'reference_pdb='}

master_phil = libtbx.phil.parse("""
input {
    mtz = None
        .type = path
        .multiple = True
    reference_pdb = None
        .type = path
    reference_mtz = None
        .type = path
}
output {
    output_directory = None
        .type = str
    labelling = none *filename foldername
        .type = choice
    keep_intermediate_files = False
        .type = bool
}
actions = *reindex_to_reference \
          *fill_missing_reflections \
          *transfer_rfree_flags \
          *refinement_pipeline
    .type = choice(multi=True)
options {
    reindex {
        tolerance = 5
            .type  = float
    }
    rfree_flag {
        input = FreeR_flag
           .type = str
        output = None
            .type = str
    }
    missing_reflections {
        fill_high = 4.0
            .type = float
        fill_low = 9999.0
            .type = float
    }
    refinement {
        pipeline = *dimple
            .type = choice(multi=True)
        col_label = IMEAN
            .type = str
    }
}
settings {
    cpus = 1
        .help = "number of cpus to use for processing"
        .type = int
}
""")

############################################################################

def raise_cmd_output_and_error(cmd):
    """Return STDOUT and STDERR from command object"""

    err_msg = ''
    err_msg += '============================>\n'
    err_msg += 'Program returned with an error ({})\n'.format(' '.join(cmd.program))
    err_msg += '============================>\n'
    err_msg += str(cmd) + '\n'
    err_msg += '============================>\n'
    err_msg += cmd.output + '\n'
    err_msg += '============================>\n'
    err_msg += cmd.error  + '\n'
    err_msg += '============================>\n'
    raise Failure(err_msg)

def print_run_and_raise_error_maybe(cmd):
    """Print settings for CMD object and raise error if fails"""
    cmd.print_settings()
    ret_code = cmd.run()
    if ret_code != 0:
        raise_cmd_output_and_error(cmd)

############################################################################

def reindex_mtz_to_reference(in_mtz, out_mtz, reference_mtz, tolerance):
    """Reindex the data in one mtz to a reference mtz"""

    print '**************************'
    print '*** Running reindexing ***'
    print '**************************'

    cmd = CommandManager('pointless')
    cmd.add_command_line_arguments(['hklin',  in_mtz,
                                    'hklref', reference_mtz,
                                    'hklout', out_mtz   ])
    cmd.add_standard_input([ 'tolerance {}'.format(tolerance) ])
    print_run_and_raise_error_maybe(cmd)
    if not os.path.exists(out_mtz):
        raise Failure('reindexing has failed -- {} does not exist'.format(out_mtz))

def fill_missing_reflections(in_mtz, out_mtz, fill_resolution_low, fill_resolution_high, delete_tmp_files=True):
    """Complete the set of miller indices in an MTZ file"""

    print '***************************'
    print '*** Filling reflections ***'
    print '***************************'

    tmp_mtz_1 = splice_ext(path=out_mtz, new='step1-truncate')
    tmp_mtz_2 = splice_ext(path=out_mtz, new='step2-uniquify')
    tmp_mtz_3 = splice_ext(path=out_mtz, new='step3-remerged')

    # Stage 1 - truncate dataset, fill missing reflections, change column name
    cmd = CommandManager('cad')
    cmd.add_command_line_arguments([ 'hklin1', in_mtz,
                                     'hklout', tmp_mtz_1 ])
    cmd.add_standard_input([ 'monitor BRIEF',
                             'labin file_number 1 ALL',
                             'resolution file 1 {} {}'.format(fill_resolution_low, fill_resolution_high) ])
    print_run_and_raise_error_maybe(cmd)
    if not os.path.exists(tmp_mtz_1):
        raise Failure('filling of missing reflections has failed -- {} does not exist'.format(tmp_mtz_1))

    print '-------------------'

    # Stage 2 - Uniqueify the file
    cmd = CommandManager('uniqueify')
    cmd.add_command_line_arguments([ '-p', '0.05',
                                     tmp_mtz_1,
                                     tmp_mtz_2 ])
    print_run_and_raise_error_maybe(cmd)
    if not os.path.exists(tmp_mtz_2):
        raise Failure('filling of missing reflections has failed -- {} does not exist'.format(tmp_mtz_2))

    print '-------------------'

    # Stage 3 - remerge the two files
    cmd = CommandManager('cad')
    cmd.add_command_line_arguments([ 'hklin1', in_mtz,
                                     'hklin2', tmp_mtz_2,
                                     'hklout', tmp_mtz_3 ])
    cmd.add_standard_input([ 'monitor BRIEF',
                             'labin file_number 1 ALL',
                             'labin file_number 2 E1=FreeR_flag',
                             'labout file_number 2 E1=dummy' ])
    print_run_and_raise_error_maybe(cmd)
    if not os.path.exists(tmp_mtz_3):
        raise Failure('filling of missing reflections has failed -- {} does not exist'.format(tmp_mtz_3))

    print '-------------------'

    # Stage 4 - remove the dummy column
    cmd = CommandManager('mtzutils')
    cmd.add_command_line_arguments([ 'hklin1', tmp_mtz_3,
                                     'hklout', out_mtz ])
    cmd.add_standard_input([ 'HEADER BRIEF',
                             'EXCLUDE 1 dummy',
                             'ONEFILE',
                             'END' ])
    print_run_and_raise_error_maybe(cmd)
    if not os.path.exists(tmp_mtz_3):
        raise Failure('filling of missing reflections has failed -- {} does not exist'.format(out_mtz))

    if delete_tmp_files:
        os.remove(tmp_mtz_1)
        os.remove(tmp_mtz_2)
        os.remove(tmp_mtz_3)

def transfer_rfree_flags(in_mtz, out_mtz, reference_mtz, input_free_r_flag, output_free_r_flag, delete_tmp_files=True):
    """Copy R-free flags from reference mtz"""

    print '*********************************'
    print '*** Transferring R-free flags ***'
    print '*********************************'

    tmp_mtz = splice_ext(path=out_mtz, new='step1-transfer')

    # Stage 1 - transfer R-free from reference
    cmd = CommandManager('cad')
    cmd.add_command_line_arguments(['hklin1', in_mtz,
                                    'hklin2', reference_mtz,
                                    'hklout', tmp_mtz   ])
    cmd.add_standard_input(['labin file_number 1 ALL',
                            'labin file_number 2 E1={}'.format(input_free_r_flag),
                            'labout file_number 2 E1={}'.format(output_free_r_flag),
                            'END'     ])
    print_run_and_raise_error_maybe(cmd)
    if not os.path.exists(tmp_mtz):
        raise Failure('transfer of R-free flags has failed -- {} does not exist'.format(tmp_mtz))

    print '*******************************'
    print '*** Completing R-free flags ***'
    print '*******************************'

    # Stage 2 - populate missing R-free values
    cmd = CommandManager('freerflag')
    cmd.add_command_line_arguments(['hklin',  tmp_mtz,
                                    'hklout', out_mtz   ])
    cmd.add_standard_input(['COMPLETE FREE={}'.format(output_free_r_flag),
                            'END'])
    print_run_and_raise_error_maybe(cmd)
    if not os.path.exists(out_mtz):
        raise Failure('expanding of R-free flags has failed -- {} does not exist'.format(out_mtz))

    if delete_tmp_files:
        os.remove(tmp_mtz)

def run_refinement_pipeline(in_mtz, ref_pdb, out_dir, program='dimple'):
    """Run refinement of the input MTZ file against a reference PDB file"""

    print '*************************'
    print '*** Running Pipelines ***'
    print '*************************'

    if program == 'dimple':
        # Define output files
        out_pdb = os.path.join(out_dir, 'final.pdb')
        out_mtz = os.path.join(out_dir, 'final.mtz')
        # Create command manager for dimple
        cmd = CommandManager('dimple')
        cmd.add_command_line_arguments([ '--jelly', '5',
                                         in_mtz,
                                         ref_pdb,
                                         out_dir ])
    else:
        raise Exception("no stop that. you're doing it wrong.")

    print_run_and_raise_error_maybe(cmd)
    if not os.path.exists(out_pdb):
        raise Failure('running refinement with {} has failed -- {} does not exist'.format(program, out_pdb))
    if not os.path.exists(out_mtz):
        raise Failure('running refinement with {} has failed -- {} does not exist'.format(program, out_mtz))

    return out_pdb, out_mtz

############################################################################

def reprocess_mtz(in_mtz, params):
    """Prepare mtz by various standard sanitisations"""

    # Determine naming for this dataset
    if params.output.labelling == 'foldername':
        label = foldername(in_mtz)
    elif params.output.labelling == 'filename':
        label = filename(in_mtz)
    else:
        label = None

    if label is not None:
        print '=======================================>>>'
        print 'Processing dataset {}'.format(label)
        print '=======================================>>>'

    # Create output directory or identify it
    if params.output.output_directory is not None:
        out_dir = easy_directory(os.path.join(params.output.output_directory, label))
        # Copy input mtz to output directory and update variable
        new_in_mtz = os.path.join(out_dir, os.path.basename(in_mtz))
        shutil.copy(in_mtz, new_in_mtz)
        # Update variable to point to new file
        in_mtz = new_in_mtz
    else:
        # Simply take current directory -- likely not used
        out_dir = None

    # Save the original mtz as variable
    orig_mtz = in_mtz

    stage = 1

    ###########################################
    #           Reindex to reference          #
    ###########################################

    if 'reindex_to_reference' in params.actions:

        out_mtz = splice_ext(path=orig_mtz, new='{}.reindexed'.format(stage))
        reindex_mtz_to_reference(in_mtz  = in_mtz,
                                 out_mtz = out_mtz,
                                 reference_mtz = params.input.reference_mtz,
                                 tolerance = params.options.reindex.tolerance)
        if (in_mtz != orig_mtz) and not params.output.keep_intermediate_files:
            os.remove(in_mtz)
        in_mtz = out_mtz
        stage += 1

    ###########################################
    #         Fill missing reflections        #
    ###########################################

    if 'fill_missing_reflections' in params.actions:

        out_mtz = splice_ext(path=orig_mtz, new='{}.filled'.format(stage))
        fill_missing_reflections(in_mtz  = in_mtz,
                                 out_mtz = out_mtz,
                                 fill_resolution_low  = params.options.missing_reflections.fill_low,
                                 fill_resolution_high = params.options.missing_reflections.fill_high,
                                 delete_tmp_files = (not params.output.keep_intermediate_files))
        if (in_mtz != orig_mtz) and not params.output.keep_intermediate_files:
            os.remove(in_mtz)
        in_mtz = out_mtz
        stage += 1

    ###########################################
    #      Transfer+complete R-free flags     #
    ###########################################

    if 'transfer_rfree_flags' in params.actions:

        out_mtz = splice_ext(path=orig_mtz, new='{}.with-free'.format(stage))
        transfer_rfree_flags(in_mtz  = in_mtz,
                             out_mtz = out_mtz,
                             reference_mtz = params.input.reference_mtz,
                             input_free_r_flag  = params.options.rfree_flag.input,
                             output_free_r_flag = params.options.rfree_flag.output,
                             delete_tmp_files = (not params.output.keep_intermediate_files))
        if (in_mtz != orig_mtz) and not params.output.keep_intermediate_files:
            os.remove(in_mtz)
        in_mtz = out_mtz
        stage += 1

    ################################################
    # Copy the current MTZ to "prepared" extension #
    ################################################

    out_mtz = splice_ext(path=orig_mtz, new='prepared')
    if params.output.keep_intermediate_files:
        shutil.copy(in_mtz, out_mtz)
    else:
        shutil.move(in_mtz, out_mtz)

    in_mtz = out_mtz

    ################################################
    #           Run refinements pipelines          #
    ################################################

    if 'refinement_pipeline' in params.actions:

        for pipeline in params.options.refinement.pipeline:

            if out_dir is None:
                ref_dir = os.path.splitext(in_mtz)[0] + '_' + pipeline
            else:
                ref_dir = os.path.join(out_dir, pipeline)

            out_pdb, out_mtz = run_refinement_pipeline(in_mtz  = in_mtz,
                                                       ref_pdb = params.input.reference_pdb,
                                                       out_dir = ref_dir,
                                                       program = pipeline)

    return out_mtz

############################################################################

success_message = 'processed successfully'

def wrapper_reprocess(args):
    params, mtz = args
    try:
        reprocess_mtz(in_mtz=mtz, params=params)
    except Exception as e:
        return (mtz, e)
    return (mtz, success_message)

############################################################################

def run(params):

    assert params.input.reference_mtz is not None, 'no reference mtz supplied (provided with reference_mtz=XXX)'
    assert len(params.input.mtz) > 0, 'no mtz files supplied'

    if params.output.output_directory:
        params.output.output_directory = easy_directory(params.output.output_directory)
        assert params.output.labelling != 'none'

    if 'reindex_to_reference' in params.actions:
        assert params.input.reference_mtz is not None, 'no reference mtz supplied (provided with reference_mtz=XXX)'

    if 'fill_missing_reflections' in params.actions:
        assert params.options.missing_reflections.fill_high < params.options.missing_reflections.fill_low

    if 'transfer_rfree_flags' in params.actions:
        assert params.input.reference_mtz is not None, 'no reference mtz supplied (provided with reference_mtz=XXX)'
        assert params.options.rfree_flag.input is not None, 'R-free flag must be specified (rfree_flag.input=XXX)'
        if params.options.rfree_flag.output is None:
            params.options.rfree_flag.output = params.options.rfree_flag.input

    if 'refinement_pipeline' in params.actions:
        assert params.input.reference_pdb is not None, 'must provide reference_pdb to run refinement_pipeline'

    #################

    arg_list = [(params, m) for m in params.input.mtz]
    print '============================>'
    print 'Running re-processing for {} mtz files'.format(len(arg_list))
    returned = libtbx.easy_mp.pool_map(fixed_func   = wrapper_reprocess,
                                       args         = arg_list,
                                       processes    = params.settings.cpus,
                                       chunksize    = 1)

    print '============================>'
    print 'Processed Successfully:'
    success = 0
    errors = []
    for mtz, message in returned:
        if message == success_message:
            print '\t{}'.format(mtz)
            success += 1
        else:
            errors.append((ret))

    if errors:
        print '============================>'
        print 'Errors:'
        print '=============>'
        for mtz, err in errors:
            print 'File: {}\n\t{}:{}'.format(mtz, type(message), str(message))
            print '============================>'

    print ''
    print '============================>'
    print 'Successfully processed:     {}'.format(success)
    print 'Errors during reprocessing: {}'.format(len(errors))

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
