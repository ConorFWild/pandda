import os

from bamboo.common.command import CommandManager
from bamboo.wrappers import allowed_refiner_args
from bamboo.ccp4_utils import remove_modres_records

class RefinerObject(object):

    def __init__(self, time=True, verbose=True):

        self.allowedArgs = allowed_refiner_args
        # Settings
        self.time = time
        self.verbose = verbose
        self.runtime = -1.0
        self.timeout = 1800
        # Custom Init
        self._custom_init()

    def run_refinement(self, inpdb, inmtz, outdir, outfile, incif=None, flags=[]):
        """Run Refinement of a Structure"""

        # Prepare standard settings
        self._prepare_to_refine_structure(inpdb, inmtz, outdir, outfile, incif, flags)
        # Check to see if already run, and if not, run the custom part of the object - different for each program
        if os.path.exists(self.outpdb) and os.path.exists(self.outmtz) and os.path.exists(self.outlog):
            if self.verbose:
                print('\tRefinement already done - DOING NOTHING')
        else:
            if self.verbose:
                print('\tRefining using {!s}.'.format(self.name))
            cmd_line_args, std_inpt_args = self._create_program_arguments(inpdb, inmtz, outdir, outfile, incif, flags)
            # Initialise CommandManager
            self.Refiner = CommandManager(self.program)
            # Set Command-Line Args
            self.Refiner.SetArguments(cmd_line_args)
            # Set Standard Input
            self.Refiner.SetInput(std_inpt_args)
            # Set Parameters
            self.Refiner.SetParameters(timeout=self.timeout)
            # RUN
            self.Refiner.Run()
            # Calculate runtime (seconds)
            self.runtime = self.Refiner.runtime

            if self.verbose:
                print('\tFinished refinement - post-processing structure.')

            try:
                # General post-process
                self._standard_post_process(inpdb, inmtz, outdir, outfile, incif, flags)

                # Could add autoprocessing here...
                self._post_process(inpdb, inmtz, outdir, outfile, incif, flags)

            finally:
                self.write_log_file()

        return self.outpdb, self.outmtz, self.outlog

    def _prepare_to_refine_structure(self, inpdb, inmtz, outdir, outfile, incif, flags):
        """Set up the generic file names"""

        # Process outputfile
        if '/' in outfile:
            raise ValueError('outfile must be a file, not a path')
        # Record Template and Outdir
        self.outtemplate = os.path.join(outdir,outfile)
        self.outdir = outdir
        # Record Filenames
        self.outpdb = self.outtemplate+'.pdb'
        self.outmtz = self.outtemplate+'.mtz'
        self.outcif = self.outtemplate+'.cif'
        self.outlog = self.outtemplate+'.log'

    def _standard_post_process(self, inpdb, inmtz, outdir, outfile, incif, flags):
        """Perform standard operations on the output files of the refinement"""

        if os.path.exists(self.outpdb):
            remove_modres_records(self.outpdb)

        return

    def write_log_file(self):
        """Write the log file"""

        if self.verbose:
            print('\tWriting refinement logfile.')

        with open(self.outlog,'w') as logfile:
            # Write out the input command
            logfile.write('\nCOMMAND\n\n')
            logfile.write('\n'.join(self.Refiner.command))
            logfile.write('\nINPUT\n\n')
            logfile.write(self.Refiner.inp)
            # Write out & err
            logfile.write('\nSTDOUT\n\n')
            logfile.write(self.Refiner.out)
            logfile.write('\nSTDERR\n\n')
            logfile.write(self.Refiner.err)
