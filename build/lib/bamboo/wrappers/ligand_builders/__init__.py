import os

from bamboo.common.command import CommandManager
from bamboo.constants import DEFAULT_LIGAND_NAMES
from bamboo.wrappers import allowed_builder_args

class BuilderObject(object):
    """Template Object for creating an Object to generate ligand restraints"""

    def __init__(self, time=True, verbose=True):

        self.allowed_args = allowed_builder_args
        self.ligname = DEFAULT_LIGAND_NAMES[0]
        # Settings
        self.time = time
        self.verbose = verbose
        self.runtime = -1.0
        self.timeout = 1800
        # Custom Init
        self._custom_init()

    def generate_ligand(self, ligsmile, outdir, outfile, flags=[]):
        """Generate ligand PDB and CIF from smile string"""

        # Prepare standard settings
        self._prepare_to_generate_ligand(ligsmile, outdir, outfile, flags)
        # Check to see if already run, and if not, run the custom part of the object - different for each program
        if os.path.exists(self.outpdb) and os.path.exists(self.outcif) and os.path.exists(self.outlog):
            if self.verbose:
                print('\tLigand already generated - DOING NOTHING')
        else:
            cmd_line_args, std_inpt_args = self._create_program_arguments(ligsmile, outdir, outfile, flags)
            # Initialise CommandManager
            self.builder = CommandManager(self.program)
            # Set Command-Line Args
            if cmd_line_args: self.builder.add_command_line_arguments(cmd_line_args)
            # Set Standard Input
            if std_inpt_args: self.builder.add_standard_input(std_inpt_args)
            # Set Parameters
            self.builder.set_timeout(timeout=self.timeout)
            # RUN
            self.builder.run()
            # Get runtime
            self.runtime = self.builder.runtime

            try:
                # Could add postprocess here
                pass

            finally:
                self.write_log_file()

        return self.outpdb, self.outcif, self.outlog

    def _prepare_to_generate_ligand(self, ligsmile, outdir, outfile, flags):
        """Set up the generic file names"""

        # Process outputfile
        if '/' in outfile:
            raise ValueError('outfile must be a file, not a path')
        # Record Template and Outdir
        self.outtemplate = os.path.join(outdir,outfile)
        self.outdir      = outdir
        # Record Filenames
        self.outpdb = self.outtemplate+'.pdb'
        self.outcif = self.outtemplate+'.cif'
        self.outlog = self.outtemplate+'.log'

    def write_log_file(self):
        """Write the log file"""

        with open(self.outlog,'w') as logfile:
            # Write out the input command
            logfile.write('\nCOMMAND\n\n')
            logfile.write('\n'.join(self.builder.cmd_line_args)+'\n')
            logfile.write('\nINPUT\n\n')
            logfile.write('\n'.join(self.builder.std_inp_lines)+'\n')
            # Write out & err
            logfile.write('\nSTDOUT\n\n')
            logfile.write(self.builder.output)
            logfile.write('\nSTDERR\n\n')
            logfile.write(self.builder.error)

