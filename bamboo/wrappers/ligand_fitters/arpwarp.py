import os, time

from bamboo.wrappers.ligand_fitters import FitterObject

class ArpwarpObject(FitterObject):

    def _custom_init(self):
        self.program = 'auto_ligand.sh'
        self.name = 'arpwarp'

    def _create_program_arguments(self, ligcif, ligpdb, mtz, apopdb, outdir, outfile, flags):
        """Run ARPwARP to fit the ligand"""

        # Run the custom parts of this function
        self.outscores = self.outtemplate+'.scores'
# FIXME        (F,P,Rf) = GetMTZRefinementHeadings(mtz)

        # Format Settings
        cmd_line_flags, std_inpt_flags = self._convert_flags(flags)
        # Form Command Line Args
        cmd_line_args = ['datafile',mtz,'protein',apopdb,'ligand',ligpdb,'workdir',outdir,'extralibrary',ligcif,'fp',F,'sigfp',P,'freer',Rf,'ligandfileout',self.outpdb]
        cmd_line_args.extend(cmd_line_flags)
        # Form Standard Input
        std_inpt_args = []
        std_inpt_args.extend(std_inpt_flags)

        return cmd_line_args, std_inpt_args

    def _post_process(self, ligcif, ligpdb, mtz, apopdb, outdir, outfile, flags):

        # Get the output models
        self._find_output_files()

        return

    def _convert_flags(self, flags):
        """Take the flags passed and generate the command string appropriately"""

        cmd_line_flags = []
        std_inpt_flags = []

        for flag in flags:
            pass

        return cmd_line_flags, std_inpt_flags

    def _find_output_files(self):
        """Reads the ARPwARP summary to find the output files"""
        pass

        return

