import os

from bamboo.wrappers.ligand_builders import BuilderObject

class GradeObject(BuilderObject):

    def _custom_init(self):
        self.program = 'grade'
        self.name = 'grade'

    def _create_program_arguments(self, ligsmile, outdir, outfile, flags):
        """Prepare the arguments for Grade to build the ligand"""

        # Format Settings
        cmd_line_flags, std_inpt_flags = self._convert_flags(flags)
        # Form Command Line Args
        cmd_line_args = [ligsmile,'-name',self.outtemplate,'-opdb',self.outpdb,'-ocif',self.outcif,'-f','-resname',self.ligname]
        cmd_line_args.extend(cmd_line_flags)
        # Form Standard Input Args
        std_inpt_args = []
        cmd_line_args.extend(std_inpt_flags)

        return cmd_line_args, std_inpt_args

    def _convert_flags(self, flags):
        """Take the flags passed and generate the command string appropriately"""

        cmd_line_flags = []
        std_inpt_flags = []

        if 'noHydrogens' in flags:
            cmd_line_flags.append('-really_noH')
        if 'mogul' not in flags:
            cmd_line_flags.append('-nomogul')

        return cmd_line_flags, std_inpt_flags
