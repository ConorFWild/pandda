import os

from bamboo.wrappers.ligand_builders import BuilderObject

class ElbowObject(BuilderObject):

    def _custom_init(self):
        self.program = 'phenix.elbow'
        self.name = 'elbow'

    def _create_program_arguments(self, ligsmile, outdir, outfile, flags):
        """Prepare the arguments for Elbow to build the ligand"""

        # Format Settings
        cmd_line_flags, std_inpt_flags = self._convert_flags(flags)
        # Form Command Line Args
        cmd_line_args = ['--smiles',ligsmile,'--output',self.outtemplate,'--name',self.ligname,'--id',self.ligname]
        cmd_line_args.extend(cmd_line_flags)
        # Form Std Input
        std_inpt_args = []
        std_inpt_args.extend(std_inpt_flags)

        return cmd_line_args, std_inpt_args

    def _convert_flags(self, flags):
        """Take the flags passed and generate the command string appropriately"""

        cmd_line_flags = []
        std_inpt_flags = []

        if 'addHydrogens' in flags:
            cmd_line_flags.extend(['--add-hydrogens','True'])
        else:
            cmd_line_flags.extend(['--add-hydrogens','False'])

        return cmd_line_flags, std_inpt_flags
