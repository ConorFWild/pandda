import os

from bamboo.wrappers.ligand_builders import BuilderObject

class WritedictObject(BuilderObject):

    def _custom_init(self):
        self.program = 'writedict'
        self.name = 'writedict'

    def _create_program_arguments(self, ligsmile, outdir, outfile, flags):
        """Prepare the arguments for Writedict to build the ligand"""

        # Create input smiles file
        self.infile = self.outtemplate+'.smi'
        # Write smile file
        with open(self.infile, 'w') as fh:
            fh.write(ligsmile)

        # Format Settings
        cmd_line_flags, std_inpt_flags = self._convert_flags(flags)
        # Form Command Line Args
        cmd_line_args = ['-in',self.infile,'-out',self.outtemplate]
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
            raise SystemExit('AAAHHH - FLAG NOT PROGRAMMED FOR WRITEDICT - `noHydrogens`')
        if 'noHydrogens' in flags:
            raise SystemExit('AAAHHH - FLAG NOT PROGRAMMED FOR WRITEDICT - `noHydrogens`')
        # Change Ligand Name
        if 'correctNames' in flags:
            pass
        else:
            cmd_line_flags.append('-nolookup')

        return cmd_line_flags, std_inpt_flags
