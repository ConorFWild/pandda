import os

from bamboo.wrappers.ligand_builders import BuilderObject

class ProdrgObject(BuilderObject):

    def _custom_init(self):
        self.program = 'cprodrg'
        self.name = 'prodrg'

    def _create_program_arguments(self, ligsmile, outdir, outfile, flags):
        """Prepare the arguments for Prodrg to build the ligand"""

        # Create input smiles file
        self.infile = self.outtemplate+'.smile'
        self.outmol = self.outtemplate+'.mol'
        # Write smile file
        with open(self.infile, 'w') as file:
            file.write(ligsmile+'\n')
            file.write('CPNAME '+self.ligname)

        # Format Settings
        cmd_line_flags, std_inpt_flags = self._convert_flags(flags)
        # Form Command Line Args
        cmd_line_args = ['XYZIN',self.infile,'LIBOUT',self.outcif,'XYZOUT',self.outpdb,'MOLOUT',self.outmol]
        cmd_line_args.extend(cmd_line_flags)
        # Form Std Input
        std_inpt_args = ['MINI YES']
        std_inpt_args.extend(std_inpt_flags)
        std_inpt_args.append('END')

        return cmd_line_args, std_inpt_args

    def _convert_flags(self, flags):
        """Take the flags passed and generate the command string appropriately"""

        cmd_line_flags = []
        std_inpt_flags = []

        if 'noHydrogens' in flags:
            raise SystemExit('AAAHHH - FLAG NOT PROGRAMMED FOR PRODRG - `noHydrogens`')

        return cmd_line_flags, std_inpt_flags
