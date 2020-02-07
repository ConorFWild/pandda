import os, time

from bamboo.wrappers.refiners import RefinerObject

class RefmacObject(RefinerObject):

    def _custom_init(self):
        self.program = 'refmac5'
        self.name = 'refmac'

    def _create_program_arguments(self, inpdb, inmtz, outdir, outfile, incif, flags):
        """Run Refmac to refine the structure"""

        # Format Settings
        cmd_line_flags, std_inpt_flags = self._convert_flags(flags)
        # Form Command Line Args
        cmd_line_args = ['XYZIN',inpdb,'HKLIN',inmtz,'XYZOUT',self.outpdb,'HKLOUT',self.outmtz]
        # Form Std Input
        std_inpt_args = []
        # Add cif files
        if incif and isinstance(incif, str):
            cmd_line_args.extend(['LIBIN',incif,'LIBOUT',self.outcif])
        elif incif and isinstance(incif, list):
            cmd_line_args.extend(['LIBIN',incif[-1],'LIBOUT',self.outcif])
        # Add processed flags
        cmd_line_args.extend(cmd_line_flags)
        std_inpt_args.extend(std_inpt_flags)
        std_inpt_args.append('END')

        return cmd_line_args, std_inpt_args

    def _post_process(self, inpdb, inmtz, outdir, outfile, incif, flags):
        pass

    def _convert_flags(self, flags):
        """Take the flags passed and generate the command string appropriately"""

        cmd_line_flags = []
        std_inpt_flags = []

        if 'keepH' in flags:
            flags.remove('keepH')
            std_inpt_flags.append('HOUT YES')
        # Check for refine b-factors
        if 'bonly' in flags:
            flags.remove('bonly')
            std_inpt_flags.append('REFI BONLY')
        # Check for noexit command for new ligands
        if 'noexit' in flags:
            flags.remove('noexit')
            std_inpt_flags.append('MAKE NEWLigand Noexit')
        # Check for occupancy refinement
        occu_residues = [f[1] for f in flags if f[0] == 'occupancy']
        # Add occupancy commands
        if occu_residues:
            # Clear up flags
            [flags.remove(f) for f in flags if f[1] in occu_residues]
            # Generate Command list for occupancy refinement
            group = 1
            for label in occu_residues:
                chain = label[0]
                resid = label[1]
                inp = 'OCCU GROUP ID {!s} CHAIN {!s} RESIDUE {!s}'.format(group, chain, resid)
                std_inpt_flags.append(inp)
                group += 1
            std_inpt_flags.append('OCCU REFINE')

        return cmd_line_flags, std_inpt_flags

