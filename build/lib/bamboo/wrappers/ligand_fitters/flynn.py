import os, time

from bamboo.wrappers.ligand_fitters import FitterObject

class FlynnObject(FitterObject):

    def _custom_init(self):
        self.program = 'flynn'
        self.name = 'flynn'

    def _custom_prepare_to_fit_ligand(self, ligcif, ligpdb, mtz, apopdb, outdir, outfile, flags):
        """Prepare custom settings before fitting the ligand"""

        # Run the custom parts of this function
        self.outscores = self.outtemplate+'.scores'

    def _create_program_arguments(self, ligcif, ligpdb, mtz, apopdb, outdir, outfile, flags):
        """Create command line and standard input to fit the ligand"""

        # Format Settings
        cmd_line_flags, std_inpt_flags = self._convert_flags(flags)
        # Form Command Line Args
        cmd_line_args = ['-in',ligpdb,'-out',self.outtemplate+'.pdb','-prot',apopdb,'-map',mtz,'-prefix',self.outtemplate+'.auto','-reportfile',self.outscores]
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

        cmd_line_flags = ['-chainid', self.preferChain, '-split']
        std_inpt_flags = []

        return cmd_line_flags, std_inpt_flags

    def _find_output_files(self):
        """Reads the Flynn summary to find the output files"""

        workdir = self.outdir
        summary = self.outscores

        if not os.path.exists(summary):
            raise FittingError('No Fitting Summary File Found: {!s}\n\t\tLog: {!s}'.format(summary, self.outlog))

        summary_contents = open(summary,'r').readlines()

        fields = []
        rows = []
        data = {}

        for i, line in enumerate(summary_contents):
            line=line.strip()
            if line.startswith('Smiles'):
                fields = line.split(',')
            elif fields:
                rows.append(line.split(','))

        if not (fields and rows):
            raise FittingError('No Fitted Files Found in Summary File: {!s}\n\t\tLog: {!s}'.format(summary, self.outlog))

        for i, row in enumerate(rows):
            # Create dictionary of score information
            rowdict = dict(zip(fields,row))
            # Create a rank (the files are already ordered by 'rank' from flynn)
            rank = i+1
            # Create filename
            filename = '{!s}_n{:0>3}_b{:0>3}_s{:0>2}_c{!s:0>3}.pdb'.format(self.outtemplate,str(rank),rowdict['Blobs'],rowdict['Stereo Variant'],int(rowdict['Conformer'])-1)
            # Add a couple more fields
            rowdict['Rank'] = rank
            rowdict['File'] = filename
            rowdict['Path'] = os.path.join(self.outdir,filename)
            # Create a common output score
            rowdict['InternalScore'] = rowdict['RSCC']
            # Check path exists
            assert os.path.exists(rowdict['Path']), 'Flynn Output File does not exist: {!s}'.format(rowdict)
            # Add to list of output model data (ordered by rank from flynn)
            data[rowdict['Path']] = rowdict
            # Add path to list of output models (ordered by rank from flynn)
            self.outmodels.append(rowdict['Path'])

        # Store scoring data
        self.fittingdata = data

        return

