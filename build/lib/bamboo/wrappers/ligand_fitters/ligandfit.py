import os, time, re

from bamboo.wrappers.ligand_fitters import FitterObject

class LigandfitObject(FitterObject):

    def _custom_init(self):
        self.program = 'phenix.ligandfit'
        self.name = 'ligandfit'

    def _custom_prepare_to_fit_ligand(self, ligcif, ligpdb, mtz, apopdb, outdir, outfile, flags):
        """Prepare custom settings before fitting the ligand"""

        # Set the ligandfit output directory
        self.workdir = os.path.join(outdir,'LigandFit_run_1_')
        # Check for existence of other LigandFit directories
        run_directory = 1
        # Increment output number until we get a new one (This will be where phenix puts the output)
        while os.path.exists(self.workdir):
            run_directory += 1
            self.workdir = os.path.join(outdir,'LigandFit_run_{!s}_'.format(run_directory))

        # Check the penultimate one. If there's a fitting summary, take that one instead (pipeline already run before)
        if (run_directory > 1) and os.path.exists(self.outpdb) and os.path.exists(self.outlog):
            test_workdir = os.path.join(outdir,'LigandFit_run_{!s}_'.format(run_directory-1))
            test_outscores = os.path.join(test_workdir,'LigandFit_summary.dat')

            if os.path.exists(test_outscores):
                self.workdir = test_workdir

        # Where the scores and output files from the program will be recorded
        self.outscores = os.path.join(self.workdir,'LigandFit_summary.dat')

        # Change directory because of the 'unique' way that phenix writes its files
        os.chdir(outdir)

    def _create_program_arguments(self, ligcif, ligpdb, mtz, apopdb, outdir, outfile, flags):
        """Create command line and standard input to fit the ligand"""

        # Format Settings
        cmd_line_flags, std_inpt_flags = self._convert_flags(flags)
        # Form Command Line Args
        cmd_line_args = ['model='+apopdb,'data='+mtz,'ligand='+ligpdb,'cif_def_file_list='+ligcif]
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

        cmd_line_flags.append('number_of_ligands=5')

        return cmd_line_flags, std_inpt_flags

    def _find_output_files(self):
        """Reads the LigandFit summary to find the output files"""

        workdir = self.workdir
        summary = self.outscores

        if not os.path.exists(summary):
            raise FittingError('No Fitting Summary File Found: {!s}\n\t\tLog: {!s}'.format(summary, self.outlog))

        summary_contents = open(summary,'r').readlines()

        fields = ['Path','RSCC','Ligandfit Score','Num Atoms']
        rows = []
        data = {}

        for i, line in enumerate(summary_contents):
            if line.startswith('SINGLE LIGAND MODELS'):
                # Extract lists
                files = eval(line.replace('SINGLE LIGAND MODELS: ',''))
                newline = summary_contents[i+1].strip()
                assert newline.startswith('SINGLE LIGAND CC VALUES: '), 'LIGANDFIT SCORE FILE NOT IN EXPECTED FORMAT! {!s}'.format(summary)
                rsccs = eval(newline.replace('SINGLE LIGAND CC VALUES: ',''))
                newline = summary_contents[i+2].strip()
                assert newline.startswith('SINGLE LIGAND SCORES   : '), 'LIGANDFIT SCORE FILE NOT IN EXPECTED FORMAT! {!s}'.format(summary)
                ligsc = eval(newline.replace('SINGLE LIGAND SCORES   : ',''))
                newline = summary_contents[i+3].strip()
                assert newline.startswith('SINGLE LIGAND PLACED   : '), 'LIGANDFIT SCORE FILE NOT IN EXPECTED FORMAT! {!s}'.format(summary)
                noats = eval(newline.replace('SINGLE LIGAND PLACED   : ',''))
                # zip into table
                rows = zip(files,rsccs,ligsc,noats)

        if not (fields and rows):
            raise FittingError('No Fitted Files Found in Summary File: {!s}\n\t\tLog: {!s}'.format(summary, self.outlog))

        # Order by Ligandfit score
        for i, row in enumerate(sorted(rows,key=lambda entr: float(entr[2]),reverse=True)):
            # Create a dictionary of the row information
            rowdict = dict(zip(fields,row))
            # Create a rank (the files are already ordered by a rhofit score)
            rank = i+1
            # Add a couple more fields to the dict
            rowdict['Rank'] = rank
            # Create a common output score
            rowdict['InternalScore'] = rowdict['Ligandfit Score']
            # Check path exists
            assert os.path.exists(rowdict['Path']), 'Ligandfit Output File does not exist: {!s}'.format(rowdict)

            # Add to list of output model data (ordered by rank from ligandfit)
            data[rowdict['Path']] = rowdict
            # Add path to list of output models (ordered by rank from ligandfit)
            self.outmodels.append(rowdict['Path'])

        # Store scoring data
        self.fittingdata = data

        return

