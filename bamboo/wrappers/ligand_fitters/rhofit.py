import os, time

from bamboo.wrappers.ligand_fitters import FitterObject

class RhofitObject(FitterObject):

    def _custom_init(self):
        self.program = 'rhofit'
        self.name = 'rhofit'

    def _custom_prepare_to_fit_ligand(self, ligcif, ligpdb, mtz, apopdb, outdir, outfile, flags):
        """Prepare custom settings before fitting the ligand"""

        # Run the custom parts of this function
        self.workdir = os.path.join(outdir,'rhofit')

    def _create_program_arguments(self, ligcif, ligpdb, mtz, apopdb, outdir, outfile, flags):
        """Create command line and standard input to fit the ligand"""

        # Format Settings
        cmd_line_flags, std_inpt_flags = self._convert_flags(flags)
        # Form Command Line Args
        cmd_line_args = ['-l', ligcif, '-lp', ligpdb, '-m', mtz, '-p', apopdb, '-d', self.workdir]
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

        if '-keepHydrogens' in flags:
            cmd_line_flags.append('-keepH')

        return cmd_line_flags, std_inpt_flags

    def _find_output_files(self):
        """Reads the Rhofit summary to find the output files"""

        workdir = self.workdir
        summary = os.path.join(self.workdir, 'results.txt')

        if not os.path.exists(summary):
            raise FittingError('No Fitting Summary File Found: {!s}\n\t\tLog: {!s}'.format(summary, self.outlog))

        summary_contents = open(summary,'r').readlines()

        fields = []
        rows = []
        data = {}

        for i, line in enumerate(summary_contents):
            # Remove newlines and whitespace
            line = line.strip()
            if line.startswith('========================'):
                assert summary_contents[i-3].strip() == 'rhofit           ligand LigProt  Poorly', 'RHOFIT SCORE FILE NOT IN EXPECTED FORMAT! {!s}'.format(summary)
                assert summary_contents[i-2].strip() == 'total   Correl  strain contact fitting   LigProt contact to residues', 'RHOFIT SCORE FILE NOT IN EXPECTED FORMAT! {!s}'.format(summary)
                assert summary_contents[i-1].strip()   == 'File               Chain    score   coeff    score   score   atoms   (% means zero weighted in score)', 'RHOFIT SCORE FILE NOT IN EXPECTED FORMAT! {!s}'.format(summary)
                fields = ['File','Chain','Rhofit Score','RSCC','Strain Score','Contact Score','Poor Atoms','Contacts']
            elif fields and line:
                values = line.split()
                if len(values)>6:
                    # Add non-blank lines to the lines to process
                    rows.append(line.split())

        if not (fields and rows):
            raise FittingError('No Fitted Files Found in Summary File: {!s}\n\t\tLog: {!s}'.format(summary, self.outlog))

        # Order by rhofit score
        for i, row in enumerate(sorted(rows, key=lambda entr: float(entr[2]))):
            # Sort out the last column (splits into many columns)
            proc_row = row[0:7]
            proc_row.append(';'.join(row[7:]))
            # Create dictionary of score information
            rowdict = dict(zip(fields,proc_row))
            # Create a rank (the files are already ordered by a rhofit score)
            rank = i+1
            # Get the cluster, dictionary number, and rank within the cluster from the filename
            cluster, dictionary_num, cluster_rank = proc_row[0].replace('.pdb','').replace('Hit_','').split('_')
            # Add a couple more fields to the dict
            rowdict['Rank'] = rank
            rowdict['Cluster'] = int(cluster)
            rowdict['Cluster Rank'] = int(cluster_rank)
            rowdict['Dictionary Number'] = int(dictionary_num)
            # Copy the File to Path for Neatness
            rowdict['Path'] = os.path.join(workdir,rowdict['File'])
            # Create a common output score
            rowdict['InternalScore'] = rowdict['Rhofit Score']
            # Check path exists
            assert os.path.exists(rowdict['Path']), 'Rhofit output File does not exist: {!s}'.format(rowdict)

            # Add to list of output model data (ordered by rank from rhofit)
            data[rowdict['Path']] = rowdict
            # Add path to list of output models (ordered by rank from rhofit)
            self.outmodels.append(rowdict['Path'])

        # Store scoring data
        self.fittingdata = data

        return
