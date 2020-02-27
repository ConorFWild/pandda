import os, shutil

from bamboo.common.command import CommandManager
from bamboo.constants import DEFAULT_OUTPUT_CHAIN
from bamboo.wrappers import allowed_fitter_args
from bamboo.ccp4_utils import map_to_reference_using_symmetry, merge_pdb_files, change_residue_chain_and_number, create_cryst_line_from_mtz

class FitterObject(object):

    def __init__(self, time=True, verbose=True):

        self.allowedArgs = allowed_fitter_args
        self.preferChain = DEFAULT_OUTPUT_CHAIN[0]
        # Settings
        self.time = time
        self.verbose = verbose
        self.runtime = -1.0
        self.timeout = 99999
        # Custom Init
        self._custom_init()

    def fit_ligand(self, ligcif, ligpdb, mtz, apopdb, outdir, outfile, flags=[]):
        """Fit Ligand to a structure"""

        # Prepare input and output filenames
        self._prepare_to_fit_ligand(ligcif, ligpdb, mtz, apopdb, outdir, outfile, flags)

        # Check to see if already run, and if not, run the custom part of the object - different for each program
        if os.path.exists(self.outpdb) and os.path.exists(self.outlog):
            if self.verbose:
                print('\tLigand already fitted - DOING NOTHING')
            self._post_process(ligcif, ligpdb, mtz, apopdb, outdir, outfile, flags)
        else:
            if self.verbose:
                print('\tFitting using {!s}.'.format(self.name))
            cmd_line_args, std_inpt_args = self._create_program_arguments(ligcif, ligpdb, mtz, apopdb, outdir, outfile, flags)
            # Initialise CommandManager
            self.Fitter = CommandManager(self.program)
            # Set Command-Line Args
            if cmd_line_args: self.Fitter.SetArguments(cmd_line_args)
            # Set Standard Input
            if std_inpt_args: self.Fitter.SetInput(std_inpt_args)
            # Set Parameters
            self.Fitter.SetParameters(timeout=self.timeout)
            # RUN
            self.Fitter.Run()
            # Calculate runtime (seconds)
            self.runtime = self.Fitter.runtime

            try:
                # Process results and generate the list of modelled ligands
                if self.verbose:
                    print('\tFinished fitting - post-processing models.')
                self._post_process(ligcif, ligpdb, mtz, apopdb, outdir, outfile, flags)

                if not self.outmodels:
                    raise FittingError('No fitted models found! Log: {!s}'.format(self.outlog))

                if not os.path.exists(self.outpdb):
                    # Map the model of the best ligand to the asu and check it's near to the protein
                    self._map_to_asu(ligcif, ligpdb, mtz, apopdb, outdir, outfile, flags)
                    # Create symbolic links to the output file
                    self._link_output_files()

                # Merge the best ligand with the apo structure
                if not os.path.exists(self.mergedpdb):
                    self.Merger = merge_pdb_files(apopdb, self.outpdb, self.mergedpdb)

            finally:
                self.write_log_file()

        return self.outpdb, self.mergedpdb, self.outlog

    def _prepare_to_fit_ligand(self, ligcif, ligpdb, mtz, apopdb, outdir, outfile, flags):
        """Set up the generic file names"""

        # Process outputfile
        if '/' in outfile:
            raise ValueError('outfile must be a file, not a path')
        # Record Template and Outdir
        self.outtemplate = os.path.join(outdir,outfile)
        self.outdir      = outdir
        # Record Filenames
        self.outpdb      = self.outtemplate+'.lig.pdb'
        self.mergedpdb   = self.outtemplate+'.pdb'
        self.outlog      = self.outtemplate+'.log'
        # Storage
        self.outmodels   = []
        self.filtmodels  = []
        self.bestmodel   = None
        self.fittingdata = None

        # Prepare custom settings
        self._custom_prepare_to_fit_ligand(ligcif, ligpdb, mtz, apopdb, outdir, outfile, flags)

        return

    def _map_to_asu(self, ligcif, ligpdb, mtz, apopdb, outdir, outfile, flags):
        """Map the best ligand file to the asu of the protein - pick the symmetry version with the most contacts"""

        candidate = None

        if self.verbose:
            print('\tMapping fitted models back to the asymmetric unit.')

        for model in self.outmodels:

            ligand_sym_file = model.replace('.pdb','.sym.pdb')

            if os.path.exists(ligand_sym_file):
                self.filtmodels.append(model)
            else:
                try:
                    # Create Cryst Line
                    shutil.move(model, model+'.temp.pdb')
                    PDBSET = create_cryst_line_from_mtz(model+'.temp.pdb', ligand_sym_file, mtz)
                    # Rename the original file, then create a symmetry equivalent ligand that is close to the protein
#                    shutil.move(model, ligand_sym_file)
#                    candidate, ordered_files = map_structure_to_asu(apopdb, mtz, ligand_sym_file, model, delete_structures=True)
                    candidate = map_to_reference_using_symmetry(refpdb=apopdb, movpdb=ligand_sym_file, pdbout=model)
                    self.filtmodels.append(candidate)
                except (DistanceError, IOError) as err:
                    print(err)
                    continue

            # If it gets here we have a candidate
            # Break here for fast (1 model per fitter max)
#            break

        if not self.filtmodels:
            raise FittingError('No fitted models within range of protein!')
        else:
            # Rename all of the filtered models
            [change_residue_chain_and_number(filtfile) for filtfile in self.filtmodels]
            # Pick the `best` one to be the output
            self.bestmodel = self.filtmodels[0]

        return

    def _link_output_files(self):
        """Make links to the output files"""

        # Create Symlinks
        os.symlink(os.path.relpath(self.bestmodel,start=self.outdir),self.outpdb)

        return

    def write_log_file(self):
        """Write the log file"""

        if self.verbose:
            print('\tWriting fitting logfile.')

        with open(self.outlog,'w') as logfile:

            if hasattr(self, 'Fitter'):
                # Write out the input command
                logfile.write('\nCOMMAND\n\n')
                logfile.write('\n'.join(self.Fitter.command))
                logfile.write('\nINPUT\n\n')
                logfile.write(self.Fitter.inp)
                # Write out & err
                logfile.write('\nFITTER STDOUT\n\n')
                logfile.write(self.Fitter.out)
                logfile.write('\nFITTER STDERR\n\n')
                logfile.write(self.Fitter.err)

            if hasattr(self, 'Merger'):
                # Write out the input command
                logfile.write('\nCOMMAND\n\n')
                logfile.write('\n'.join(self.Merger.command))
                logfile.write('\nINPUT\n\n')
                logfile.write(self.Merger.inp)
                # Write out & err
                logfile.write('\nMERGER STDOUT\n\n')
                logfile.write(self.Merger.out)
                logfile.write('\nMERGER STDERR\n\n')
                logfile.write(self.Merger.err)

        return
