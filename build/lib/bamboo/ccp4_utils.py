import os, sys, re, shutil, tempfile

from bamboo.common.command import CommandManager
from bamboo.common.file import FileObj

######################################
#    Misc - move to other module?    #
######################################

def select_factors_and_phases_for_map(cols, maptype='2FOFC'):
    """Takes the Column Headings and returns a dictionary for creating a map type"""

    # Check for Valid Map
    if maptype not in ['2FOFC', 'FOFC']:
        raise Exception('MapType not in allowed map types: {!s} not in {!s}'.format(maptype,MAP_TYPES))
    # Dictionary to contain headings
    dict = {'F':'','PHI':''}
    # Match
    if maptype == '2FOFC':
        if ('2FOFCWT' in cols) and ('PH2FOFCWT' in cols):
            dict['F'] = '2FOFCWT'
            dict['P'] = 'PH2FOFCWT'
        elif ('FWT' in cols) and ('PHWT' in cols):
            dict['F'] = 'FWT'
            dict['P'] = 'PHWT'
        elif ('FWT' in cols) and ('PHFWT' in cols):
            dict['F'] = 'FWT'
            dict['P'] = 'PHFWT'
        else:
            raise Exception('Failed to select STRUCTURE FACTOR Columns for Map Type {!s}. \nColumns : {!s}'.format(maptype,cols))
    elif maptype == 'FOFC':
        if ('FOFCWT' in cols) and ('PHFOFCWT' in cols):
            dict['F'] = 'FOFCWT'
            dict['P'] = 'PHFOFCWT'
        elif ('DELFWT' in cols) and ('DELPHWT' in cols):
            dict['F'] = 'DELFWT'
            dict['P'] = 'DELPHWT'
        elif ('DELFWT' in cols) and ('PHDELWT' in cols):
            dict['F'] = 'DELFWT'
            dict['P'] = 'PHDELWT'
        else:
            raise Exception('Failed to select PHASE Columns: {!s}'.format(cols))
    else:
        raise Exception('MapType does not have functionality programmed: {!s}'.format(maptype))
    # Return
    return dict

######################################
#          PDB functions             #
######################################

def isolate_residue(inpdb, outpdb, resname):
    """Extract the residues identified by resname using pdbcur - i.e. 'UNL'"""

    PDBCUR = CommandManager('pdbcur')
    PDBCUR.add_command_line_arguments('XYZIN',inpdb,'XYZOUT',outpdb)
    PDBCUR.add_standard_input(['lvresidue /*/*/({!s})'.format(resname),'END'])
    PDBCUR.run()

    return PDBCUR

def remove_residue(inpdb, outpdb, resname):
    """Delete the residues identified by resname using pdbcur - i.e. 'UNL'"""

    PDBCUR = CommandManager('pdbcur')
    PDBCUR.add_command_line_arguments('XYZIN',inpdb,'XYZOUT',outpdb)
    PDBCUR.add_standard_input(['delresidue /*/*/({!s})'.format(resname),'END'])
    PDBCUR.run()

    return PDBCUR

def isolate_residue_by_res_id(inpdb, outpdb, chain, resnum, model='*', inscode=''):
    """Isolate the residues identified by residue ids using pdbcur - i.e. '0/A/54.A/'"""

    if inscode:
        selection = '/{!s}/{!s}/{!s}.{!s}'.format(model,chain,resnum,inscode)
    else:
        selection = '/{!s}/{!s}/{!s}'.format(model,chain,resnum)

    PDBCUR = CommandManager('pdbcur')
    PDBCUR.add_command_line_arguments('XYZIN',inpdb,'XYZOUT',outpdb)
    PDBCUR.add_standard_input(['lvresidue {!s}'.format(selection),'END'])
    PDBCUR.run()

    return PDBCUR

def remove_residue_by_res_id(inpdb, outpdb, chain, resnum, model='*', inscode='', removeSolvent=False):
    """Remove the residues identified by res info using pdbcur - i.e. 'UNL'"""

    if inscode:
        selection = '/{!s}/{!s}/{!s}.{!s}'.format(model,chain,resnum,inscode)
    else:
        selection = '/{!s}/{!s}/{!s}'.format(model,chain,resnum)

    std_input = ['delresidue {!s}'.format(selection)]+(removeSolvent)*['delsolvent']+['END']

    PDBCUR = CommandManager('pdbcur')
    PDBCUR.add_command_line_arguments('XYZIN',inpdb,'XYZOUT',outpdb)
    PDBCUR.add_standard_input(std_input)
    PDBCUR.run()

    return PDBCUR

def merge_pdb_files(pdb1, pdb2, pdbout):
    """Merge two PDB Files using CCP4s pdb_merge"""

    # Initialise Commander
    MERGER = CommandManager('pdb_merge')
    # Set command arguments
    MERGER.add_command_line_arguments('xyzin1', pdb1, 'xyzin2', pdb2, 'xyzout', pdbout)
    # Set inputs
    MERGER.add_standard_input('END')
    # run!
    MERGER.run()

    return MERGER

def reset_pdb_file(pdbin, pdbout):
    """Resets the B-Factors in a file and removes anisotropy"""

    if not pdbout.endswith('.pdb'):
        pdbout = pdbout+'.pdb'
    if not os.path.exists(pdbin):
        raise IOError('pdbin does not exist! {!s}'.format(pdbin))
    if os.path.exists(pdbout):
        raise Exception('pdbout already exists! {!s}'.format(pdbout))

    pdbtemp = pdbout.replace('.pdb','.temp.pdb')

    # Initialise Commander
    PDBCUR = CommandManager('pdbcur')
    # Set Command Arguments
    PDBCUR.add_command_line_arguments('XYZIN',pdbin,'XYZOUT',pdbtemp)
    # Set inputs
    PDBCUR.add_standard_input(['NOANISOU','DELSOLVENT','END'])
    # run!
    PDBCUR.run()

    if not os.path.exists(pdbtemp):
        raise ExternalProgramError('PDBCUR has failed to remove anisotropy and delete solvent. {!s}\nOUT: {!s}\nERR: {!s}'.format(pdbin, PDBCUR.out, PDBCUR.err))

    # Initialise Commander
    PDBSET = CommandManager('pdbset')
    # Set Command Arguments
    PDBSET.add_command_line_arguments('XYZIN',pdbtemp,'XYZOUT',pdbout)
    # Set inputs
    PDBSET.add_standard_input(['BFACTOR','END'])
    # run!
    PDBSET.run()

    if not os.path.exists(pdbout):
        raise ExternalProgramError('PDBSET has failed to reset B-factors. {!s}\nOUT: {!s}\nERR: {!s}'.format(pdbtemp, PDBSET.out, PDBSET.err))
    else:
        os.remove(pdbtemp)

    return pdbout

def remove_modres_records(pdb):
    """greps out lines beginning with MODRES from a file"""

    assert os.path.exists(pdb), 'input file must exist'

    origpdb = pdb + '.origfile'
    assert not os.path.exists(origpdb), 'temporary file must not exist'

    shutil.move(pdb, origpdb)

    cmd = "grep -v '^MODRES' {!s} > {!s}".format(origpdb, pdb)
    os.system(cmd)

    if not os.path.exists(pdb):
        raise ExternalProgramError('Failed to remove modres records from {!s}'.format(pdb))

    return

def create_alpha_carbon_backbone(pdbin, pdbout):
    """Takes a pdb files and removes eveything except for the alpha carbons"""

    if not pdbout.endswith('.pdb'):
        pdbout = pdbout+'.pdb'
    if not os.path.exists(pdbin):
        raise IOError('pdbin does not exist! {!s}'.format(pdbin))
    if os.path.exists(pdbout):
        raise Exception('pdbout already exists! {!s}'.format(pdbout))

    # Initialise Commander
    PDBCUR = CommandManager('pdbcur')
    # Set Command Arguments
    PDBCUR.add_command_line_arguments('XYZIN',pdbin,'XYZOUT',pdbout)
    # Set inputs
    PDBCUR.add_standard_input(['lvatom "CA[C]:*"','END'])
    # run!
    PDBCUR.run()

    if not os.path.exists(pdbout):
        raise ExternalProgramError('PDBCUR has failed to create carbon backbone. {!s}\nOUT: {!s}\nERR: {!s}'.format(pdbin, PDBCUR.out, PDBCUR.err))

    return pdbout

def create_cryst_line(pdbin, pdbout, sg, cell):
    """Adds a cryst line to pdbin"""

    if not pdbout.endswith('.pdb'):
        pdbout = pdbout+'.pdb'
    if not os.path.exists(pdbin):
        raise IOError('pdbin does not exist! {!s}'.format(pdbin))
    if os.path.exists(pdbout):
        raise Exception('pdbout already exists! {!s}'.format(pdbout))

    # Initialise Commander
    PDBSET = CommandManager('pdbset')
    # Set Command Arguments
    PDBSET.add_command_line_arguments('XYZIN',os.path.abspath(pdbin),'XYZOUT',os.path.abspath(pdbout))
    # Set Stdin
    PDBSET.add_standard_input(['SPACEGROUP {!s}'.format(sg),'CELL {!s}'.format(' '.join(map(str,cell)))])
    # run!
    PDBSET.run()

    if not os.path.exists(pdbout):
        raise ExternalProgramError('PDBSET has failed to create cryst line for {!s}\nOUT: {!s}\nERR: {!s}'.format(pdbin, PDBSET.out, PDBSET.err))

    return PDBSET

######################################
#              Symmetry              #
######################################

def generate_symmetry_mates(pdbin, pdbout, sgno, cell):
    """Takes the input pdb and generates the unit cell from the point group symmetry"""

    if not pdbout.endswith('.pdb'):
        pdbout = pdbout+'.pdb'
    if not os.path.exists(pdbin):
        raise IOError('pdbin does not exist! {!s}'.format(pdbin))
    if os.path.exists(pdbout):
        raise Exception('pdbout already exists! {!s}'.format(pdbout))

    assert isinstance(sgno, int), 'SPACEGROUP MUST BE AN INTEGER! {!s}'.format(sgno)
    assert isinstance(cell, list), 'CELL MUST BE A LIST! {!s}'.format(cell)

    # Initialise Commander
    PDBSET = CommandManager('pdbset')
    # Set Command Arguments
    PDBSET.add_command_line_arguments('XYZIN',os.path.abspath(pdbin),'XYZOUT',os.path.abspath(pdbout))
    # Set inputs
    PDBSET.add_standard_input(['SYMGEN {!s}'.format(sgno),'CELL {!s}'.format(' '.join(map(str,cell))), 'END'])
    # run!
    PDBSET.run()

    if not os.path.exists(pdbout):
        raise ExternalProgramError('PDBSET has failed to generate SYMMETRY mates. {!s}\nCOM: {!s}\nOUT: {!s}\nERR: {!s}'.format(pdbin,PDBSET.command,PDBSET.out,PDBSET.err))

    return pdbout

def map_to_reference_using_symmetry(refpdb, movpdb, pdbout, conrad=5):
    """Transforms `movpdb` to the closest symmetry site to `refpdb` using csymmatch - Symmetry info must be contained in the header of the pdb file"""

    if not pdbout.endswith('.pdb'):
        pdbout = pdbout+'.pdb'
    if not os.path.exists(refpdb):
        raise IOError('refpdb does not exist! {!s}'.format(refpdb))
    if not os.path.exists(movpdb):
        raise IOError('movpdb does not exist! {!s}'.format(movpdb))
    if os.path.exists(pdbout):
        raise Exception('pdbout already exists! {!s}'.format(pdbout))

    # Initialise Commander
    CSYMMATCH = CommandManager('csymmatch')
    # Set Command Arguments
    CSYMMATCH.add_command_line_arguments('-pdbin-ref',os.path.abspath(refpdb),'-pdbin',os.path.abspath(movpdb),'-pdbout',os.path.abspath(pdbout),'-connectivity-radius',str(conrad))
    # run!
    CSYMMATCH.run()

    if not os.path.exists(pdbout):
        raise ExternalProgramError('CSYMMATCH has failed to map {!s} to {!s}.\nERR: {!s}'.format(movpdb,refpdb,CSYMMATCH.err))

    return pdbout

def generate_full_unit_cell(pdbin, pdbout, mtz=None, sgno=None, cell=None):
    """Map pdbin to the whole of the unit cell using symmetry operations"""

    if not os.path.exists(pdbin):
        raise IOError('pdbin does not exist! {!s}'.format(pdbin))
    if not pdbout.endswith('.pdb'):
        pdbout += '.pdb'
    if os.path.exists(pdbout):
        raise Exception('File already exists: {!s}'.format(pdbout))

    if (sgno != None) and (cell != None):
        pass
    elif (mtz != None):
        # Get the spacegroup number from the mtzfile
        mtzsummary = get_mtz_summary_dict(mtz)
        sgno = mtzsummary['spacegroupno']
        cell = mtzsummary['cell']
    else:
        raise Exception('mtz, or sgno AND cell, must be given!')

    # Generate the symmetry equivalent molecules
    pdbout = generate_symmetry_mates(pdbin=pdbin, pdbout=pdbout, sgno=sgno, cell=cell)

    return pdbout

######################################
#         MTZ/MAP functions          #
######################################

class MtzMeta(object):
    """Object to hold the meta data for an MTZ file"""

    def __init__(self):
        self.reslow = None
        self.reshigh = None
        self.spacegroup = None
        self.spacegroupno = None
        self.cell = None
        self.cryst = None
    def __str__(self):
        return str(self.__dict__)

class ColumnLabels(object):
    """Object to contain the column labels of an MTZ file"""

    def __init__(self):
        self.f = None
        self.sigf = None
        self.phas = None
        self.free = None
        self.comp_f = None
        self.comp_p = None
        self.diff_f = None
        self.diff_p = None
    def __str__(self):
        return str(self.__dict__)

class MtzSummary(object):
    """Class for summarising an MTZ file"""

    def __init__(self, mtz_file):
        self.path = os.path.abspath(mtz_file)
        self.file = FileObj(mtz_file)
        # Read the mtz_file and create summary
        self.summary = get_mtz_summary_dict(mtz_file)
        # Record information from summary
        self.label = ColumnLabels()
        self.label.f    = self.summary['f_labels'][0]    if self.summary['f_labels']    else None
        self.label.sigf = self.summary['sigf_labels'][0] if self.summary['sigf_labels'] else None
        self.label.i    = self.summary['i_labels'][0]    if self.summary['i_labels']    else None
        self.label.sigi = self.summary['sigi_labels'][0] if self.summary['sigi_labels'] else None
        self.label.phas = self.summary['p_labels'][0]    if self.summary['p_labels']    else None
        self.label.free = self.summary['r_flags'][0]     if self.summary['r_flags']     else None
        self.label.comp_f = self.summary['wtmap_f_comp'][0] if self.summary['wtmap_f_comp'] else None
        self.label.comp_p = self.summary['wtmap_p_comp'][0] if self.summary['wtmap_p_comp'] else None
        self.label.diff_f = self.summary['wtmap_f_diff'][0] if self.summary['wtmap_f_diff'] else None
        self.label.diff_p = self.summary['wtmap_p_diff'][0] if self.summary['wtmap_p_diff'] else None
        # Get some Meta
        self.data = MtzMeta()
        self.data.reslow        = self.summary['reslow']        if self.summary['reslow']       else None
        self.data.reshigh       = self.summary['reshigh']       if self.summary['reshigh']      else None
        self.data.spacegroup    = self.summary['spacegroup']    if self.summary['spacegroup']   else None
        self.data.spacegroupno  = self.summary['spacegroupno']  if self.summary['spacegroupno'] else None
        self.data.cell          = self.summary['cell']
        self.data.cryst         = dict(zip(['A','B','C','ALPHA','BETA','GAMMA','SG'],self.data.cell+[self.data.spacegroup]))
        self.data.a,     self.data.b,    self.data.c     = self.summary['cell'][0:3]
        self.data.alpha, self.data.beta, self.data.gamma = self.summary['cell'][3:6]

def get_mtz_summary_dict(mtz_file):
    """Get an MTZ Summary"""

    # Extract Contents of MTZ
    MTZDMP = CommandManager('mtzdmp')
    MTZDMP.add_command_line_arguments(mtz_file)
    MTZDMP.run()
    # Check for errors
    if MTZDMP.process.returncode != 0:
        raise RuntimeError('mtzdmp failed to read file {!s}:\nReturn: {!s}\nOut: {!s}'.format(mtz_file,MTZDMP.process.returncode,MTZDMP.output))

    # Create empty dict to contain the summary
    summary = {}

    # Get the resolution range
    regex = re.compile('\*  Resolution Range :.*\n.*\n.*\((.*)A \)\n')
    matches = regex.findall(MTZDMP.output)
    assert matches, 'No Resolution Range found in MTZFile {!s}'.format(mtz_file)
    assert len(matches)==1, 'Too many matching lines found for Resolution Range in MTZFile {!s}\n\t{!s}'.format(mtz_file,matches)
    summary['reslow'],summary['reshigh'] = map(float,matches[0].replace(' ','').split('-'))

    # Get the Number of Columns
    regex = re.compile('\* Number of Columns =(.*)\n')
    matches = regex.findall(MTZDMP.output)
    assert matches, 'No Number of Columns found for {!s}'.format(mtz_file)
    assert len(matches)==1, 'Too many matching lines found for Number of Columns in MTZFile {!s}\n\t{!s}'.format(mtz_file,matches)
    summary['numcols'] = int(matches[0].strip())

    # Get the Number of Reflections
    regex = re.compile('\* Number of Reflections =(.*)\n')
    matches = regex.findall(MTZDMP.output)
    assert matches, 'No Number of Reflections found for {!s}'.format(mtz_file)
    assert len(matches)==1, 'Too many matching lines found for Number of Reflections in MTZFile {!s}\n\t{!s}'.format(mtz_file,matches)
    summary['numreflections'] = int(matches[0].strip())

    # Get the Column Labels
    regex = re.compile('\* Column Labels :.*\n.*\n(.*)\n')
    matches = regex.findall(MTZDMP.output)
    assert matches, 'No Column Labels found for {!s}'.format(mtz_file)
    assert len(matches)==1, 'Too many matching lines found for Column Headings in MTZFile {!s}\n\t{!s}'.format(mtz_file,matches)
    summary['colheadings'] = matches[0].strip().split()

    # Get the Column Types
    regex = re.compile('\* Column Types :.*\n.*\n(.*)\n')
    matches = regex.findall(MTZDMP.output)
    assert matches, 'No Column Types found for {!s}'.format(mtz_file)
    assert len(matches)==1, 'Too many matching lines found for Column Types in MTZFile {!s}\n\t{!s}'.format(mtz_file,matches)
    summary['coltypes'] = matches[0].strip().split()

    # Get the different datasets
    regex = re.compile('\* Associated datasets :.*\n.*\n(.*)\n')
    matches = regex.findall(MTZDMP.output)
    assert matches, 'No Dataset Numbers found for {!s}'.format(mtz_file)
    assert len(matches)==1, 'Too many matching lines found for Dataset Numbers in MTZFile {!s}\n\t{!s}'.format(mtz_file,matches)
    summary['coldatasets'] = matches[0].strip().split()

    # Get the Spacegroup
    regex = re.compile('\* Space group =.*\'(.*)\'.\(number(.*)\)')
    matches = regex.findall(MTZDMP.output)
    assert matches, 'No Space Group found for {!s}'.format(mtz_file)
    assert len(matches)==1, 'Too many matching lines found for Spacegroup in MTZFile {!s}\n\t{!s}'.format(mtz_file,matches)
    summary['spacegroup'] = matches[0][0].strip()
    summary['spacegroupno'] = int(matches[0][1].strip())

    # Get the Cell Dimensions
    regex = re.compile('\* Cell Dimensions :.*\n.*\n(.*)\n')
    matches = regex.findall(MTZDMP.output)
    assert matches, 'No Cell Dimensions found for {!s}'.format(mtz_file)
    assert len(matches)==1, 'Too many matching lines found for Cell Dimensions in MTZFile {!s}\n\t{!s}'.format(mtz_file,matches)
    summary['cell'] = map(float,matches[0].split())

    # Get the Cell Dimensions

    # Sort the flags
    # F: Amplitudes, H: Miller Indices, I: Integers, P: Phases, Q: Standard Deviations, J: Intensities
    sfac = [summary['colheadings'][i] for i,type in enumerate(summary['coltypes']) if type=='F']
    inty = [summary['colheadings'][i] for i,type in enumerate(summary['coltypes']) if type=='J']
    mill = [summary['colheadings'][i] for i,type in enumerate(summary['coltypes']) if type=='H']
    ints = [summary['colheadings'][i] for i,type in enumerate(summary['coltypes']) if type=='I']
    phas = [summary['colheadings'][i] for i,type in enumerate(summary['coltypes']) if type=='P']
    sdev = [summary['colheadings'][i] for i,type in enumerate(summary['coltypes']) if type=='Q']

    # 2FOFC cols
    wt_f_map_opts = ['2FOFCWT','FWT']
    wt_p_map_opts = ['PH2FOFCWT','PHWT','PHFWT']
    # FOFC cols
    wt_f_dif_opts = ['FOFCWT','DELFWT']
    wt_p_dif_opts = ['PHFOFCWT','DELPHWT','PHDELWT']
    # Record 2FOFC pair (composite map)
    summary['wtmap_f_comp'] = [s for s in sfac if s in wt_f_map_opts]
    summary['wtmap_p_comp'] = [p for p in phas if p in wt_p_map_opts]
    # Record FOFC pair (different map)
    summary['wtmap_f_diff'] = [s for s in sfac if s in wt_f_dif_opts]
    summary['wtmap_p_diff'] = [p for p in phas if p in wt_p_dif_opts]

    # Find the main F col
    summary['f_labels'] = [s for s in sfac if (s in ['F','FP','FCTR','FOSC','FOBS'] \
                            or s.startswith('F_') or s.startswith('FP_') \
                            or s.upper().startswith('F-OBS') or s.upper().startswith('FOUT_'))]
    if not summary['f_labels']:
        summary['f_labels'] = [s for s in sfac if s.startswith('FP')]
    # Find the main SIGF col
    summary['sigf_labels'] = [s for s in sdev if (s in ['SIGF','SIGFP','SIGFOBS'] \
                            or s.startswith('SIGF_') or s.startswith('SIGFP_') \
                            or s.upper().startswith('SIGF-OBS') or s.upper().startswith('SIGFOUT_'))]
    if not summary['sigf_labels']:
        summary['sigf_labels'] = [s for s in sdev if s.startswith('SIGFP')]
    # Find the I cols
    summary['i_labels'] = [s for s in inty if s.startswith('I')]
    # Find the SIGI cols
    summary['sigi_labels'] = [s for s in sdev if s.startswith('SIGI')]
    # Find the F_calc
    summary['f_calcs'] = [s for s in sfac if (s=='FC' or s=='FMODEL' or s.startswith('FC_') or s.upper().startswith('F-MODEL'))]
    # Find the PHI_calcs
    summary['p_calcs'] = [p for p in phas if (p=='PHIC' or p=='PHIFMODEL' or p.startswith('PHIC_') or p.upper().startswith('PHIF-MODEL'))]
    # Find the main phase col
    summary['p_labels'] = summary['p_calcs']
    # Find the RFree Flag
    summary['r_flags'] = [r for r in ints if ('RFREE' in r.upper() or 'FREER' in r.upper() or r=='FREE' or r.upper().startswith('R-FREE'))]

    # XXX Will probably delete this later
    summary['unknown'] = [c for c in summary['colheadings'] if c not in summary['f_labels']+summary['sigf_labels']+summary['f_calcs']+summary['p_calcs']+summary['p_labels']+summary['r_flags']+summary['wtmap_f_comp']+summary['wtmap_p_comp']+summary['wtmap_f_diff']+summary['wtmap_p_diff']]

    return summary

def get_mtz_resolution(mtz_file):
    """Gets the max resolution from the file"""

    # Extract Contents of MTZ
    MTZDMP = CommandManager('mtzdmp')
    MTZDMP.add_command_line_arguments(mtz_file)
    MTZDMP.run()
    # Check for errors
    if MTZDMP.process.returncode != 0:
        raise RuntimeError('mtzdmp failed to read file {!s}:\nReturn: {!s}\nOut: {!s}'.format(mtz_file,MTZDMP.process.returncode,MTZDMP.output))
    # Search for the Column Headings
    regex = re.compile('\*  Resolution Range :.*\n.*\n.*\((.*)A \)\n')
    matches = regex.findall(MTZDMP.output)
    # Check for validity of matches
    assert matches, 'No Resolution Range found in MTZFile {!s}'.format(mtz_file)
    assert len(matches)==1, 'Too many matching lines found for Column Headings in MTZFile {!s}\n\t{!s}'.format(mtz_file,matches)
    # Return
    return map(float,matches[0].replace(' ','').split('-'))

def convert_intensities_to_amplitudes(mtzin, mtzout):
    """Takes an input mtz and converts the intensities to structure factors for model building"""

    mtzobj = MtzSummary(mtzin)

    # Copy data columns from new mtz file
    I, SIGI = mtzobj.label.i, mtzobj.label.sigi
    if not (I and SIGI):
        raise LabelError('No Intensities found in {!s}'.format(mtzin))

    # TODO If rfree is present, retain rfree from the reference mtz
#    RFree = mtzobj.label.free

    # Initialise Commander
    CTRUNC = CommandManager('ctruncate')
    # Set command arguments
    CTRUNC.add_command_line_arguments('-mtzin', mtzin, '-mtzout', mtzout, '-colin', '/*/*/[{!s},{!s}]'.format(I,SIGI))
    # Run!
    CTRUNC.run()

    if not os.path.exists(mtzout):
        raise ExternalProgramError('CTRUNCATE has failed to convert intensities to SFs. {!s}\nOUT: {!s}\nERR: {!s}'.format(mtzin, CTRUNC.out, CTRUNC.err))

    return CTRUNC

def apply_rfree_set(refmtz, mtzin, mtzout):
    """Takes an input mtz and a reference mtz and transplants the rfree flags"""

    tempmtz = mtzin.replace('.mtz','.temp.mtz')

    refobj = MtzSummary(refmtz)
    newobj = MtzSummary(mtzin)

    # Copy data columns from new mtz file
    F1, SIGF1 = newobj.label.f, newobj.label.sigf
    if not (F1 and SIGF1):
        raise LabelError('No Amplitudes found in {!s}'.format(mtzin))

    # Get the spacegroup and the rfree from the reference mtz
    sgno = refobj.data.spacegroupno
    RFree2 = refobj.label.free
    if not sgno:
        raise LabelError('No Spacegroup information found in {!s}'.format(refmtz))
    if not RFree2:
        raise LabelError('No RFreeFlags found in {!s}'.format(refmtz))

    # Initialise Commander
    CAD = CommandManager('cad')
    # Set command arguments
    CAD.add_command_line_arguments('hklin1', mtzin, 'hklin2', refmtz, 'hklout', tempmtz)
    # Set inputs
    CAD.add_standard_input(['symmetry {!s}'.format(sgno),'labin file_number 1 E1={!s} E2={!s}'.format(F1, SIGF1), \
                                               'labout file_number 1 E1={!s} E2={!s}'.format(F1, SIGF1), \
                                               'labin file_number 2 E1={!s}'.format(RFree2), \
                                               'labout file_number 2 E1={!s}'.format(RFree2),'END'])
    # Run!
    CAD.run()

    if not os.path.exists(tempmtz):
        raise ExternalProgramError('CAD has failed to transplant RFree Flags. {!s}\nOUT: {!s}\nERR: {!s}'.format(mtzin, CAD.out, CAD.err))

    # Now complete the free r set

    # Initialise Commander
    FREE = CommandManager('freerflag')
    # Set command arguments
    FREE.add_command_line_arguments('hklin', tempmtz, 'hklout', mtzout)
    # Set inputs
    FREE.add_standard_input(['COMPLETE FREE={!s}'.format(RFree2), 'END'])
    # Run!
    FREE.run()

    if not os.path.exists(mtzout):
        raise ExternalProgramError('freerflag has failed to complete the RFree Flag set. {!s}\nOUT: {!s}\nERR: {s}'.format(tempmtz, CAD.out, CAD.err))

    os.remove(tempmtz)

    return CAD, FREE

def create_cryst_line_from_mtz(pdbin, pdbout, mtz):
    """Takes information from the mtz file and creates a cryst line in pdbin"""

    # Get the spacegroup number from the mtzfile
    mtzsummary = get_mtz_summary_dict(mtz)
    sg   = mtzsummary['spacegroup'].replace(' ','')
    cell = mtzsummary['cell']

    pdbset = create_cryst_line(pdbin, pdbout, sg, cell)

    return pdbset

def fft_mtz_to_map(mtz_file, map_file, cols):
    """Converts an MTZ Format File to a MAP File (using fft as default)"""

    # Initialise
    writer = CommandManager('fft')
    # Set Program Arguments
    writer.add_command_line_arguments('hklin',mtz_file,'mapout',map_file)
    # Set Program Input
    writer.add_standard_input(['LABIN F1={!s} PHI={!s}'.format(cols['F'],cols['P']),'END'])
    # RUN!
    writer.run()
    # Check Output
    if writer.process.returncode != 0:
        print('\nOUT\n\n'+writer.out)
        print('\nERR\n\n'+writer.err)
        raise RuntimeError('fft failed to generate map from {!s}'.format(mtz_file))

    return writer

def mask_map(mapin, maskpdb, mapout, border=1):
    """Takes mapin and masks around atoms in maskpdb"""

    # Masking object
    masker = CommandManager('mapmask')
    # Set input files
    masker.add_command_line_arguments(['mapin',mapin,'mapout',mapout,'xyzin',maskpdb])
    # Set stdin
    masker.add_standard_input(['BORDER {!s}'.format(border),'END'])
    # Run!
    masker.run()
    # Report errors
    if masker.process.returncode!=0:
        raise RuntimeError('mapmask failed to mask map {!s}'.format(mapin))

    # Return Command Managers for flexible handling of out & err
    return masker

def create_asu_map(mapin, mapout):
    """Takes mapin and masks to the asymmetric unit"""

    # Masking object
    masker = CommandManager('mapmask')
    # Set input files
    masker.add_command_line_arguments(['mapin',mapin,'mapout',mapout])
    # Set stdin
    masker.add_standard_input(['XYZLIM ASU','END'])
    # Run!
    masker.run()
    # Report errors
    if masker.process.returncode!=0:
        raise RuntimeError('mapmask failed to create asu map from {!s}'.format(mapin))

    # Return Command Managers for flexible handling of out & err
    return masker

def create_masked_map(mtzin, maskpdb, mapout=None, maptype='2FOFC', border=2):
    """Creates density from mtzin, and masks around maskpdb - NOTE THAT MASKPDB MUST BE ALIGNED WITH MTZIN"""

    # Create mapfile if not given
    if not mapout:
        mapout = os.path.splitext(maskpdb)[0]+'.masked.'+maptype+'.ccp4'
    # Temporary (unmasked) map
    tempmap = mapout+'.nomask'
    # Convert the inital MAT to a MAP file
    fft = convert_mtz_to_map(mtzfile=mtzin, mapfile=tempmap, maptype=maptype)
    # Mask the map
    masker = mask_map(mapin=tempmap,maskpdb=maskpdb,mapout=mapout,border=border)
    # Remove the temporary mapfile
    os.remove(tempmap)
    # Return Command Managers for flexible handling of out & err
    return fft, masker

def convert_mtz_to_map(mtzfile, mapfile=None, method='fft', maptype='2FOFC', force=False):
    """Converts an MTZ Format File to a MAP File (using fft as default)"""

    # Create mapfile if not given
    if not mapfile:
        mapfile = os.path.splitext(mtzfile)[0]+'.'+maptype+'.ccp4'
    # Check for file validity
    if os.path.exists(mapfile):
        if force:
            os.remove(mapfile)
        else:
            return None

    # Read the Column Headings from the MTZ File
    cols = get_mtz_summary_dict(mtzfile)['colheadings']
    selcols = select_factors_and_phases_for_map(cols,maptype)

    if method == 'fft':
        writer = fft_mtz_to_map(mtzfile, mapfile, selcols)
    else:
        raise Exception('Program not recognised: {!s}'.format(method))

    return writer

######################################
#          CIF functions             #
######################################

def merge_cif_libraries(incifs, outcif):
    """Take a list of cifs and merge into one cif"""

    assert isinstance(incifs, list), "'incifs' is not a list!"
    assert len(incifs) > 1, "'incifs' must be two or more files!"

    current = incifs.pop(0)

    to_be_deleted = []

    for additional in incifs:

        # Create a file handle and path for the output
        temp_handle, temp_path = tempfile.mkstemp(suffix='.lib', prefix='libcheck_')
        to_be_deleted.append(temp_path)

        # Merge current and additional to temp_path
        LIBCHK = CommandManager('libcheck')
        LIBCHK.add_standard_input('_DOC N','_FILE_L {!s}'.format(current),'_FILE_L2 {!s}'.format(additional),'_FILE_O {!s}'.format(temp_path.replace('.lib','')),'_END')
        LIBCHK.run()

        current = temp_path

    shutil.copy(current, outcif)

    for file in to_be_deleted:
        os.remove(file)

    assert os.path.exists(outcif), 'OUTPUT CIF DOES NOT EXIST! {!s}'.format(outcif)

    return outcif

######################################
#              Bfactors              #
######################################

def get_pdb_bfactor_summary(pdb_file):
    """Analyse the b-factors in a pdb file and return a dictionary with all of the averages, rms deviations and Z-scores for the b-factors in a structure"""

    program, contents = extract_bfactor_statistics(pdb_file)
    # Process the output table from baverage
    processed_contents = process_baverage_output(contents)
    # Process the stdout from the program
    processed_stdout = process_baverage_stdout(program.out)

    return program, processed_contents

def process_baverage_output(contents):
    """Process baverage output"""

    chain = ''
    resid = ''

    processed = {}

    for line in contents:
        if not line:
            continue
        elif line.strip().startswith('#'):
            chain = line.split()[2]
            assert len(chain)==1, 'Chain ID is more than one letter!'
        else:
            if not chain:
                raise Exception('No Chain ID has been assigned!')

            resid, m_av, m_rms, s_av, s_rms, a_av, a_rms = map(float, line.split())
            resid = int(resid)

            processed[(chain, resid)] = bFactorSummary(chain, resid, m_av, m_rms, s_av, s_rms, a_av, a_rms)

    return processed

def extract_bfactor_statistics(pdb_file):
    """Analyse the b-factors in a pdb file and return a dictionary with all of the averages, rms deviations and Z-scores for the b-factors in a structure"""

    # Create a temporary file for the output summary table and pdb_file
    temp_handle, temp_path = tempfile.mkstemp(suffix='.table', prefix='baverage_')
    temp_handle2, temp_path2 = tempfile.mkstemp(suffix='.pdb', prefix='baverage_')

    BAVERAGE = CommandManager('baverage')
    BAVERAGE.add_command_line_arguments('XYZIN',pdb_file,'RMSTAB',temp_path,'XYZOUT',temp_path2)
    BAVERAGE.add_standard_input(['END'])
    BAVERAGE.run()

    if not os.path.exists(temp_path):
        raise ExternalProgramError('BAVERAGE has failed to calculate b-factor summary statistics for {!s}'.format(pdb_file))

    # Process Table and Standard Out
    table_contents = open(temp_path, 'r').read().split('\n')

    if not table_contents:
        raise ExternalProgramError('BAVERAGE has failed to calculate b-factor summary statistics for {!s}'.format(pdb_file))
    else:
        os.remove(temp_path)
        os.remove(temp_path2)

    return BAVERAGE, table_contents

