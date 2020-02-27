
# ==================== #
# PDB Format Constants #
# ==================== #

PDB_HEADERS_TO_KEEP = ('HEADER','TITLE','REMARK','CRYST','SCALE')
PDB_FOOTERS_TO_KEEP = ('CONECT',)

# =================== #
# Ligand Constants    #
# =================== #

DEFAULT_LIGAND_NAMES = ('UNL',)
DEFAULT_OUTPUT_CHAIN = ('X',)
DEFAULT_OUTPUT_RESNUM = (666,)

# =================== #
# Peptide Constants   #
# =================== #

LONGEST_PEPTIDE = 20
LONGEST_OLIGO = 6

# =================== #
# MTZ Column Headings #
# =================== #

MAP_TYPES = ('2FOFC','FOFC')

# =================== #
# Compound Constants  #
# =================== #

COMMON_COMPOUNDS = ('[Zn]','[Ni]','K','[K]','[KH]','Mg','[Mg]','[MgH2]','[Na]','[NaH]','I','Cl','NCS','CC(C)O','CS(C)O','OCCO','CC(O)O','OS(O)(O)O','O','O=CC(O)CO','OCC(O)CO','OP(O)(O)O','OCO')

# =================== #
# Ion Constants       #
# =================== #

ION_NAMES = ('FE','BR','CL','CA','CO','CU','ZN','MG','MN','CD','F','NA','HG','PB','IN','NI','IOD','SR','YB','AL','IUM','K','LI','RB','FE2','BA','SM','CS','MN3','CU1','TL','PT','CR','AG','GA')

WATER_NAMES = ('HOH','H2O')

SOLVENT_NAMES = ('ACT', 'EDO')

# =================== #
# Protein Constants   #
# =================== #

AA_3_TO_1_CODES_DICT = {'ALA':'A','GLY':'G','LEU':'L','SER':'S','VAL':'V','THR':'T','LYS':'K','ASP':'D','ILE':'I','ASN':'N','GLU':'E','PRO':'P','ARG':'R','PHE':'F','GLN':'Q','TYR':'Y','HIS':'H','CYS':'C','MET':'M','TRP':'W'}

# Epigenetic modifications
AA_MODIFICATIONS_DICT = {  'CC(NCCCC[C@@H](C=O)N)=O'               : ('K(ac)', 'K(Ac)', '{KAc}', 'Kac'),
                           'CC(NCCCC[C@@H](C(=O)O)N)=O'            : ('K(ac)', 'K(Ac)', '{KAc}', 'Kac'),
                           'CNCCCC[C@@H](C=O)N'                    : ('K(me)'),
                           'CN(C)CCCC[C@@H](C=O)N'                 : ('K(me2)', 'K(2me)'),
                           'C[N+](C)(C)CCCC[C@@H](C=O)N'           : ('K(me3)', 'K(3me)'),
                           'CN(C)C(=N)NCCC[C@@H](C=O)N'            : ('R(ame2)', 'R(Asym2me)', 'R(ame)'),
                           'CNC(=N)NCCC[C@@H](C=O)N'               : ('R(me)'),
                           'CNC(NCCC[C@@H](C=O)N)=NC'              : ('R(me2)'),
                           'CC([C@@H](C=O)N)OP(O)(O)=O'            : ('T(p)'),
                           'C([C@@H](C=O)N)OP(O)(O)=O'             : ('S(p)'),
                           'C(c1ccc(cc1)OP(O)(O)=O)[C@@H](C=O)N'   : ('Y(p)')
                           }

MISC_SMILES_DICT = {'C[C@@H](OC1CC(C2CNN(C)C2)CNC1N)[C@@H]1C(Cl)CCC(F)C1Cl'   :   ('HEME')}

ATOM_NUMBER_DICT = {1:'H',2:'He',3:'Li',4:'Be',5:'B',6:'C',7:'N',8:'O',9:'F',10:'Ne',11:'Na',12:'Mg',13:'Al',14:'Si',15:'P',16:'S',17:'Cl',18:'Ar',19:'K',20:'Ca',21:'Sc',22:'Ti',23:'V',24:'Cr',25:'Mn',26:'Fe',27:'Co',28:'Ni',29:'Cu',30:'Zn',31:'Ga',32:'Ge',33:'As',34:'Se',35:'Br',36:'Kr',37:'Rb',38:'Sr',39:'Y',40:'Zr',41:'Nb',42:'Mo',43:'Tc',44:'Ru',45:'Rh',46:'Pd',47:'Ag',48:'Cd',49:'In',50:'Sn',51:'Sb',52:'Te',53:'I',54:'Xe',55:'Cs',56:'Ba',57:'La',58:'Ce',59:'Pr',60:'Nd',61:'Pm',62:'Sm',63:'Eu',64:'Gd',65:'Tb',66:'Dy',67:'Ho',68:'Er',69:'Tm',70:'Yb',71:'Lu',72:'Hf',73:'Ta',74:'W',75:'Re',76:'Os',77:'Ir',78:'Pt',79:'Au',80:'Hg',81:'Tl',82:'Pb',83:'Bi',84:'Po',85:'At',86:'Rn',87:'Fr',88:'Ra',89:'Ac',90:'Th',91:'Pa',92:'U',93:'Np',94:'Pu',95:'Am',96:'Cm',97:'Bk',98:'Cf',99:'Es',100:'Fm',101:'Md',102:'No',103:'Lr',104:'Rf',105:'Db',106:'Sg',107:'Bh',108:'Hs',109:'Mt',110:'Ds',111:'Rg',112:'Cn',113:'Uut',114:'Fl',115:'Uup',116:'Lv',117:'Uus',118:'Uuo'}

# ==================== #
# Nucleotide Constants #
# ==================== #

NUCLEOTIDE_DICT = {'DC':'C','DA':'A','DG':'G','DT':'T','Cd':'C','Ad':'A','Gd':'G','Td':'T','U':'U','G':'G','A':'A','C':'C'}


