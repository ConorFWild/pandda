# Lists/Constants

# =================== #
# Refiner Settings    #
# =================== #

refiners = ['refmac','phenix.refine']
allowed_refiner_args = []

# =================== #
# Ligand Generation   #
# =================== #

ligand_builders = ['elbow','grade','writedict']
#ligand_builders = ['elbow','grade','prodrg','writedict']
allowed_builder_args = []

# =================== #
# Ligand Fitting      #
# =================== #

ligand_fitters = ['ligandfit','rhofit','flynn']
#ligand_fitters = ['phenix.ligandfit','arpwarp','rhofit','coot','jligand','flynn','afitt']
allowed_fitter_args = []

# =================== #

def get_ligand_builder(program):
    """Get Ligand Builder Object for `program`"""

    from bamboo.wrappers.ligand_builders.elbow import ElbowObject
    from bamboo.wrappers.ligand_builders.grade import GradeObject
    from bamboo.wrappers.ligand_builders.prodrg import ProdrgObject
    from bamboo.wrappers.ligand_builders.writedict import WritedictObject

    # Check it's a valid program
    if program not in ligand_builders:      raise LigandBuilderSelectionError('{!s} not in {!s}'.format(program, ', '.join(ligand_builders)))

    if program in ['elbow','phenix.elbow']:     return ElbowObject()
    elif program == 'grade':                    return GradeObject()
    elif program == 'prodrg':                   return ProdrgObject()
    elif program in ['writedict','afitt']:      return WritedictObject()
    else:
        raise LigandBuilderSelectionError('MY BAD. The Code for {!s} has not been written yet...'.format(program))

def get_ligand_fitter(program):
    """Get Ligand Fitter Object for `program`"""

    from bamboo.wrappers.ligand_fitters.ligandfit import LigandfitObject
    from bamboo.wrappers.ligand_fitters.arpwarp import ArpwarpObject
    from bamboo.wrappers.ligand_fitters.rhofit import RhofitObject
    from bamboo.wrappers.ligand_fitters.flynn import FlynnObject
    from bamboo.wrappers.ligand_fitters.coot import CootObject

    # Check it's a valid program
    if program not in ligand_fitters:       raise LigandFitterSelectionError('{!s} not in {!s}'.format(program, ', '.join(ligand_fitters)))

    if program in ['ligandfit','phenix.ligandfit']:     return LigandfitObject()
    elif program == 'arpwarp':                          return ArpwarpObject()
    elif program == 'rhofit':                           return RhofitObject()
    elif program in ['flynn','afitt']:                  return FlynnObject()
    elif program == 'coot':                             return CootObject()
    else:
        raise LigandFitterSelectionError('MY BAD. The Code for {!s} has not been written yet...'.format(program))

def get_refiner(program):
    '''Returns Initialised Refiner Object'''

    from bamboo.wrappers.refiners.phenix import PhenixrefineObject
    from bamboo.wrappers.refiners.refmac import RefmacObject

    # Check it's a valid program
    if program not in refiners:
        raise RefinerSelectionError('{!s} not in {!s}'.format(program, ', '.join(refiners)))

    if program.lower()=='phenix.refine':                return PhenixrefineObject()
    elif program.lower()=='refmac':                     return RefmacObject()
    else:
        raise RefinerSelectionError('MY BAD. The Code for {!s} has not been written yet...'.format(program))

