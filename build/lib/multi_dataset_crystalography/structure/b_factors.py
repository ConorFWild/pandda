import iotbx.pdb

from scitbx.array_family import flex
from scitbx.math import basic_statistics

from giant.maths.geometry import is_within
from giant.structure.select import non_h, protein, backbone, sidechains
from giant.structure.formatting import ShortLabeller

########################################################################################

class BfactorStatistics(object):
    """B-factor statstics for a structure"""
    def __init__(self, all=None, protein=None, backbone=None, sidechain=None):
        self.all       = all
        self.protein   = protein
        self.backbone  = backbone
        self.sidechain = sidechain

    def to_z_score(self, b_vals, method='backbone'):
        """Convert a set of B-factors to Z-scores"""
        assert method in ['all','protein','backbone','sidechain']
        if method == 'all':       stats = self.all
        if method == 'protein':   stats = self.protein
        if method == 'backbone':  stats = self.backbone
        if method == 'sidechain': stats = self.sidechain
        if stats.biased_standard_deviation == 0.0:
            # No variation - return list of zeroes
            return b_vals-b_vals

        return (b_vals - stats.mean)/stats.biased_standard_deviation

    @classmethod
    def from_pdb(cls, pdb_input=None, pdb_hierarchy=None):
        """Calculate the b-factor statistics of a model"""

        assert [pdb_input, pdb_hierarchy].count(None)==1,'Provide pdb_input OR pdb_hierarchy'
        if pdb_input: pdb_hierarchy = pdb_input.construct_hierarchy()

        cache = pdb_hierarchy.atom_selection_cache()

        all_b       = non_h(hierarchy=pdb_hierarchy, cache=cache, copy=True).atoms().extract_b()
        protein_b   = protein(hierarchy=pdb_hierarchy, cache=cache, copy=True).atoms().extract_b()
        backbone_b  = backbone(hierarchy=pdb_hierarchy, cache=cache, copy=True).atoms().extract_b()
        sidechain_b = sidechains(hierarchy=pdb_hierarchy, cache=cache, copy=True).atoms().extract_b()

        return cls(all       = basic_statistics(all_b),
                   protein   = basic_statistics(protein_b),
                   backbone  = basic_statistics(backbone_b),
                   sidechain = basic_statistics(sidechain_b) )

    def show(self):
        if self.all:
            print '================>'
            print 'All Atoms:'
            print self.all.show()
        if self.protein:
            print '================>'
            print 'Protein Atoms:'
            print self.protein.show()
        if self.backbone:
            print '================>'
            print 'Backbone Atoms:'
            print self.backbone.show()
        if self.sidechain:
            print '================>'
            print 'Sidechain Atoms:'
            print self.sidechain.show()


########################################################################################

def normalise_b_factors_to_z_scores(pdb_input=None, pdb_hierarchy=None, b_factor_statistics=None, method='backbone'):
    """Calculate the b-factor statistics of the model"""

    assert [pdb_input, pdb_hierarchy].count(None)==1,'Provide pdb_input OR pdb_hierarchy'
    if pdb_input: pdb_hierarchy = pdb_input.construct_hierarchy()

    if not b_factor_statistics: b_factor_statistics = BfactorStatistics.from_pdb(pdb_hierarchy=pdb_hierarchy)
    new_b = b_factor_statistics.to_z_score(b_vals=pdb_hierarchy.atoms().extract_b(), method=method)

    output_h = pdb_hierarchy.deep_copy()
    output_h.atoms().set_b(new_b)

    return output_h

def occupancy_weighted_average_b_factor(atoms):
    return flex.mean_weighted(atoms.extract_b(), atoms.extract_occ())

def calculate_residue_group_bfactor_ratio(residue_group, hierarchy, data_table=None, rg_label=None, column_suffix=''):
    """Calculate bfactor quality metrics of the residue to surrounding residues"""

    rg_sel = residue_group
    # Set defaults
    if rg_label is None:    rg_label = (rg_sel.unique_resnames()[0]+'-'+rg_sel.parent().id+'-'+rg_sel.resseq+rg_sel.icode).replace(' ','')
    if data_table is None:  data_table = pandas.DataFrame(index=[rg_label], column=[])
    # Check validity
    if len(rg_sel.unique_resnames()) != 1: raise Exception(rg_label+': More than one residue name associated with residue group -- cannot process')

    # Extract rg_sel objects
    rg_sel_ags    = rg_sel.atom_groups()
    rg_sel_atoms  = rg_sel.atoms()
    rg_sel_coords = rg_sel_atoms.extract_xyz()
    # Select nearby atom_group for scoring
    near_ags = [ag for ag in hierarchy.atom_groups() if (is_within(4, rg_sel_coords, ag.atoms().extract_xyz()) and (ag not in rg_sel_ags))]
    if near_ags:
        # Extract the names of the nearby residues
        near_ag_names = ':'.join(sorted(set([ShortLabeller.format(ag) for ag in near_ags])))
        data_table.set_value(   index = rg_label,
                                col   = 'Surrounding Residues Names'+column_suffix,
                                value = near_ag_names )
        # Extract atoms from nearby groups
        near_ats = iotbx.pdb.hierarchy.af_shared_atom()
        [near_ats.extend(ag.detached_copy().atoms()) for ag in near_ags]
        # Calculate B-factors of the residue
        res_mean_b = occupancy_weighted_average_b_factor(atoms=rg_sel_atoms)
        data_table.set_value(   index = rg_label,
                                col   = 'Average B-factor (Residue)'+column_suffix,
                                value = res_mean_b )
        # Calculate B-factors of the surrounding atoms
        sch_mean_b = occupancy_weighted_average_b_factor(atoms=near_ats)
        data_table.set_value(   index = rg_label,
                                col   = 'Average B-factor (Surroundings)'+column_suffix,
                                value = sch_mean_b )
        # Store the ratio of the b-factors
        data_table.set_value(   index = rg_label,
                                col   = 'Surroundings B-factor Ratio'+column_suffix,
                                value = res_mean_b/sch_mean_b )

    return data_table

########################################################################################


