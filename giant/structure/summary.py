import iotbx.pdb

from scitbx.math import basic_statistics

from giant.xray.edstats import Edstats
from giant.structure.b_factors import BfactorStatistics, normalise_b_factors_to_z_scores
from giant.structure.html import write_html_summary

class StructureSummary(object):

    def __init__(self, pdb_input=None, pdb_hierarchy=None, pdb_file=None, mtz_file=None):
        """Summarise some characteristics of a structure"""

        assert [pdb_file, pdb_input, pdb_hierarchy].count(None)==2,'Provide pdb_file or pdb_input OR pdb_hierarchy'
        if pdb_file:  pdb_input     = iotbx.pdb.input(file_name = pdb_file)
        if pdb_input: pdb_hierarchy = pdb_input.construct_hierarchy()
        # Structure + Data
        self.mtz_file       = mtz_file
        self.pdb_input      = pdb_input
        self.pdb_hierarchy  = pdb_hierarchy
        # B-factors
        self.b_factors = BfactorStatistics.from_pdb(pdb_hierarchy=pdb_hierarchy)
        # Edstats
        if mtz_file: self.edstats = Edstats(mtz_file=mtz_file, pdb_file=pdb_file)
        else:        self.edstats = None

    def normalise_b_factors(self, root=None):
        if root is None: root = self.pdb_hierarchy
        return normalise_b_factors_to_z_scores(pdb_hierarchy=root, b_factor_summary=self.b_factors)

    def _get_atom_group_edstats_scores(self, atom_group):
        resname = atom_group.resname
        chain   = atom_group.parent().parent().id
        res_i   = atom_group.parent().resseq_as_int()
        inscode = ' '
        return self.edstats.scores[(resname, chain, res_i, inscode)]

    def get_atom_group_summary(self, atom_group):
        return atomGroupSummary( atom_group = atom_group,
                                 global_b_factor_statistics = self.b_factors )
    def get_edstats_scores(self, atom_group):
        return self._get_atom_group_edstats_scores(atom_group=atom_group)

    def to_html(self, fname):
        write_html_summary(fname=fname, atom_group_summaries=[self.get_atom_group_summary(ag) for ag in self.pdb_hierarchy.atom_groups()])

class AtomGroupSummary(object):

    def __init__(self, atom_group, global_b_facs_stats=None):

        self.atom_group = atom_group

        # Meta
        self.residue_class = iotbx.pdb.common_residue_names_get_class(atom_group.resname)

        # B-Factors
        self.b_facs         = atom_group.atoms().extract_b()
        self.b_facs_stats   = basic_statistics(atom_group.atoms().extract_b())
        # B-Factors Z-Statistics
        if global_b_facs_stats is not None:
            self.b_facs_z       = global_b_facs_stats.to_z_score(b_vals=atom_group.atoms().extract_b())
            self.b_facs_z_stats = basic_statistics(self.b_factors_z)
        else:
            self.b_facs_z       = None
            self.b_facs_z_stats = None

        # Occupancies
        self.occies         = atom_group.atoms().extract_occ()
        self.occies_stats   = basic_statistics(self.occies)

        # Coordinates
        self.centroid       = atom_group.atoms().extract_xyz().mean()

    def to_html(self, fname):
        write_html_summary(fname=fname, atom_group_summaries=[self])

