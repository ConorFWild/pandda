from __future__ import print_function

import os, sys, glob, time, re

#################################
try:
    import matplotlib
    matplotlib.use('Agg')
    matplotlib.interactive(0)
    from matplotlib import pyplot
    pyplot.style.use('ggplot')
except:
    pass
#################################

import numpy, pandas

import iotbx.pdb

from scitbx.array_family import flex

from bamboo.common import Meta, Info
from bamboo.plot import simple_histogram, simple_bar, simple_scatter

from bamboo.maths.angles import mean_std_angles

from giant.structure import make_label
from giant.structure.dihedrals import get_all_phi_psi_for_hierarchy
from giant.structure.b_factors import BfactorStatistics, normalise_b_factors_to_z_scores
from giant.structure.iterators import residues_via_conformers, conformers_via_residue_groups

import giant.structure.select as s_select

class StructureCollection(object):
    """Class to hold information about a set of related structures"""

    def __init__(self, hierarchy_inputs=None, hierarchies=None, inputs=None, labels=None):

        # Object to stores all of the structure objects
        self.structures = Info(['hierarchies','inputs','labels'])
        # Unpack and store structure objects
        if hierarchy_inputs:
            self.structures.hierarchies = [i.hierarchy for i in hierarchy_inputs]
            self.structures.inputs      = [i.input     for i in hierarchy_inputs]
        else:
            self.structures.hierarchies = hierarchies
            self.structures.inputs      = inputs
        # Unpack and store labels
        if not labels: labels = range(len(self.structures.hierarchies))
        self.structures.labels = labels

        #Initialise output tables
        self.initialise_tables()

    @classmethod
    def from_files(cls, filenames, labels=None):
        """Initialise the object from files"""
        if not labels: labels=filenames
        return cls(hierarchy_inputs=[iotbx.pdb.hierarchy.input(f) for f in filenames], labels=labels)

    ######################################################################################

    def initialise_tables(self):
        """Populate the index of the tables"""
        # Object to contain all of the output data tables
        self.tables = Info(['structures','residues','atoms'])
        self._initialise_structure_table()
        self._initialise_residue_table()

    def _initialise_structure_table(self):
        """Initialise the tables.structures object with structure labels"""
        self.tables.structures = pandas.DataFrame( data = None,
                                                   index   = pandas.Index(data=self.structures.labels, name='structure label'),
                                                   columns = []      )

    def _initialise_residue_table(self):
        """Initialise the tables.residues object with residue labels"""
        residue_labels = [make_label(c) for c in conformers_via_residue_groups(self.structures.hierarchies[0])]
        # Multiple observations of the same residue
        self.tables.residue_observations = pandas.Panel( data       = None,
                                                         items      = pandas.Index(data=residue_labels, name=['chain','residue','altloc'], tupleize_cols=True),
                                                         major_axis = pandas.Index(data=self.structures.labels, name='structure label'),
                                                         minor_axis = pandas.Index(data=[], name='variables')     )
        # Condensed statistics for each residue
        self.tables.residue_statistics = pandas.DataFrame( data = None,
                                                           index = pandas.Index(data=residue_labels, name=['chain','residue','altloc'], tupleize_cols=True),
                                                           columns = []      )

    ######################################################################################

    def load_all(self):
        self.load_all_structure_data()
        self.load_all_residue_data()
    def load_all_structure_data(self):
        self.extract_structure_info()
        self.calculate_structure_mean_b_factors()
    def load_all_residue_data(self):
        self.extract_residue_info()
        self.calculate_residue_mean_b_factors()
        self.calculate_residue_mean_normalised_b_factors()
        self.calculate_residue_phi_psi_angles()
        self.reduce_residue_bfactor_observations()
        self.reduce_residue_phi_psi_observations()

    def write_all(self, out_dir):
        self.write_tables(out_dir=out_dir)
        self.write_residue_statistics_as_structures(out_dir=out_dir)
        self.write_normalised_structures(out_dir=out_dir)
        self.write_graphs(out_dir=out_dir)

    ######################################################################################

    def write_tables(self, out_dir):
        """Dump all data as csv files in out_dir"""
        if not os.path.exists(out_dir): os.mkdir(out_dir)
        print('------------------------------------>')
        # Write csv for all global structure variables
        filename = os.path.join(out_dir, 'all_structures.csv')
        print('Writing {}'.format(filename))
        self.tables.structures.to_csv(filename)
        # Write csv for all global residue variables
        filename = os.path.join(out_dir, 'all_residues.csv')
        print('Writing {}'.format(filename))
        self.tables.residue_statistics.to_csv(filename)
        # Write csv for each variable for residues
        for variable in self.tables.residue_observations.minor_axis:
            filename = os.path.join(out_dir, 'all_'+variable+'.csv')
            print('Writing {}'.format(filename))
            self.tables.residue_observations.xs(variable, axis=2).T.to_csv(filename)

    def write_residue_statistics_as_structures(self, out_dir):
        """Write out structures with different residue statistics as the b-factors"""
        if not os.path.exists(out_dir): os.mkdir(out_dir)
        print('------------------------------------>')
        blank_h = self.structures.hierarchies[0].deep_copy()
        blank_h.atoms().set_b(flex.double([0]*blank_h.atoms_size()))
        cache = blank_h.atom_selection_cache()
        for variable in self.tables.residue_statistics.columns:
            out_pdb = os.path.join(out_dir, 'variation-{}.pdb'.format(variable))
            print('Writing {} as b-factors: {}'.format(variable, out_pdb))
            # Create new copy of the structure
            var_h = blank_h.deep_copy()

            for res_lab in self.tables.residue_statistics.index:
                # Get the value of the statistic for the residue
                statistic = self.tables.residue_statistics.get_value(res_lab, variable)
                if numpy.isnan(statistic): continue
                # Get the residue atoms for the residue
                selection = 'chain "{0}" and resid "{1}" and (altloc " " or altloc "{2}")'.format(*res_lab)
                residue = var_h.select(cache.selection(selection))
                # Set the b-factors of the residue to the statistic - WILL OVERRIDE OTHER CONFORMERS
                residue.atoms().set_b(flex.double([statistic]*residue.atoms_size()))
                # Write the output structure
                var_h.write_pdb_file(out_pdb)

    def write_normalised_structures(self, out_dir):
        """Write out the structures with normalised B-factors"""
        if not os.path.exists(out_dir): os.mkdir(out_dir)
        print('------------------------------------>')
        for lab_h, pdb_h in zip(self.structures.labels, self.structures.hierarchies):
            out_pdb = os.path.join(out_dir, '{}-normalised-b.pdb'.format(lab_h))
            print('Writing Normalised B-Factors: {}'.format(lab_h))
            pdb_h_z = normalise_b_factors_to_z_scores(pdb_hierarchy=pdb_h, method='protein')
            pdb_h_z.write_pdb_file(file_name=out_pdb)

    def write_graphs(self, out_dir):
        """Dump histograms in out_dir"""
        if not os.path.exists(out_dir): os.mkdir(out_dir)
        print('------------------------------------>')
        # Iterate through variables and plot for all of the structures
        for variable in self.tables.structures.columns:
            # Extract observations for this variable
            col_data = self.tables.structures[variable]
            filename = os.path.join(out_dir,'structure-{}.png'.format(variable))
            if col_data.any() and (col_data.nunique()>1):
                print('Writing {}'.format(filename))
                simple_histogram( filename = filename,
                                  data  = col_data.dropna(),
                                  title = 'histogram for {!s}'.format(variable),
                                  x_lab = variable  )
        print('------------------------------------>')
        # Iterate through variables and plot for all of the residues
        for variable in self.tables.residue_statistics.columns:
            col_data = self.tables.residue_statistics[variable]
            filename = os.path.join(out_dir,'residue-{}.png'.format(variable))
            if col_data.any() and (col_data.nunique()>1):
                col_data = col_data.dropna()
                print('Writing {}'.format(filename))
                simple_bar( filename = filename,
                            y_vals   = col_data.values,
                            x_labels = col_data.index,
                            title = 'bar for {!s}'.format(variable),
                            x_lab = 'residue',
                            y_lab = variable  )
        print('------------------------------------>')
        # Plot residue variables against each other for all of the residues
        for var_1 in self.tables.residue_statistics.columns:
            col_data_1 = self.tables.residue_statistics[var_1]
            if not col_data_1.any(): continue
            for var_2 in self.tables.residue_statistics.columns:
                if var_1 == var_2: break
                col_data_2 = self.tables.residue_statistics[var_2]
                if not col_data_2.any(): continue
                filename = os.path.join(out_dir,'residue-{}-{}.png'.format(var_1, var_2))
                print('Writing {}'.format(filename))
                simple_scatter( filename = filename,
                                x_vals = col_data_1.values,
                                y_vals = col_data_2.values,
                                title = 'bar for {!s} against {!s}'.format(var_1, var_2),
                                x_lab = var_1,
                                y_lab = var_2  )
        print('------------------------------------>')
        # Iterate through variables for observations of the residues
        for variable in self.tables.residue_observations.minor_axis:
            # Extract observations for this variable
            var_data = self.tables.residue_observations.xs(variable, 2)
            # Iterate through residues and plot for variable values
            for res_lab in self.tables.residue_observations.items:
                filename = os.path.join(out_dir,'{}-{}.png'.format(variable, '_'.join(res_lab)))
                col_data = var_data[res_lab]
                if col_data.any() and (col_data.nunique()>1):
                    print('Writing {}'.format(filename))
                    simple_histogram( filename = filename,
                                      data  = col_data.dropna(),
                                      title = '{} histogram for {!s}'.format(variable, res_lab),
                                      x_lab = variable  )

    ######################################################################################

    def extract_structure_info(self):
        """Extract information from the input objects, if given"""
        if not self.structures.inputs: raise Exception('No input objects given')
        # ----------------------------------------------------->
        self.tables.structures['r-work'] = numpy.nan
        self.tables.structures['r-free'] = numpy.nan
        self.tables.structures['res-high'] = numpy.nan
        self.tables.structures['res-low']  = numpy.nan
        # ----------------------------------------------------->
        print('------------------------------------>')
        for lab_h, pdb_i in zip(self.structures.labels, self.structures.inputs):
            print('Extracting Global PDB Metrics: {}'.format(lab_h))
            # Extract global statistics
            r_info = pdb_i.get_r_rfree_sigma()
            # Populate tables
            self.tables.structures.set_value(lab_h, 'r-work', r_info.r_work)
            self.tables.structures.set_value(lab_h, 'r-free', r_info.r_free)
            self.tables.structures.set_value(lab_h, 'res-high', r_info.high)
            self.tables.structures.set_value(lab_h, 'res-low',  r_info.low)

    def calculate_structure_mean_b_factors(self):
        """Calculate global B-factors for a structure"""
        # ----------------------------------------------------->
        self.tables.structures['mean-b-all']       = numpy.nan
        self.tables.structures['mean-b-protein']   = numpy.nan
        self.tables.structures['mean-b-backbone']  = numpy.nan
        self.tables.structures['mean-b-sidechain'] = numpy.nan
        # ----------------------------------------------------->
        self.tables.structures['rmsd-b-all']       = numpy.nan
        self.tables.structures['rmsd-b-protein']   = numpy.nan
        self.tables.structures['rmsd-b-backbone']  = numpy.nan
        self.tables.structures['rmsd-b-sidechain'] = numpy.nan
        # ----------------------------------------------------->
        print('------------------------------------>')
        for lab_h, pdb_h in zip(self.structures.labels, self.structures.hierarchies):
            print('Calculating Global Mean B-Factors: {}'.format(lab_h))
            # Extract B-factor statistics
            b_stats = BfactorStatistics.from_pdb(pdb_hierarchy=pdb_h.deep_copy())
            # Populate tables
            self.tables.structures.set_value(lab_h, 'mean-b-all',       b_stats.all.mean        )
            self.tables.structures.set_value(lab_h, 'mean-b-protein',   b_stats.protein.mean    )
            self.tables.structures.set_value(lab_h, 'mean-b-backbone',  b_stats.backbone.mean   )
            self.tables.structures.set_value(lab_h, 'mean-b-sidechain', b_stats.sidechain.mean  )
            # ----------------------------------------------------->
            self.tables.structures.set_value(lab_h, 'rmsd-b-all',       b_stats.all.biased_standard_deviation        )
            self.tables.structures.set_value(lab_h, 'rmsd-b-protein',   b_stats.protein.biased_standard_deviation    )
            self.tables.structures.set_value(lab_h, 'rmsd-b-backbone',  b_stats.backbone.biased_standard_deviation   )
            self.tables.structures.set_value(lab_h, 'rmsd-b-sidechain', b_stats.sidechain.biased_standard_deviation  )

    ######################################################################################

    def extract_residue_info(self):
        """Extract information about each of the residues"""
        # ----------------------------------------------------->
        self.tables.residue_observations.loc[:,:,'num_conformers'] = numpy.nan
        self.tables.residue_observations.loc[:,:,'num_atoms']      = numpy.nan
        self.tables.residue_observations.loc[:,:,'occupancy'] = numpy.nan
        # ----------------------------------------------------->
        print('------------------------------------>')
        for lab_h, pdb_h in zip(self.structures.labels, self.structures.hierarchies):
            print('Extracting Residue Info: {}'.format(lab_h))
            for c in conformers_via_residue_groups(pdb_h):
                res_lab = make_label(c)
                # Extract counts for residue group
                self.tables.residue_observations.set_value(res_lab, lab_h, 'num_conformers', len(c.parent().conformers()))
                self.tables.residue_observations.set_value(res_lab, lab_h, 'num_atoms', c.atoms_size())
                # Extract occupancy
                p_occ = [o for o in c.atoms().extract_occ() if o<1.0]
                if c.altloc and (len(set(p_occ)) == 1):
                    self.tables.residue_observations.set_value(res_lab, lab_h, 'occupancy', p_occ[0])

    def calculate_residue_mean_b_factors(self):
        """Extract Mean-B values in each of the structures"""
        # ----------------------------------------------------->
        self.tables.residue_observations.loc[:,:,'mean-b-all']       = numpy.nan
        self.tables.residue_observations.loc[:,:,'mean-b-backbone']  = numpy.nan
        self.tables.residue_observations.loc[:,:,'mean-b-sidechain'] = numpy.nan
        # ----------------------------------------------------->
        print('------------------------------------>')
        for lab_h, pdb_h in zip(self.structures.labels, self.structures.hierarchies):
            print('Calculating Local Mean B-Factors: {}'.format(lab_h))
            cache = pdb_h.atom_selection_cache()
            # Non-Hydrogens
            for c in conformers_via_residue_groups(s_select.non_h(hierarchy=pdb_h, cache=cache)):
                res_lab = make_label(c)
                res_mean_b = flex.mean_weighted(c.atoms().extract_b(), c.atoms().extract_occ())
                self.tables.residue_observations.set_value(res_lab, lab_h, 'mean-b-all', res_mean_b)
            # Backbone Atoms
            for c in conformers_via_residue_groups(s_select.backbone(hierarchy=pdb_h, cache=cache)):
                res_lab = make_label(c)
                res_mean_b = flex.mean_weighted(c.atoms().extract_b(), c.atoms().extract_occ())
                self.tables.residue_observations.set_value(res_lab, lab_h, 'mean-b-backbone', res_mean_b)
            # Sidechain Atoms
            for c in conformers_via_residue_groups(s_select.sidechains(hierarchy=pdb_h, cache=cache)):
                res_lab = make_label(c)
                res_mean_b = flex.mean_weighted(c.atoms().extract_b(), c.atoms().extract_occ())
                self.tables.residue_observations.set_value(res_lab, lab_h, 'mean-b-sidechain', res_mean_b)

    def calculate_residue_mean_normalised_b_factors(self):
        """Extract Mean-B values in each of the structures"""
        # ----------------------------------------------------->
        self.tables.residue_observations.loc[:,:,'mean-bz-all']       = numpy.nan
        self.tables.residue_observations.loc[:,:,'mean-bz-backbone']  = numpy.nan
        self.tables.residue_observations.loc[:,:,'mean-bz-sidechain'] = numpy.nan
        # ----------------------------------------------------->
        print('------------------------------------>')
        for lab_h, pdb_h in zip(self.structures.labels, self.structures.hierarchies):
            print('Calculating Local Normalised Mean B-Factors: {}'.format(lab_h))
            # Normalise the b-factors of the structure
            pdb_h_z = normalise_b_factors_to_z_scores(pdb_hierarchy=pdb_h, method='protein')
            cache = pdb_h_z.atom_selection_cache()
            # Non-Hydrogens
            for c in conformers_via_residue_groups(s_select.non_h(hierarchy=pdb_h_z, cache=cache)):
                res_lab = make_label(c)
                res_mean_b = flex.mean_weighted(c.atoms().extract_b(), c.atoms().extract_occ())
                self.tables.residue_observations.set_value(res_lab, lab_h, 'mean-bz-all', res_mean_b)
            # Backbone Atoms
            for c in conformers_via_residue_groups(s_select.backbone(hierarchy=pdb_h_z, cache=cache)):
                res_lab = make_label(c)
                res_mean_b = flex.mean_weighted(c.atoms().extract_b(), c.atoms().extract_occ())
                self.tables.residue_observations.set_value(res_lab, lab_h, 'mean-bz-backbone', res_mean_b)
            # Sidechain Atoms
            for c in conformers_via_residue_groups(s_select.sidechains(hierarchy=pdb_h_z, cache=cache)):
                res_lab = make_label(c)
                res_mean_b = flex.mean_weighted(c.atoms().extract_b(), c.atoms().extract_occ())
                self.tables.residue_observations.set_value(res_lab, lab_h, 'mean-bz-sidechain', res_mean_b)

    def calculate_normalised_normalised_residue_b_factors(self):
        """Calculate twice-normalised residue b-factors"""

        pass

    def reduce_residue_bfactor_observations(self):
        """Calculate statistics of the bfactor variation"""
        # ----------------------------------------------------->
        self.tables.residue_statistics['mean-mean-bz-all'] = numpy.nan
        self.tables.residue_statistics['rmsd-mean-bz-all'] = numpy.nan
        # ----------------------------------------------------->
        self.tables.residue_statistics['mean-mean-bz-backbone'] = numpy.nan
        self.tables.residue_statistics['rmsd-mean-bz-backbone'] = numpy.nan
        # ----------------------------------------------------->
        self.tables.residue_statistics['mean-mean-bz-sidechain'] = numpy.nan
        self.tables.residue_statistics['rmsd-mean-bz-sidechain'] = numpy.nan
        # ----------------------------------------------------->
        print('------------------------------------>')
        for res_lab in self.tables.residue_observations.items:
            print('Calculating Mean and RMS B-factors: {}'.format(res_lab))
            # Extract the residue data
            res_data = self.tables.residue_observations[res_lab]
            res_mean_bz_all  = res_data['mean-bz-all']
            res_mean_bz_back = res_data['mean-bz-backbone']
            res_mean_bz_side = res_data['mean-bz-sidechain']
            # All atoms
            if res_mean_bz_all.any():
                self.tables.residue_statistics.set_value(res_lab, 'mean-mean-bz-all', numpy.mean(res_mean_bz_all))
                self.tables.residue_statistics.set_value(res_lab, 'rmsd-mean-bz-all', numpy.std(res_mean_bz_all))
            # Backbone atoms
            if res_mean_bz_back.any():
                self.tables.residue_statistics.set_value(res_lab, 'mean-mean-bz-backbone', numpy.mean(res_mean_bz_back))
                self.tables.residue_statistics.set_value(res_lab, 'rmsd-mean-bz-backbone', numpy.std(res_mean_bz_back))
            # Sidechain atoms
            if res_mean_bz_side.any():
                self.tables.residue_statistics.set_value(res_lab, 'mean-mean-bz-sidechain', numpy.mean(res_mean_bz_side))
                self.tables.residue_statistics.set_value(res_lab, 'rmsd-mean-bz-sidechain', numpy.std(res_mean_bz_side))

    def calculate_residue_phi_psi_angles(self):
        """Extract phi-psi angles for each of the structures"""
        if not self.structures.hierarchies: raise Exception('No input objects given')
        # ----------------------------------------------------->
        self.tables.residue_observations.loc[:,:,'phi'] = numpy.nan
        self.tables.residue_observations.loc[:,:,'psi'] = numpy.nan
        # ----------------------------------------------------->
        print('------------------------------------>')
        for lab_h, pdb_h in zip(self.structures.labels, self.structures.hierarchies):
            print('Extracting Phi-Psi Angles: {}'.format(lab_h))
            # Note these are residue objects so model->chain->conformer->residue
            for res, (phi, psi) in get_all_phi_psi_for_hierarchy(pdb_h=pdb_h):
                # Create residue label - need to account for fact that all single confs are duplicated
                if res.is_pure_main_conf: res_lab = (res.parent().parent().id, res.resid(), '')
                else:                     res_lab = (res.parent().parent().id, res.resid(), res.parent().altloc)
                # Check to see if already done so only populate once
                if numpy.isnan(self.tables.residue_observations.get_value(res_lab, lab_h, 'phi')):
                    self.tables.residue_observations.set_value(res_lab, lab_h, 'phi', phi)
                # Populate values
                if numpy.isnan(self.tables.residue_observations.get_value(res_lab, lab_h, 'psi')):
                    self.tables.residue_observations.set_value(res_lab, lab_h, 'psi', psi)

    def reduce_residue_phi_psi_observations(self):
        """Calculate statistics of the phi-psi variation"""
        # ----------------------------------------------------->
        self.tables.residue_statistics['mean-phi'] = numpy.nan
        self.tables.residue_statistics['mean-psi'] = numpy.nan
        self.tables.residue_statistics['rmsd-phi'] = numpy.nan
        self.tables.residue_statistics['rmsd-psi'] = numpy.nan
        self.tables.residue_statistics['rms-rmsd-phi-rmsd-psi'] = numpy.nan
        self.tables.residue_statistics['max-rmsd-phi-rmsd-psi'] = numpy.nan
        # ----------------------------------------------------->
        print('------------------------------------>')
        for res_lab in self.tables.residue_observations.items:
            # Extract the residue data
            res_data = self.tables.residue_observations[res_lab]
            res_phi  = res_data['phi']
            res_psi  = res_data['psi']
            if (not res_phi.any()) and (not res_phi.any()): continue
            print('Calculating RMS Phi-Psi Angles: {}'.format(res_lab))
            # Analyse phi angles
            phi_m, phi_s = mean_std_angles(angles=res_phi)
            self.tables.residue_statistics.set_value(res_lab, 'mean-phi', phi_m)
            self.tables.residue_statistics.set_value(res_lab, 'rmsd-phi', phi_s)
            # Analyse psi angles
            psi_m, psi_s = mean_std_angles(angles=res_psi)
            self.tables.residue_statistics.set_value(res_lab, 'mean-psi', psi_m)
            self.tables.residue_statistics.set_value(res_lab, 'rmsd-psi', psi_s)
            # Combine into one value
            self.tables.residue_statistics.set_value(res_lab, 'rms-rmsd-phi-rmsd-psi', numpy.sqrt(psi_s**2 + phi_s**2))
            self.tables.residue_statistics.set_value(res_lab, 'max-rmsd-phi-rmsd-psi', numpy.max([psi_s, phi_s]))

