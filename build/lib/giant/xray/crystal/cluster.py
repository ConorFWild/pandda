
import itertools
import numpy, pandas

import iotbx.mtz
import scipy.cluster.hierarchy

from bamboo.stats.cluster import generate_group_idxs

from giant.xray.crystal import CrystalSummary
from giant.xray.unit_cell import UnitCellVariation, lcv_from_unit_cells, pairwise_lcv, pairwise_ecv

class CrystalGroup(object):
    def __init__(self, crystals):
        """Class to hold multiple related crystals"""
        self.crystals = crystals
        self.uc_stats = UnitCellVariation([c.unit_cell for c in self.crystals])
        self.space_groups = list(set([c.space_group.info().symbol_and_number() for c in self.crystals]))
    def size(self):
        return len(self.crystals)
    def dendrogram(self, fname, method='lcv', **kwargs):
        """Print a dendrogram of the variation in the unit cells"""
        assert method in ['ecv','lcv'], 'method not recognised'
        if   method == 'lcv': link_func = pairwise_lcv
        elif method == 'ecv': link_func = pairwise_ecv
        return unit_cell_dendrogram(fname      = fname,
                                    unit_cells = [c.unit_cell for c in self.crystals],
                                    link_func  = link_func,
                                    labels     = [c.id for c in self.crystals],
                                    **kwargs)
    ######################################################
    @classmethod
    def by_space_group(cls, crystals):
        """Group crystals by space group and return CrystalGroup for each group"""
        sg_func = lambda c: c.space_group.info().symbol_and_number()
        return [cls(list(x[1])) for x in itertools.groupby(sorted(crystals, key=sg_func), key=sg_func)]
    @classmethod
    def by_unit_cell(cls, crystals, method='lcv', cutoff=0.5):
        """Cluster crystals by unit cell and return CrystalGroup for each cluster"""
        if len(crystals) == 1: return [cls(crystals)]
        assert method in ['lcv'], 'method not recognised'
#        # Method 1
#        if   method == 'lcv': link_func = lambda a,b: lcv_from_unit_cells(a.unit_cell, b.unit_cell)
#        hierarchy = libtbx.cluster.HierarchicalClustering(crystals, link_func)
#        clusters = hierarchy.getlevel(cutoff)
#        return [cls(c) for c in clusters]
        # Method 2
        if   method == 'lcv': link_func = pairwise_lcv
        dist_mat = link_func(unit_cells=[c.unit_cell for c in crystals])
        link_mat = scipy.cluster.hierarchy.linkage(dist_mat, method='single', metric='euclidean')
        clusters = scipy.cluster.hierarchy.fcluster(link_mat, t=cutoff, criterion='distance')
        return [cls([crystals[idx] for idx in g]) for i_g,g in generate_group_idxs(clusters)]

def dendrogram(fname, link_mat, labels=None, ylab=None, xlab=None, ylim=None, annotate_y_min=0.25, num_nodes=20):
    from matplotlib import pyplot
    fig = pyplot.figure()
    ax1 = pyplot.subplot(1,1,1)
    dend = scipy.cluster.hierarchy.dendrogram(link_mat, p=num_nodes, truncate_mode='lastp', labels=labels)
    # Change labels if requested
#    if labels: ax1.set_xticklabels([labels[i] for i in dend['leaves']])
    if xlab:   ax1.set_xlabel(xlab)
    if ylab:   ax1.set_ylabel(ylab)
    if ylim:   ax1.set_ylim(ylim)
    # Make sure the labels are rotated
    xlocs, xlabels = pyplot.xticks()
    pyplot.setp(xlabels, rotation=90)
    for i, d in zip(dend['icoord'], dend['dcoord']):
        x = 0.5 * sum(i[1:3])
        y = d[1]
        if y < annotate_y_min: continue
        pyplot.plot(x, y, 'ro')
        pyplot.annotate("%.3g" % y, (x, y), xytext=(0, -8), textcoords='offset points', va='top', ha='center')
    pyplot.tight_layout()
    fig.savefig(fname)
    return fig

def unit_cell_dendrogram(fname, unit_cells, link_func, labels=None, **kwargs):
    """Calculate a dendrogram for a set of mtz files/object, clustering by unit_cell parameters"""
    dist_mat = link_func(unit_cells=unit_cells)
    link_mat = scipy.cluster.hierarchy.linkage(dist_mat, method='single', metric='euclidean')
    dendrogram(fname=fname, link_mat=link_mat, labels=labels, **kwargs)

def sort_pdbs_by_space_group(pdb_files=None, pdb_inputs=None, labels=None):
    """Group the pdbfiles by spacegroup"""
    assert [pdb_files, pdb_inputs].count(None) == 1, 'Provide pdb_files OR pdb_inputs'
    if pdb_files: pdb_inputs = [None]*len(pdb_files)
    else:         pdb_files  = [None]*len(pdb_inputs)
    # Check or create labels
    if labels: assert len(labels) == len(pdb_files), 'labels must be the same length as pdb_files/pdb_objects'
    else:      labels = range(len(pdb_files))
    # Create summary objects
    pdb_summaries = [CrystalSummary.from_pdb(pdb_file=p_f, pdb_input=p_o, id=lab) for p_f,p_o,lab in zip(pdb_files, pdb_inputs, labels)]
    # Group by SpaceGroup
    return CrystalGroup.by_space_group(pdb_summaries)

def sort_mtzs_by_spacegroup(mtz_files=None, mtz_objects=None, labels=None):
    """Group the mtzfiles by spacegroup"""
    assert [mtz_files, mtz_objects].count(None) == 1, 'Provide mtz_files OR mtz_objects'
    if mtz_files: mtz_objects = [None]*len(mtz_files)
    else:         mtz_files   = [None]*len(mtz_objects)
    # Check or create labels
    if labels: assert len(labels) == len(mtz_files), 'labels must be the same length as mtz_files/mtz_objects'
    else:      labels = range(len(mtz_files))
    # Create summary objects
    mtz_summaries = [CrystalSummary.from_mtz(mtz_file=m_f, mtz_object=m_o, id=lab) for m_f,m_o,lab in zip(mtz_files, mtz_objects, labels)]
    # Group by SpaceGroup
    return CrystalGroup.by_space_group(mtz_summaries)

