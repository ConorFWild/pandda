import os, copy

from collections import OrderedDict

import iotbx.pdb, iotbx.mtz, iotbx.ccp4_map

import scitbx.matrix

from bamboo.common import Meta, Info
from bamboo.common.file import FileManager
from bamboo.common.path import easy_directory

from giant.io.pdb import strip_pdb_to_input
from giant.xray.data import extract_structure_factors
from giant.xray.crystal import CrystalSummary
from giant.xray.symmetry import get_crystal_contact_operators, apply_symmetry_operators, combine_hierarchies
from giant.structure.align import align_structures_rigid, align_structures_flexible


class MultiCrystalDatasetPlain:
    def __init__(self, datasets=None, partitions=None):

        # Instantiate datasets
        self.datasets = datasets

        # Generate partitions on datasets
        if partitions is None:
            self.partitions = OrderedDict()
        else:
            self.partitions = partitions


        # Get max resolution
        self.max_res = max([d.data.summary.high_res
                            for dtag, d
                            in self.datasets.items()])

    def new_from_datasets(self, datasets):
        # Make a clone
        clone = MultiCrystalDataset(datasets=datasets,
                                    partitions=self.partitions.copy())

        # Remove datasets from partitions that aren't there anymore
        valid_dtags = [dtag for dtag in datasets]
        for partition in clone.partitions:
            clone.partitions[partition] = [dtag for dtag in clone.partitions[partition] if dtag in valid_dtags]
        return clone

    def set_partition(self, name, dtags):
        self.partitions[name] = dtags

    def get_partition(self, name):
        return self.partitions[name]

    def partition_datasets(self, partition):
        return {dtag: d
                for dtag, d
                in self.datasets.items()
                if dtag in self.partitions[partition]}




class MultiCrystalDataset:

    def __init__(self, datasets=None, dataloader=None, sample_loader=None, partitions=None):

        # self.trace = OrderedDict()

        # Instantiate datasets
        self.dataloader = dataloader
        if datasets is None:
            assert dataloader is not None
            self.datasets = self.dataloader()
            # self.trace[self.dataloader.name] = self.dataloader.log()
            # print("######### Dataloader summary #########")
            # for block, sub_block in self.dataloader.log().items():
            #     print("######### {} ##########".format(block))
            #     for sub_block, log_str in sub_block.items():
            #         print("# # {} # #".format(sub_block))
            #         print(log_str)
        else:
            self.datasets = datasets

        # Generate partitions on datasets
        if partitions is None:
            self.partitions = OrderedDict()
        else:
            self.partitions = partitions

        # Define sample laoder
        self.sample_loader = sample_loader

        # Set aside variable to store samples
        self.samples = {}

        # Get max resolution
        self.max_res = max([d.data.summary.high_res
                            for dtag, d
                            in self.datasets.items()])

    def new_from_datasets(self, datasets):
        # Make a clone
        clone = MultiCrystalDataset(datasets=datasets,
                                    sample_loader=self.sample_loader,
                                    partitions=self.partitions.copy())

        # Remove datasets from partitions that aren't there anymore
        valid_dtags = [dtag for dtag in datasets]
        for partition in clone.partitions:
            clone.partitions[partition] = [dtag for dtag in clone.partitions[partition] if dtag in valid_dtags]
        return clone

    def load_samples(self, res):
        self.samples = self.sample_loader(res, self.datasets)

    def load_samples_from_partition(self, res, partition):
        samples = self.sample_loader(res, self.partition_datasets(partition))
        self.samples.update(samples)

    def dump_samples(self):
        for dtag, sample in self.samples.items():
            del self.samples[dtag]
        del self.samples

    def set_partition(self, name, dtags):
        self.partitions[name] = dtags

    def get_partition(self, name):
        return self.partitions[name]

    def partition_datasets(self, partition):
        return {dtag: d
                for dtag, d
                in self.datasets.items()
                if dtag in self.partitions[partition]}

    def partition_samples(self, partition):
        return {dtag: samples
                for dtag, samples
                in self.samples.items()
                if dtag in self.partitions[partition]}


class _DatasetObj(object):


    def label(self, num=-1, tag=None):
        self.num = num
        if tag: self.tag = str(tag)
        return self


class AtomicModel(_DatasetObj):


    def __init__(self, input, hierarchy):
        self.input = input
        self.hierarchy = hierarchy
        self.filename = None
        self.alignment = None

    @classmethod
    def from_file(cls, filename):
#        ih = iotbx.pdb.hierarchy.input(filename)
        ih = strip_pdb_to_input(filename, remove_ter=True)
        c = cls(input=ih.input, hierarchy=ih.hierarchy)
        c.filename = filename
        return c

    @classmethod
    def from_other(cls, other):
        return cls(input=other.input, hierarchy=other.hierarchy)

    def align_to(self, other_hierarchy, method='local', **kwargs):
        """Align this model to another"""
        assert not self.alignment, 'Already aligned!'
        assert isinstance(other_hierarchy, iotbx.pdb.hierarchy.root)
        assert method in ['global','local'], 'alignment method not supported'
        if method == 'global':
            self.alignment = align_structures_rigid(mov_hierarchy=self.hierarchy, ref_hierarchy=other_hierarchy, **kwargs)
        else:
            self.alignment = align_structures_flexible(mov_hierarchy=self.hierarchy, ref_hierarchy=other_hierarchy, **kwargs)
        return self.alignment


class CrystallographicModel(AtomicModel):


    def __init__(self, input, hierarchy):
        super(CrystallographicModel, self).__init__(input=input, hierarchy=hierarchy)
        self.crystal_symmetry = self.input.crystal_symmetry()
        self.unit_cell = self.crystal_symmetry.unit_cell()
        self.space_group = self.crystal_symmetry.space_group()

        self._crystal_contacts_operators = None

    def crystal_contact_operators(self, distance_cutoff=10):
        """Return the symmetry operations to obtain the crystallographic copies within buffer distance of the atomic model"""
        return get_crystal_contact_operators(
                            hierarchy=self.hierarchy,
                            crystal_symmetry=self.crystal_symmetry,
                            distance_cutoff=distance_cutoff)

    def crystal_contacts(self, distance_cutoff=10, combine_copies=False):
        """Return the crystallographic symmetry copies within buffer distance of the atomic model"""
        ops = self.crystal_contact_operators(distance_cutoff=distance_cutoff)
        sym_hierarchies, chain_mappings = apply_symmetry_operators(
                                                    hierarchy=self.hierarchy,
                                                    crystal_symmetry=self.crystal_symmetry,
                                                    sym_ops_mat=ops)

        if combine_copies: sym_hierarchies = combine_hierarchies(sym_hierarchies)
        return sym_hierarchies


class ExperimentalData(_DatasetObj):


    @classmethod
    def from_file(cls, filename):
        if filename.endswith('.mtz'):
            return XrayData.from_file(filename=filename)


class XrayData(ExperimentalData):


    def __init__(self, mtz_object):
        self._mtz_object = mtz_object
        self.filename = None
        self.summary = CrystalSummary.from_mtz(mtz_object=mtz_object)
        self.miller_arrays = {}
        self.fft_maps = {}

    def _remove_mtz_objects(self):
        self._mtz_object = None
        self.summary._mtz_object = None

    @classmethod
    def from_file(cls, filename):
        assert filename.endswith('.mtz'), 'Given filename is not an mtz file'
        c = cls(mtz_object=iotbx.mtz.object(filename))
        c.filename = filename
        # Have the filename so no need to hold the object in memory -- allows pickling of object
        c._remove_mtz_objects()
        return c

    def mtz_object(self):
        if self._mtz_object:
            return self._mtz_object
        elif self.filename:
            return iotbx.mtz.object(self.filename)
        else:
            raise Exception('No filename to load data from')

    def get_structure_factors(self, columns):
        """Extract a let of structure factors from the mtz_object"""
        assert columns.count(',') == 1
        return extract_structure_factors(self.mtz_object(), ampl_label=columns.split(',')[0], phas_label=columns.split(',')[1])


class ModelAndData(_DatasetObj):


    def __init__(self, model, data):
        """Convenience object to hold model-data pairs"""
        self.model = model
        self.data  = data
        self.num = self.tag = None
        self.file_manager = None
        self.meta = Meta()
        self.parent = None
        self.child = None
        self.children = []

    @classmethod
    def from_file(cls, model_filename=None, data_filename=None):
        model = data = None
        if model_filename:
            assert os.path.exists(model_filename), 'Model file does not exist!'
            model = CrystallographicModel.from_file(filename=model_filename)
        if data_filename:
            assert os.path.exists(data_filename),  'Experimental Data file does not exist!'
            data = ExperimentalData.from_file(filename=data_filename)
        return cls(model=model, data=data)

    def initialise_output_directory(self, dir):
        """Initialise a dataset output directory"""
        # Create a file and directory organiser
        self.file_manager = FileManager(rootdir=easy_directory(dir))

    def get_pickle_copy(self):
        """Get copy of self that can be pickled - some cctbx objects cannot be pickled..."""

        return self


class PanddaDataset(ModelAndData):


    def __init__(self, model, data):
        """Subclass of ModelAndData used for PanDDA Analysis"""

        super(PanddaDataset, self).__init__(model=model, data=data)

        self.child = None
        self.events = []


class PanddaReferenceDataset(ModelAndData):


    _origin_shift = None

    def __init__(self, model, data):

        super(PanddaReferenceDataset, self).__init__(model=model, data=data)

        self.set_origin_shift((0.,0.,0.))
        self.child = None

    def set_origin_shift(self, shift):
        """Creates an alignment corresponding to an origin shift"""

        self._origin_shift = shift
        r = scitbx.matrix.rec([1,0,0,0,1,0,0,0,1], (3,3))
        t = scitbx.matrix.rec(shift, (3,1))
        rt = scitbx.matrix.rt((r,t))
        self.model.alignment = GlobalAlignment(alignment_mx=rt, alignment_sites=calphas(self.model.hierarchy).atoms().extract_xyz(),  id='ref')
        return self.model.alignment

    def origin_shift(self):
        return self._origin_shift

    def ref2grid(self, *args, **kwargs):
        return self.model.alignment.nat2ref(*args, **kwargs)
    def grid2ref(self, *args, **kwargs):
        return self.model.alignment.ref2nat(*args, **kwargs)

    def copy(self):
        return copy.deepcopy(self)
