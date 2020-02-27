from __future__ import print_function

from collections import OrderedDict

import traceback
import sys, copy, gc
import time
import logging

logger = logging.getLogger(__name__)

import numpy
from scipy import spatial

from joblib import Parallel, delayed

import iotbx.pdb, iotbx.mtz, iotbx.ccp4_map
import cctbx.maptbx, cctbx.uctbx

from libtbx import easy_mp
from libtbx.utils import Sorry, Failure
from libtbx.math_utils import ifloor, iceil

import cctbx
import cctbx_uctbx_ext
import scitbx.matrix
from scitbx.array_family import flex

from giant.xray.maps.scale import scale_map_to_reference
from multi_dataset_crystalography.grid import Grid, GridPartition
from multi_dataset_crystalography.grid.masks import AtomicMask, GridMask
from giant.structure.select import calphas, protein, sel_altloc, non_water

import joblib as jl

from bamboo.common import Meta, Info
from bamboo.common.holders import HolderList

from multi_dataset_crystalography.functions import wrapper_run

from multi_dataset_crystalography.dataset.reference import PanddaReferenceDataset

import dask
# from dask.distributed import worker_client


class MapLoaderDask(object):

    def __init__(self, verbose, resolution_factor, density_scaling):

        self.verbose = verbose
        self.resolution_factor = resolution_factor
        self.density_scaling = density_scaling

    def __call__(self, dataset, grid, reference_map, map_resolution):

        assert reference_map.is_sparse(), 'Reference map is not in sparse form'

        # ============================================================================>
        # Create map handler in the native coordinate frame
        # ============================================================================>
        # Extract the map data
        # TODO: make sure new version works
        # fft_map = dataset.data.fft_maps['truncated']
        dataset.data.fft_maps['truncated'] = dataset.data.miller_arrays['truncated'].fft_map(
            resolution_factor=float(self.resolution_factor),
            d_min=float(map_resolution),
            symmetry_flags=cctbx.maptbx.use_space_group_symmetry)
        fft_map = dataset.data.fft_maps['truncated']
        gc.collect()
        # Scale the map
        if self.density_scaling == 'none':
            pass
        elif self.density_scaling == 'sigma':
            fft_map.apply_sigma_scaling()
        elif self.density_scaling == 'volume':
            fft_map.apply_volume_scaling()
        # Create map object
        # native_map_true = ElectronDensityMap.from_fft_map(fft_map).as_map()
        native_map_true = ElectronDensityMap.from_fft_map(fft_map)

        # ============================================================================>
        # Morph the map to the reference frame
        # ============================================================================>
        # Extract the map sites from the grid partition
        point_mappings_grid = grid.partition.nn_groups[grid.global_mask().outer_mask_indices()]
        assert sum(point_mappings_grid == -1) == 0
        sites_cart_map = grid.grid2cart(grid.global_mask().outer_mask(),
                                        origin_shift=True,
                                        )
        # Translate the grid partition mappings to the dataset alignment mappings
        mappings_grid2dataset = get_interpolated_mapping_between_coordinates(query_list=grid.partition.sites_cart,
                                                                             ref_list=dataset.model.alignment.reference_sites,
                                                                             tol=0.01,
                                                                             )
        point_mappings_dataset = numpy.array([mappings_grid2dataset[i] for i in point_mappings_grid])
        assert sum(point_mappings_dataset == -1) == 0
        sites_cart_map_d = dataset.model.alignment.ref2nat(coordinates=sites_cart_map,
                                                           mappings=point_mappings_dataset,
                                                           )
        morphed_map_data = native_map_true.get_cart_values(sites_cart_map_d)

        # Scale map to reference
        scale_mask = grid.index_on_other(query=grid.global_mask().inner_mask_indices(),
                                         other=grid.global_mask().outer_mask_indices())
        scaled_map_data = scale_map_to_reference(ref_vals=reference_map.data,
                                                 vals=morphed_map_data,
                                                 mask_idxs=flex.size_t(scale_mask))
        # Create map holder
        morphed_map = reference_map.new_from_template(map_data=scaled_map_data, sparse=reference_map.is_sparse())
        morphed_map.meta.num = dataset.num
        morphed_map.meta.tag = dataset.tag
        morphed_map.meta.type = 'observed map'
        morphed_map.meta.resolution = reference_map.meta.resolution
        morphed_map.meta.map_uncertainty = None
        morphed_map.meta.obs_map_mean = morphed_map_data.min_max_mean().mean
        morphed_map.meta.obs_map_rms = morphed_map_data.standard_deviation_of_the_sample()
        morphed_map.meta.scl_map_mean = scaled_map_data.min_max_mean().mean
        morphed_map.meta.scl_map_rms = scaled_map_data.standard_deviation_of_the_sample()

        # Print a running row of dots
        print('>', end='');
        sys.stdout.flush()

        return morphed_map.make_sparse()

    def repr(self):
        repr = OrderedDict()
        return repr


class DefaultSampleLoader:

    def __init__(self, resolution_factor, density_scaling, cpus, verbose, grid_getter, reference=None,
                 multiprocessing="dask", grid=None, ref_map=None):
        self.resolution_factor = resolution_factor
        self.density_scaling = density_scaling

        self.cpus = cpus
        self.verbose = verbose
        self.grid_getter = grid_getter
        self.grid = grid
        self.reference = reference

        self.ref_map = ref_map
        self.truncated_reference = None
        self.truncated_datasets = None
        self.ref_map = ref_map
        self.multiprocessing = multiprocessing

    def __call__(self, cut_resolution, datasets):
        # ============================================================================>
        # Load maps for characterisation datasets
        # ============================================================================>
        gc.collect()
        pandda_load_and_morph_maps = PanddaLoadAndMorphMaps(self.resolution_factor, self.density_scaling,
                                                            self.cpus,
                                                            self.verbose,
                                                            multiprocessing=self.multiprocessing
                                                            )
        sample = pandda_load_and_morph_maps(datasets=self.truncated_datasets,
                                            ref_map=self.ref_map,
                                            map_resolution=cut_resolution,
                                            grid=self.grid)

        return sample

    def truncate_datasets(self, datasets):
        pandda_diffraction_data_truncater = PanddaDiffractionDataTruncater()
        truncated_reference, truncated_datasets = pandda_diffraction_data_truncater(datasets=datasets,
                                                                                    reference=self.reference)
        self.truncated_reference = truncated_reference
        self.truncated_datasets = truncated_datasets

    def get_reference(self, cut_resolution):
        # ============================================================================>
        # Load the reference map so that we can re-scale the individual maps to this
        # ============================================================================>
        pandda_reference_map_loader = PanddaReferenceMapLoader(self.resolution_factor, self.density_scaling)
        ref_map = pandda_reference_map_loader(self.reference,
                                              cut_resolution,
                                              self.grid)
        self.ref_map = ref_map

    def get_grid(self, reference):
        self.grid = self.grid_getter(reference)
        return self.grid

    def get_sample(self, cut_resolution, dataset):
        arg = MapLoader(dataset=dataset, grid=self.grid, reference_map=self.ref_map, verbose=self.verbose,
                        map_resolution=cut_resolution, resolution_factor=self.resolution_factor,
                        density_scaling=self.density_scaling)
        sample = wrapper_run(arg)

        return sample

    def instantiate(self, grid, ref_map):
        return DefaultSampleLoader(self.resolution_factor,
                                   self.density_scaling,
                                   self.cpus,
                                   self.verbose,
                                   self.grid_getter,
                                   reference=self.reference,
                                   multiprocessing="dask",
                                   grid=grid,
                                   ref_map=ref_map)

    def __repr__(self):
        repr = {"resolution_factor": self.resolution_factor,
                "density_scaling": self.density_scaling,
                }
        return repr


class PanddaDiffractionDataTruncater:

    def __init__(self):
        pass

    def __call__(self, datasets, reference):
        """Truncate data at the same indices across all the datasets"""

        # ==============================>
        # Find how many reflections are present in the reference dataset
        # ==============================>
        ref_cols = reference.meta.column_labels
        ref_size = reference.data.miller_arrays[ref_cols].set().size()

        # ==============================>
        # Truncate miller indices to the common set (not including the reference dataset)
        # ==============================>
        datasets = [d for dtag, d in datasets.items()]
        common_set = datasets[0].data.miller_arrays['scaled'].set()
        for dataset in datasets[1:]:
            common_set = common_set.common_set(dataset.data.miller_arrays['scaled'], assert_is_similar_symmetry=False)

        # ==============================>
        # Truncate diffraction data for all of the datasets (including the reference dataset)
        # ==============================>
        reference.data.miller_arrays['truncated'] = reference.data.miller_arrays[
            ref_cols].common_set(common_set, assert_is_similar_symmetry=False)
        for dataset in datasets:
            dataset.data.miller_arrays['truncated'] = dataset.data.miller_arrays['scaled'].common_set(common_set,
                                                                                                      assert_is_similar_symmetry=False)

        # TODO: Make sure that the reference and datasets modified in this whole class are copies rather than references
        truncated_reference = reference
        truncated_datasets = {d.tag: d for d in datasets}

        return truncated_reference, truncated_datasets

    def repr(self):
        repr = OrderedDict()
        return repr


class ElectronDensityMap(object):

    def __init__(self, map_data, unit_cell, map_indices=None, map_size=None, map_origin=(0.0, 0.0, 0.0), sparse=False,
                 meta=None, parent=None, children=None):

        assert isinstance(map_data, flex.double)
        assert isinstance(unit_cell, cctbx.uctbx.unit_cell) or isinstance(unit_cell, cctbx_uctbx_ext.unit_cell)

        if sparse:
            assert map_data.nd() == 1, 'Map data must be 1-dimensional when sparse=True'
            assert [map_indices, map_size].count(None) == 0, 'Must provide map_indices and map_size when sparse=True'
            assert len(map_data) == len(
                map_indices), 'map_data and map_indices must be the same length when sparse=True ({} != {})'.format(
                len(map_data), len(map_indices))
            assert max(map_indices) < numpy.prod(map_size), 'indices are not compatible with map_size ({} > {})'.format(
                max(map_indices), numpy.prod(map_size))
            if not isinstance(map_indices, flex.size_t):
                map_indices = flex.size_t(map_indices)
        else:
            if map_size is None:
                assert map_data.nd() == 3, 'map_data must be 3-dimension if map_size is not given'
                map_size = map_data.all()
            assert len(map_size) == 3, 'map_size must be three dimensional'
            assert map_indices is None, 'Do not provide map_indices for non-sparse matrices'
            assert numpy.prod(map_size) == map_data.size()
            # Reshape the map data if necessary
            if map_data.nd() == 1:
                map_data = map_data.deep_copy()
                map_data.reshape(flex.grid(map_size))
            assert map_data.all() == map_size, 'map_data is not the same shape as map_size ({} != {})'.format(
                map_data.all(), map_size)

        self.data = map_data
        self._map_size = map_size
        self._map_indices = map_indices
        self._map_origin = map_origin
        self.unit_cell = unit_cell
        self.meta = meta if meta else Meta()
        self.parent = parent
        self.children = children if children else []

        assert len(self._map_size) == 3, 'map_size must be tuple of length 3'
        assert sparse == self.is_sparse()

    @classmethod
    def from_fft_map(cls, fft_map):
        return cls(map_data=fft_map.real_map(), unit_cell=fft_map.unit_cell())

    def new_from_template(self, map_data, sparse=False, copy_meta=False, same_parent=False):
        """Create a new ElectronDensityMap using this map as a template"""

        # map_data is in sparse form
        if sparse:
            # Make sparse copy if template is not sparse
            if not self.is_sparse():
                self = self.copy().make_sparse()
            # Check input map data is compatible
            assert map_data.nd() == 1, 'map_data must 1-dimensional'
            assert map_data.size() == self._map_indices.size()
            # Extract parameters for sparseness
            map_size = self._map_size
            map_indices = self._map_indices
        else:
            assert map_data.size() == numpy.prod(self._map_size)
            map_size = self._map_size
            map_indices = None

        if copy_meta:
            meta = copy.deepcopy(self.meta)
        else:
            meta = None
        if same_parent:
            parent = self.parent
        else:
            parent = None

        return ElectronDensityMap(map_data=map_data,
                                  unit_cell=self.unit_cell,
                                  map_indices=map_indices,
                                  map_size=map_size,
                                  map_origin=self._map_origin,
                                  meta=meta,
                                  parent=parent,
                                  sparse=sparse)

    def copy(self):
        return self.new_from_template(map_data=self.data.deep_copy(),
                                      sparse=self.is_sparse(),
                                      copy_meta=True,
                                      same_parent=True)

    def normalised_copy(self):
        """Perform rms scaling on map data"""

        # Create output copy and make sparse (always calculate the mean and rms from the sparse values)
        result = self.copy().make_sparse()
        map_data = result.get_map_data(sparse=True)
        # Apply normalisation
        result.data = (result.data - numpy.mean(map_data)) * (1.0 / numpy.std(map_data))

        # Return the modified map
        if self.is_sparse():
            return result.make_sparse()
        else:
            return result.make_dense()

    def _check_compatibility(self, other):
        assert self.is_sparse() is other.is_sparse()

    def __add__(self, other):
        if isinstance(other, ElectronDensityMap):
            self._check_compatibility(other=other)
            return self.__add__(other.data)
        else:
            return self.new_from_template(map_data=self.data + other, sparse=self.is_sparse())

    def __sub__(self, other):
        if isinstance(other, ElectronDensityMap):
            self._check_compatibility(other=other)
            return self.__sub__(other.data)
        else:
            return self.new_from_template(map_data=self.data - other, sparse=self.is_sparse())

    def __mul__(self, other):
        if isinstance(other, ElectronDensityMap):
            self._check_compatibility(other=other)
            return self.__mul__(other.data)
        else:
            return self.new_from_template(map_data=self.data * other, sparse=self.is_sparse())

    def __div__(self, other):
        if isinstance(other, ElectronDensityMap):
            self._check_compatibility(other=other)
            return self.__div__(other.data)
        else:
            return self.new_from_template(map_data=self.data * (1.0 / other), sparse=self.is_sparse())

    def __rdiv__(self, other):
        return self.__div__(other)

    def is_sparse(self):
        return (self._map_indices is not None)

    def embed(self, map_data):
        """Embed map data relative to the real map origin, rather than (0,0,0)"""
        if self._map_origin == (0.0, 0.0, 0.0): return map_data
        return cctbx.maptbx.rotate_translate_map(unit_cell=self.unit_cell,
                                                 map_data=map_data,
                                                 rotation_matrix=scitbx.matrix.rec([1, 0, 0, 0, 1, 0, 0, 0, 1],
                                                                                   (3, 3)).elems,
                                                 translation_vector=(
                                                         -1.0 * scitbx.matrix.rec(self._map_origin, (3, 1))).elems)

    def as_map(self):
        map_data = self.get_map_data(sparse=False)
        return basic_map(
            cctbx.maptbx.basic_map_unit_cell_flag(),
            self.embed(map_data),
            map_data.focus(),
            self.unit_cell.orthogonalization_matrix(),
            cctbx.maptbx.out_of_bounds_clamp(0).as_handle(),
            self.unit_cell)

    def get_cart_values(self, cart_points):
        assert not self.is_sparse(), 'map must not be in sparse format for sampling'
        # Shift input points to the grid frame -- TODO implement in function so that rotations can be automated integrated
        cart_points = (cart_points - self._map_origin)
        frac_values = self.unit_cell.fractionalize(cart_points)
        # Get the map data with the correct origin
        map_data = self.get_map_data(sparse=False)
        map_vals = map(map_data.eight_point_interpolation, frac_values)
        return flex.double(map_vals)

    def to_file(self, filename, space_group):
        map_data = self.get_map_data(sparse=False)
        iotbx.ccp4_map.write_ccp4_map(
            file_name=filename,
            unit_cell=self.unit_cell,
            space_group=space_group,
            map_data=self.embed(map_data),
            labels=flex.std_string(['Output map from giant/pandda']))

    def get_map_data(self, sparse):
        """Get the map data as sparse/dense without altering state of master object"""
        if sparse is not self.is_sparse():
            result = self.copy()
            if sparse:
                result.make_sparse()
            else:
                result.make_dense()
        else:
            result = self

        return result.data

    def make_sparse(self):
        """Convert the map data into sparse form"""
        if self.is_sparse(): return self

        data_flat = self.data.as_1d()
        data_mask = (data_flat != 0.0)
        sparse_idxs = data_mask.iselection()
        sparse_data = data_flat.select(data_mask)

        self.data = sparse_data
        self._map_indices = sparse_idxs

        return self

    def make_dense(self):
        """Convert the map data into dense form"""
        if not self.is_sparse(): return self

        data_bulk = numpy.zeros(numpy.product(self._map_size))
        data_bulk.put(indices=self._map_indices, values=self.data)
        data_bulk = flex.double(data_bulk)
        data_bulk.reshape(flex.grid(self._map_size))

        self.data = data_bulk
        self._map_indices = None

        return self


class PanddaReferenceMapLoader:

    def __init__(self, resolution_factor, density_scaling):
        self.resolution_factor = resolution_factor
        self.density_scaling = density_scaling

    def __call__(self, reference, map_resolution, grid):
        """Load the reference map, and calculate some map statistics"""

        # ==============================>
        # Take the scaled diffraction data for the reference dataset and create fft
        # ==============================>
        # ref_dataset = self.datasets.reference()
        ref_dataset = reference
        fft_map = ref_dataset.data.miller_arrays['truncated'].fft_map(
            resolution_factor=float(self.resolution_factor),
            d_min=float(map_resolution),
            symmetry_flags=cctbx.maptbx.use_space_group_symmetry)
        # ==============================>
        # Scale the map
        # ==============================>
        if self.density_scaling == 'none':
            pass
        elif self.density_scaling == 'sigma':
            fft_map.apply_sigma_scaling()
        elif self.density_scaling == 'volume':
            fft_map.apply_volume_scaling()
        # ==============================>
        # Transform to the reference frame
        # ==============================>
        # Extract the points for the map (in the grid frame)
        masked_cart = grid.grid2cart(grid.global_mask().outer_mask(), origin_shift=True)
        # Create map handler in the native frame and extract the map values
        ref_map_true = ElectronDensityMap.from_fft_map(fft_map)
        masked_vals = ref_map_true.get_cart_values(masked_cart)
        # ==============================>
        # Create a new electron density map object for the "grid map"
        # ==============================>
        ref_map = ElectronDensityMap(map_data=masked_vals, unit_cell=grid.unit_cell(),
                                     map_indices=grid.global_mask().outer_mask_indices(),
                                     map_size=grid.grid_size(),
                                     map_origin=grid.cart_origin(),
                                     sparse=True)
        # Store the map as a child of the dataset
        ref_dataset.child = ref_map
        # ==============================>
        # Add some meta for debugging, etc
        # ==============================>
        ref_map.meta.type = 'reference-map'
        ref_map.meta.resolution = map_resolution
        ref_map.meta.map_mean = ref_map.get_map_data(sparse=True).min_max_mean().mean
        ref_map.meta.map_rms = ref_map.get_map_data(sparse=True).standard_deviation_of_the_sample()

        return ref_map

    def repr(self):
        repr = OrderedDict()
        repr["resolution_factor"] = self.resolution_factor
        repr["density_scaling"] = self.density_scaling
        return repr


class MapHolderList(HolderList):
    _holder_class = ElectronDensityMap

    def _get_num(self, item):
        return item.meta.num

    def _get_tag(self, item):
        return item.meta.tag


class PanddaLoadAndMorphMaps:

    def __init__(self, resolution_factor, density_scaling, cpus, verbose, multiprocessing):

        self.resolution_factor = resolution_factor
        self.density_scaling = density_scaling
        self.cpus = cpus
        self.verbose = verbose
        self.multiprocessing = multiprocessing

    def __call__(self, datasets, ref_map, map_resolution, grid):
        """Create map from miller arrays. Transform map into the reference frame by sampling at the given points."""

        assert ref_map.is_sparse(), 'Reference map is not in sparse form'

        # ==============================>
        # Create holder for the output map objects
        # ==============================>
        sample = {}
        # ==============================>
        # Return empty list if no datasets
        # ==============================>
        if not datasets: return sample

        # ==============================>
        # Load maps in parallel
        # ==============================>
        print('Loading maps (using {!s} cores)'.format(self.cpus))
        arg_list = [
            MapLoader(dataset=d, grid=grid, reference_map=ref_map, verbose=self.verbose,
                      map_resolution=map_resolution, resolution_factor=self.resolution_factor,
                      density_scaling=self.density_scaling)
            for dtag, d in datasets.items()]
        res = arg_list[0].run()

        # Print a sort of progress bar
        print('1' + ''.join(['{:<5}'.format(i) for i in range(0, len(arg_list) + 5, 5)])[2:])
        print(' ' * len(arg_list) + '|\r', end='')
        sys.stdout.flush()
        gc.collect()
        if self.multiprocessing == "dask":
            with worker_client(timeout=120, separate_thread=False) as client:
                dataset_maps_futures = client.map(wrapper_run, arg_list)
                dataset_maps = client.gather(dataset_maps_futures)

            # results = []
            # for arg in arg_list:
            #     y = dask.delayed(wrapper_run)(arg)
            #     results.append(y)
            # dataset_maps = dask.compute(results)
            # print(dask.distributed.get_worker())
            #
            # client = dask.distributed.get_client()
            # map_futures = client.map(wrapper_run, arg_list)
            # dask.distributed.secede()
            # dataset_maps = client.gather(map_futures)
            # dask.distributed.rejoin()

        else:
            dataset_maps = jl.Parallel(n_jobs=self.cpus,
                                       verbose=5)(jl.delayed(wrapper_run)(arg)
                                                  for arg
                                                  in arg_list)
        # ==============================>
        # Managed
        # ==============================>
        print('|')
        sample = {m.meta.tag: m
                  for m
                  in dataset_maps}
        # ==============================>
        # Clear fft map data to save memory
        # ==============================>
        # for dtag, m in sample.items():
        #     # TODO: is this the best way of handling this now?
        #     map_dataset = datasets[m.meta.tag]
        #     map_dataset.data.fft_maps['truncated'] = None

        return sample


class MapLoader(object):

    def __init__(self, dataset, grid, reference_map, verbose, map_resolution, resolution_factor, density_scaling):
        """
        The main object for loading the maps for PanDDA.
        Constructed so that the class can be initialised and then called within a multiprocessing function.

        Usage:
            new    = MapLoader(...)
            output = new.run()
        or:
            output = MapLoader.process(...)
        """

        self.data = (dataset, grid, reference_map, verbose, map_resolution, resolution_factor, density_scaling)

    # @classmethod
    # def process(cls, dataset, grid, reference_map, args, verbose):
    #     """Process the dataset immediately and return output"""
    #     return cls(dataset=dataset, grid=grid, reference_map=reference_map, args=args, verbose=verbose)

    def run(self):

        dataset, grid, reference_map, verbose, map_resolution, resolution_factor, density_scaling = self.data

        # log_file = dataset.file_manager.get_file('dataset_log')

        assert reference_map.is_sparse(), 'Reference map is not in sparse form'

        # ============================================================================>
        # Create map handler in the native coordinate frame
        # ============================================================================>
        # Extract the map data
        # TODO: make sure new version works
        # fft_map = dataset.data.fft_maps['truncated']
        dataset.data.fft_maps['truncated'] = dataset.data.miller_arrays['truncated'].fft_map(
            resolution_factor=float(resolution_factor),
            d_min=float(map_resolution),
            symmetry_flags=cctbx.maptbx.use_space_group_symmetry)
        fft_map = dataset.data.fft_maps['truncated']
        gc.collect()
        # Scale the map
        if density_scaling == 'none':
            pass
        elif density_scaling == 'sigma':
            fft_map.apply_sigma_scaling()
        elif density_scaling == 'volume':
            fft_map.apply_volume_scaling()
        # Create map object
        # native_map_true = ElectronDensityMap.from_fft_map(fft_map).as_map()
        native_map_true = ElectronDensityMap.from_fft_map(fft_map)

        # ============================================================================>
        # Morph the map to the reference frame
        # ============================================================================>
        # Extract the map sites from the grid partition
        point_mappings_grid = grid.partition.nn_groups[grid.global_mask().outer_mask_indices()]
        assert sum(point_mappings_grid == -1) == 0
        sites_cart_map = grid.grid2cart(grid.global_mask().outer_mask(), origin_shift=True)
        # Translate the grid partition mappings to the dataset alignment mappings
        mappings_grid2dataset = get_interpolated_mapping_between_coordinates(query_list=grid.partition.sites_cart,
                                                                             ref_list=dataset.model.alignment.reference_sites,
                                                                             tol=0.01)
        point_mappings_dataset = numpy.array([mappings_grid2dataset[i] for i in point_mappings_grid])
        assert sum(point_mappings_dataset == -1) == 0
        sites_cart_map_d = dataset.model.alignment.ref2nat(coordinates=sites_cart_map, mappings=point_mappings_dataset)
        morphed_map_data = native_map_true.get_cart_values(sites_cart_map_d)

        # Scale map to reference
        scale_mask = grid.index_on_other(query=grid.global_mask().inner_mask_indices(),
                                         other=grid.global_mask().outer_mask_indices())
        scaled_map_data = scale_map_to_reference(ref_vals=reference_map.data,
                                                 vals=morphed_map_data,
                                                 mask_idxs=flex.size_t(scale_mask))
        # Create map holder
        morphed_map = reference_map.new_from_template(map_data=scaled_map_data, sparse=reference_map.is_sparse())
        morphed_map.meta.num = dataset.num
        morphed_map.meta.tag = dataset.tag
        morphed_map.meta.type = 'observed map'
        morphed_map.meta.resolution = reference_map.meta.resolution
        morphed_map.meta.map_uncertainty = None
        morphed_map.meta.obs_map_mean = morphed_map_data.min_max_mean().mean
        morphed_map.meta.obs_map_rms = morphed_map_data.standard_deviation_of_the_sample()
        morphed_map.meta.scl_map_mean = scaled_map_data.min_max_mean().mean
        morphed_map.meta.scl_map_rms = scaled_map_data.standard_deviation_of_the_sample()

        # Print a running row of dots
        print('>', end='');
        sys.stdout.flush()

        return morphed_map.make_sparse()


class PanDDAGridSetup:

    def __init__(self, cpus, mask_pdb, align_mask_to_reference, alignment_method,
                 outer_mask, inner_mask, inner_mask_symmetry, grid_spacing, padding, verbose, mask_selection_string):

        self.cpus = cpus

        self.mask_pdb = mask_pdb
        self.align_mask_to_reference = align_mask_to_reference
        self.alignment_method = alignment_method

        self.outer_mask = outer_mask
        self.inner_mask = inner_mask
        self.inner_mask_symmetry = inner_mask_symmetry

        self.grid_spacing = grid_spacing
        self.padding = padding
        self.verbose = verbose

        self.mask_selection_string = mask_selection_string

        self.grid = None

    def __call__(self, reference=None):
        """Generate the grid objects for the analysis"""

        # ============================================================================>
        #####
        # Create Sampling Grid (for generated maps)
        #####
        # Create reference grid based on the reference structure
        # ============================================================================>
        if self.grid is None:
            # Which dataset to be used to mask the grid
            if bool(self.mask_pdb):
                mask_dataset = PanddaReferenceDataset.from_file(model_filename=self.mask_pdb).label(
                    tag='masking')
                if self.align_mask_to_reference:
                    try:
                        mask_dataset.model.alignment = None
                        mask_dataset.model.align_to(other_hierarchy=reference.model.hierarchy,
                                                    method=self.alignment_method,
                                                    require_hierarchies_identical=False)
                    except:
                        msg = traceback.format_exc()
                        msg += '\n------------------>>>'
                        msg += '\n\nFailed to align masking pdb ({}) to the reference structure.'.format(
                            self.mask_pdb)
                        msg += '\nIf the masking structure does not need alignment, rerun with params.masks.align_mask_to_reference=False'
                        raise Failure(msg)
                else:
                    mask_dataset.set_origin_shift((0.0, 0.0, 0.0))
            else:
                mask_dataset = reference.copy()

            # Create the grid using the masking dataset (for determining size and extent of grid)
            print("Creating referene grid")
            self.create_reference_grid(dataset=mask_dataset, grid_spacing=self.grid_spacing, reference=reference)
            print("Masking reference grid")
            self.mask_reference_grid(dataset=mask_dataset, selection=self.mask_selection_string)
            # Store the transformation to shift the reference dataset to the "grid frame", where the grid origin is (0,0,0)
            print("shifting reference origin")
            reference.set_origin_shift([-1.0 * a for a in self.grid.cart_origin()])
            # Partition the grid with the reference dataset (which grid points use which coordinate transformations)
            # print("Partitioning reference grid")
            self.partition_reference_grid(dataset=reference)

        return self.grid

    def create_reference_grid(self, dataset, grid_spacing, reference):
        """Create a grid over the given dataset"""

        # ============================================================================>
        # Extract sites and calculate grid extent
        # ============================================================================>
        sites_cart = dataset.model.alignment.nat2ref(dataset.model.hierarchy.atoms().extract_xyz())
        # Calculate the extent of the grid
        buffer = self.outer_mask + self.padding
        grid_min = flex.double([s - buffer for s in sites_cart.min()])
        grid_max = flex.double([s + buffer for s in sites_cart.max()])

        # ============================================================================>
        # Create main grid object
        # ============================================================================>
        self.grid = Grid(grid_spacing=grid_spacing,
                         origin=tuple(grid_min),
                         approx_max=tuple(grid_max),
                         verbose=self.verbose)
        # ==============================>
        # Calculate alignment between reference dataset and grid
        # ==============================>
        ref_dataset = reference
        ref_dataset.set_origin_shift(-1.0 * grid_min)
        # Write out masks if selected
        # if self.args.output.developer.write_grid_frame_masks:
        #     f_name = splice_ext(self.file_manager.get_file('reference_structure'), 'grid', position=-1)
        #     if not os.path.exists(f_name):
        #         tmp_h = ref_dataset.model.hierarchy.deep_copy()
        #         tmp_h.atoms().set_xyz(ref_dataset.ref2grid(tmp_h.atoms().extract_xyz()))
        #         tmp_h.write_pdb_file(f_name)
        #     f_name = splice_ext(self.file_manager.get_file('reference_symmetry'), 'grid', position=-1)
        #     if not os.path.exists(f_name):
        #         tmp_h = ref_dataset.model.crystal_contacts(distance_cutoff=self.args.params.masks.outer_mask+5, combine_copies=True)
        #         tmp_h.atoms().set_xyz(ref_dataset.ref2grid(tmp_h.atoms().extract_xyz()))
        #         tmp_h.write_pdb_file(f_name)

        return self.grid

    def mask_reference_grid(self, dataset, selection=None):
        """Create masks for the reference grid based on distances from atoms in the reference structure"""

        # ============================================================================>
        # Get main and neighbouring symmetry copies of the masking structure
        # ============================================================================>

        ref_h = dataset.model.hierarchy
        sym_h = dataset.model.crystal_contacts(distance_cutoff=self.outer_mask + 5.0, combine_copies=True)
        # ============================================================================>
        # Apply mask (protein=default if selection is not given)
        # ============================================================================>
        if selection:
            ref_h = ref_h.select(ref_h.atom_selection_cache().selection(selection), copy_atoms=True)
        else:
            ref_h = protein(ref_h)
        # ============================================================================>
        # Always generate symmetry mask using all non-water atoms - TODO also allow custom definitions? TODO
        # ============================================================================>
        sym_h = non_water(sym_h)
        # ============================================================================>
        # Check that these contain atoms
        # ============================================================================>
        if len(ref_h.atoms()) == 0: raise Sorry('Zero atoms have been selected to mask the grid')
        if len(sym_h.atoms()) == 0: raise Sorry('Zero atoms have been selected to mask the grid')
        # ============================================================================>
        # Extract coordinates
        # ============================================================================>
        ref_sites_cart = dataset.model.alignment.nat2ref(ref_h.atoms().extract_xyz())
        sym_sites_cart = dataset.model.alignment.nat2ref(sym_h.atoms().extract_xyz())
        # ============================================================================>
        # Global mask used for removing points in the bulk solvent regions
        # ============================================================================>
        if self.grid.global_mask() is None:
            global_mask = AtomicMask(parent=self.grid,
                                     sites_cart=ref_sites_cart,
                                     max_dist=self.outer_mask,
                                     min_dist=self.inner_mask)
            self.grid.set_global_mask(global_mask)
        # ============================================================================>
        # Global mask used for removing points close to symmetry copies of the protein
        # ============================================================================>
        if self.grid.symmetry_mask() is None:
            symmetry_mask = GridMask(parent=self.grid,
                                     sites_cart=sym_sites_cart,
                                     max_dist=self.outer_mask,
                                     min_dist=self.inner_mask_symmetry)
            self.grid.set_symmetry_mask(symmetry_mask)
        # ============================================================================>
        # Write masked maps
        # ============================================================================>
        # # Write protein masked map
        # indices = self.grid.global_mask().total_mask_indices()
        # f_name = self.file_manager.get_file('reference_dataset').replace('.mtz','.totalmask.ccp4')
        # if self.args.output.developer.write_grid_frame_masks:
        #     self.grid.write_indices_as_map(indices=indices, f_name=splice_ext(f_name, 'grid', position=-1), origin_shift=False)
        # if 1 or self.args.output.developer.write_reference_frame_common_masks_and_maps:
        #     self.grid.write_indices_as_map(indices=indices, f_name=splice_ext(f_name, 'ref', position=-1), origin_shift=True)
        #
        # # Write symmetry masked map
        # indices = self.grid.symmetry_mask().total_mask_indices()
        # f_name = self.file_manager.get_file('reference_dataset').replace('.mtz','.symmask.ccp4')
        # if self.args.output.developer.write_grid_frame_masks:
        #     self.grid.write_indices_as_map(indices=indices, f_name=splice_ext(f_name, 'grid', position=-1), origin_shift=False)
        # if 1 or self.args.output.developer.write_reference_frame_common_masks_and_maps:
        #     self.grid.write_indices_as_map(indices=indices, f_name=splice_ext(f_name, 'ref', position=-1), origin_shift=True)

        return self.grid

    def partition_reference_grid(self, dataset, altlocs=['', 'A']):

        # ============================================================================>
        # Select the sites for generating the voronoi alignments (calphas)
        # ============================================================================>
        partition_h = calphas(sel_altloc(dataset.model.hierarchy, altlocs=altlocs))
        site_cart_ca = partition_h.atoms().extract_xyz()
        # ============================================================================>
        # Create voronoi cells based on these atoms
        # ============================================================================>
        t1 = time.time()
        # self.grid.create_grid_partition(sites_cart=site_cart_ca)
        # self.grid.partition.partition(mask  = self.grid.global_mask(),
        #                               cpus  = self.cpus)
        self.grid.partition = partition_grid(self.grid,
                                             reference_dataset=dataset,
                                             executor=ExecutorJoblib(cpus=21),
                                             )
        t2 = time.time()
        # ============================================================================>
        # Print cell-by-cell summary or the partitioning
        # ============================================================================>
        # self.log.bar(True, False)
        # self.log('Partition Summary:')
        # self.log.bar(False, True)
        # voronoi_counts = dict(zip(*numpy.unique(self.grid.partition.nn_groups, return_counts=True)))
        # # Cell-by-Cell summary of the voronoi cells
        # self.log.bar()
        # self.log('CHN - RES -  RESID  - ATOM - ALT :    VORONOI VOLUME')
        # self.log.bar()
        # for i_atom, atom in enumerate(partition_h.atoms_with_labels()):
        #     self.log('{:<3} - {:<3} - {:<7} - {:<4} - {:<3} : {:>10} points'.format(atom.chain_id, atom.resname, atom.resid(), atom.name, atom.altloc, voronoi_counts.get(i_atom,0)))
        # self.log.bar()
        # self.log('Unpartitioned space: {} points'.format(voronoi_counts.get(-1,0)))
        # self.log.bar()
        # # Chain-by-chain summary of the voronoi cells
        # for c in partition_h.chains():
        #     self.log('Chain {:1} - {:5} regions - ({:5} residues)'.format(c.id, len(c.atoms()), len(c.residue_groups())))
        # self.log.bar()
        # self.log('Total: {} regions ({} chains, {} residues)'.format(len(partition_h.atoms()), len(list(partition_h.chains())), len(list(partition_h.residue_groups()))))

        return self.grid

    def __repr__(self):
        repr = {}
        return repr


def partition_grid(grid, reference_dataset, executor, mask=None, altlocs=['', 'A']):
    logger.info("Partitioning Calphas")
    partition_h = calphas(sel_altloc(reference_dataset.model.hierarchy, altlocs=altlocs))
    logger.info("Getting site carts")
    site_cart_ca = partition_h.atoms().extract_xyz()

    logger.info("making partition")
    grid_partition = GridPartition(grid,
                                   site_cart_ca,
                                   )

    # assert isinstance(cpus, int) and (cpus > 0)

    # Sites that we are partitioning
    logger.info("querying sites")
    if mask:
        query_sites = flex.vec3_double(mask.outer_mask())
    else:
        query_sites = flex.vec3_double(grid.grid_points())
    # Find the nearest grid_site for each query_site (returns index of the grid site)
    print("STARTING MULTIPROCESSING")
    if executor.cpus == 1:
        output = [find_sites((grid_partition.sites_grid, query_sites))]
    else:
        # Chunk the points into groups
        chunk_size = iceil(1.0 * len(query_sites) / executor.cpus)
        chunked_points = [query_sites[i:i + chunk_size]
                          for i
                          in range(0, len(query_sites), chunk_size)]
        assert sum(map(len, chunked_points)) == len(query_sites)
        assert len(chunked_points) == executor.cpus
        # Map to cpus
        # arg_list = [(self.sites_grid, chunk) for chunk in chunked_points]
        # output = easy_mp.pool_map(fixed_func=find_sites, args=arg_list, processes=cpus)
        funcs = []
        for i, chunk in enumerate(chunked_points):
            f = FindSites(sites_grid=grid_partition.sites_grid,
                          chunk=chunk,
                          )
            funcs.append(f
                         )

        output = executor(funcs)

    assert len(output) == executor.cpus, '{!s} != {!s}'.format(len(output), executor.cpus)

    logger.info("finding sites")
    # output = executor(find_sites,
    #                   [(grid_partition.sites_grid, [query_site])
    #                    for query_site
    #                    in query_sites
    #                    ],
    #                   )
    # funcs = [lambda: find_sites((grid_partition.sites_grid,
    #                              [query_site],
    #                              )
    #                             )
    #          for query_site
    #          in query_sites
    #          ]
    # output = executor(funcs)

    # assert len(output) == cpus, '{!s} != {!s}'.format(len(output), cpus)
    # Extract the indices of the mapped points
    logger.info("nn_groups")
    nn_groups = []
    [nn_groups.extend(o) for o in output]
    nn_groups = numpy.array(nn_groups)
    logger.info("assering")
    logger.info(output)

    assert len(query_sites) == len(nn_groups)
    logger.info("assertation passed")

    # Reformat into full grid size
    if mask:
        grid_partition.nn_groups = -1 * numpy.ones(grid.grid_size_1d(), dtype=int)
        grid_partition.nn_groups.put(mask.outer_mask_indices(), nn_groups)
    else:
        grid_partition.nn_groups = nn_groups

    logger.info("returning from nn groups")

    return grid_partition


def query_by_grid_indices(self, idxs):
    """Return the atom label for a grid site index"""
    assert self.nn_groups is not None, 'Grid not yet partitioned'
    return numpy.array([self.nn_groups[i] for i in idxs])


def query_by_grid_points(self, gps):
    """Return the atom label for a grid point"""
    assert self.nn_groups is not None, 'Grid not yet partitioned'
    indxr = self.grid.indexer()
    return numpy.array([self.nn_groups[indxr(g)] for g in gps])


def query_by_cart_points(self, sites_cart):
    """Dynamically calculate the nearest atom site to the input points"""
    tree = spatial.KDTree(data=self.sites_cart)
    nn_dists, nn_groups = tree.query(sites_cart)
    return numpy.array(nn_groups)


def find_sites(sites_tuple):
    ref_sites, query_sites = sites_tuple
    tree = spatial.KDTree(data=ref_sites)
    nn_dists, nn_groups = tree.query(query_sites)
    return nn_groups


def get_interpolated_mapping_between_coordinates(query_list, ref_list, tol=0.01):
    """
    Take each of query_list and find the sites in ref_list (within tolerance).
    Missing sites will be interpolated to the closest neighbouring site.
    Return list of indices mapping the site in one to the closest site in the other.
    """
    ref_list = flex.vec3_double(ref_list)
    tmp_idxs_q_to_r = [closest_point_within_tolerance(query=q, ref_list=ref_list, tol=tol) for q in query_list]
    assert tmp_idxs_q_to_r.count(-1) != len(tmp_idxs_q_to_r), 'no matching sites found between mappings'
    out_idxs_q_to_r = copy.copy(tmp_idxs_q_to_r)
    l = len(tmp_idxs_q_to_r)
    # Populate the missing values with the nearest value
    for i in range(l):
        d_i = 0
        while out_idxs_q_to_r[i] == -1:
            d_i += 1;
            p_i = i + d_i;
            n_i = i - d_i
            if (p_i < l) and (tmp_idxs_q_to_r[p_i] != -1):
                out_idxs_q_to_r[i] = out_idxs_q_to_r[p_i]
            elif (n_i >= 0) and (tmp_idxs_q_to_r[n_i] != -1):
                out_idxs_q_to_r[i] = out_idxs_q_to_r[n_i]

    return out_idxs_q_to_r


def closest_point_within_tolerance(query, ref_list, tol):
    dist_sq = list((ref_list - query).dot())
    dmin_sq = min(dist_sq)
    if dmin_sq > tol ** 2:
        return -1
    return dist_sq.index(dmin_sq)


class ExecutorEasyMP:
    def __init__(self,
                 cpus=21,
                 ):
        self.cpus = cpus

    def __call__(self, funcs):
        results = easy_mp.pool_map(func=wrapper_run,
                                   args=funcs,
                                   processes=self.cpus,
                                   chunksize=1,
                                   )

        return results


class ExecutorJoblib:
    def __init__(self,
                 cpus=21,
                 ):
        self.cpus = cpus

    def __call__(self, func_list):
        results = Parallel(n_jobs=self.cpus,
                           verbose=8,
                           )(delayed(func)() for func in func_list)
        return results


class FindSites:
    def __init__(self,
                 sites_grid,
                 chunk,
                 ):
        self.sites_grid = sites_grid
        self.chunk = chunk

    def __call__(self):
        return find_sites((self.sites_grid,
                           self.chunk,
                           )
                          )

    def repr(self):
        repr = OrderedDict()
        repr["len_chunk"] = len(self.chunk)
        return repr
