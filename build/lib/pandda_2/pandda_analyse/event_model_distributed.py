from __future__ import print_function

import os, sys, glob, time, re, gc

module_path = os.path.abspath(os.path.join('.'))
if module_path not in sys.path:
    sys.path.append(module_path)
b = sys.path
sys.path = [module_path] + b

import pathlib as p
import copy, warnings
from collections import OrderedDict
import time

import numpy

import joblib as jl

from libtbx.math_utils import ifloor, iceil
from scitbx.array_family import flex

from bamboo.common.status import status_bar, status_bar_2

from pandda.analyse.classes import PanddaStatMapList, MapHolderList
from pandda.analyse.functions import DatasetAligner, MapLoader, DensityStatistics, UncertaintyCalculator, wrapper_run

from dask.distributed import worker_client


def load(model):
    # Get truncated datasets
    model.dataset.sample_loader.truncate_datasets(model.dataset.datasets)
    model.dataset.datasets = model.dataset.sample_loader.truncated_datasets

    # Get cur resolution for dataset
    test_datasets = model.dataset.partition_datasets("test")
    train_datasets = model.dataset.partition_datasets("train")
    test_datasets.update(train_datasets)
    resolutions = [d.data.summary.high_res for dtag, d in test_datasets.items()]
    max_res = max(resolutions)

    # Get ref map
    model.dataset.sample_loader.get_reference(max_res)

    # Attach grid to
    model.map_maker.attach_grid(model.grid)

    # Attach grid to event table maker
    model.event_table_maker.grid = model.grid

    return model


def fit(model):

    model.fit()

    return model


def evaluate(model):

    model.evaluate()

    return model


def criticise(model):

    event_table = model.criticise()

    return event_table


class PanDDAEventModelDistributed:

    def __init__(self, statistical_model, clusterer, event_finder, bdc_calculator=None, statistics=[], dataset=None,
                 reference=None, grid=None, map_maker=None, event_table_maker=None, cpus=1, tree=None, name=None):

        self.trace = None
        self.cpus = cpus
        self.tree = tree
        self.name = name

        self.statistical_model = statistical_model
        self.clusterer = clusterer
        self.event_finder = event_finder
        self.bdc_calculator = bdc_calculator
        self.statistics = statistics
        self.map_maker = map_maker
        self.event_table_maker = event_table_maker

        self.parameters = None

        self.statistical_maps = None
        self.clusters = None
        self.events = None
        self.bdcs = None
        self.dataset = None
        self.grid = None

        self.dataset = dataset
        self.reference = reference
        self.grid = grid

    def instantiate(self, reference, tree):
        self.reference = reference
        self.tree = tree

        self.grid = self.dataset.sample_loader.get_grid(reference=self.reference)

        # Get truncated datasets
        self.dataset.sample_loader.truncate_datasets(self.dataset.datasets)
        self.dataset.datasets = self.dataset.sample_loader.truncated_datasets

        # Get cur resolution for dataset
        test_datasets = self.dataset.partition_datasets("test")
        train_datasets = self.dataset.partition_datasets("train")
        test_datasets.update(train_datasets)
        resolutions = [d.data.summary.high_res for dtag, d in test_datasets.items()]
        max_res = max(resolutions)

        # Get ref map
        self.dataset.sample_loader.get_reference(max_res)

        # Attach grid to
        self.map_maker.attach_grid(self.grid)

        # Attach grid to event table maker
        self.event_table_maker.grid = self.grid

        return self

    def clone(self, name=None, dataset=None):

        # Create a clone with new dataset and name
        model = PanDDAEventModelDistributed(self.statistical_model,
                                            self.clusterer,
                                            self.event_finder,
                                            dataset=dataset,
                                            reference=self.reference,
                                            grid=self.grid,
                                            bdc_calculator=self.bdc_calculator,
                                            statistics=self.statistics,
                                            map_maker=self.map_maker,
                                            event_table_maker=self.event_table_maker,
                                            cpus=self.cpus,
                                            tree=self.tree,
                                            name=name)

        model.tree.update({str(model.name): {"dummy": None}})

        # Attach grid to clusterer and event finder
        model.clusterer.grid = self.grid
        model.event_finder.grid = self.grid
        model.bdc_calculator.grid = self.grid

        return model

    def fit(self, samples_train, samples_test):

        self.parameters = self.statistical_model.fit(samples_train, samples_test)

        return self

    def evaluate_single(self, sample, dataset, ref_map):

        print("Making statistical map")
        statistical_map = self.statistical_model.evaluate(sample)

        print("Clustering")
        clusters = self.clusterer(dataset, statistical_map)

        print("evaluating")
        events = self.event_finder(dataset, clusters[0], clusters[1])

        print("Finding bdcs")
        bdcs = self.bdc_calculator(dataset, sample, ref_map, events)

        return statistical_map, clusters, events, bdcs

    def evaluate(self):
        dtags = set(self.dataset.partition_datasets("test").keys()
                    + self.dataset.partition_datasets("train").keys()
                    )
        truncated_datasets = self.dataset.sample_loader.truncated_datasets
        res = max([d.data.summary.high_res for dtag, d in self.dataset.datasets.items()])
        print(res)
        sample_loaders = {dtag: lambda d: self.dataset.sample_loader.get_sample(res, d)
                          for dtag
                          in dtags}
        gc.collect()
        results = jl.Parallel(n_jobs=int(self.cpus),
                              verbose=10)(jl.delayed(self.evaluate_single)(sample_loaders[dtag],
                                                                           truncated_datasets[dtag],
                                                                           self.dataset.sample_loader.ref_map)
                                          for dtag
                                          in dtags)
        # with worker_client() as client:
        #     results = client.map(self.evaluate_single, [(sample_loaders[dtag],
        #                                                  truncated_datasets[dtag],
        #                                                  self.dataset.sample_loader.ref_map)
        #                                                 for dtag
        #                                                 in dtags])


        # self.statistical_maps = {dtag: res[0]
        #                          for dtag, res
        #                          in zip(dtags, results)}
        self.clusters = {dtag: res[1]
                         for dtag, res
                         in zip(dtags, results)}
        self.events = {dtag: res[2]
                       for dtag, res
                       in zip(dtags, results)}

        self.bdcs = {dtag: res[3]
                       for dtag, res
                       in zip(dtags, results)}

        return self

    def criticise_single(self, sample_loader, truncated_dataset, ref_map, events, bdcs, tree):

        dataset_path = p.Path(tree(("processed_datasets", truncated_dataset.tag))[0])

        self.map_maker.process_single(sample_loader, truncated_dataset, ref_map, events, bdcs, dataset_path)

    def criticise_all(self, tree):
        # dtags = set(self.dataset.partition_datasets("test").keys()
        #             + self.dataset.partition_datasets("train").keys()
        #             )
        #
        # res = max([d.data.summary.high_res for dtag, d in self.dataset.datasets.items()])
        # sample_loaders = {dtag: lambda d: self.dataset.sample_loader.get_sample(res, d)
        #                   for dtag
        #                   in dtags}
        #
        # self.map_maker.statistical_model = self.statistical_model
        #
        # gc.collect()
        # # jl.Parallel(n_jobs=int(self.cpus),
        # #             verbose=10)(jl.delayed(self.criticise_single)(sample_loaders[dtag],
        # #                                                           self.dataset.sample_loader.truncated_datasets[dtag],
        # #                                                           self.dataset.sample_loader.ref_map,
        # #                                                           self.events[dtag],
        # #                                                           self.bdcs[dtag],
        # #                                                           self.tree)
        # #                         for dtag
        # #                         in dtags)

        # Produce maps that are shared by iteration
        dir_path = p.Path(tree([str(self.name)])[0])
        dir_path_string = str(dir_path)
        self.map_maker.process_shell(self.dataset.sample_loader.reference,
                                     self.dataset.sample_loader.ref_map,
                                     dir_path_string)

        # Produce the event table
        event_table_path = dir_path / "event_table.csv"

        event_table = self.event_table_maker(self.dataset.partition_datasets("test"),
                                             self.events,
                                             event_table_path)

        return event_table

    def log(self):
        log = OrderedDict()
        return log


class PanDDANormalModel:

    def __init__(self, method="adjusted+uncertainty", cpus=1):

        self.method = method
        self.cpus = cpus

        self.mu = 0
        self.sigma_uncertainty = {}
        self.sigma_adjusted = 0

        # TODO: legacy requirement
        self.statistical_maps = PanddaStatMapList()

    def fit(self, samples_train, samples_test):

        # TODO: move into fit
        map_data_size = 0
        for dtag, sample in samples_train.items():
            map_data_size = sample.data.size()
            break

        characterisation_maps = [sample for dtag, sample in samples_train.items()]

        analysis_maps = [sample for dtag, sample in samples_test.items()]

        self.mu = self.fit_mu(characterisation_maps, map_data_size)

        self.sigma_uncertainty = self.fit_sigma_uncertainty(analysis_maps, map_data_size, cpus=self.cpus)

        self.sigma_adjusted = self.fit_sigma_adjusted(analysis_maps, self.sigma_uncertainty, map_data_size,
                                                      cpus=self.cpus)

    def fit_mu(self, dataset_maps, map_data_size):
        """Calculate the average map from all of the different observations"""
        print("\t### Fitting mu!")

        # Extract the maps to be used for averaging
        if len(dataset_maps) == 1:

            # Extract the map from the list
            m = dataset_maps[0]
            # Mean and median are simply the map value -- copy directly to the statistical maps
            mean_map_vals = medn_map_vals = numpy.array(m.data)

        else:

            # Chunk the points into groups - Compromise between cpu time and memory usage - ~200 dataset -> chunksize of 5000
            chunk_size = 500 * iceil(1000.0 / len(dataset_maps))
            chunk_idxs = [i for i in range(0, map_data_size, chunk_size)]
            num_chunks = len(chunk_idxs)

            t1 = time.time()

            mean_map_vals = numpy.zeros(map_data_size)
            medn_map_vals = numpy.zeros(map_data_size)

            for i_chunk, chunk_start in enumerate(chunk_idxs):
                status_bar_2(n=i_chunk, n_max=num_chunks)

                tmp_map_vals = numpy.array([m.data[chunk_start:chunk_start + chunk_size] for m in dataset_maps])

                # Check that the output values are the expected dimensions
                if i_chunk + 1 < num_chunks:
                    assert len(tmp_map_vals) == len(dataset_maps)
                    assert len(tmp_map_vals.T) == chunk_size

                tmp_map_means = numpy.mean(tmp_map_vals, axis=0)
                mean_map_vals[chunk_start:chunk_start + chunk_size] = tmp_map_means
                tmp_map_medns = numpy.median(tmp_map_vals, axis=0)
                medn_map_vals[chunk_start:chunk_start + chunk_size] = tmp_map_medns

            status_bar_2(n=num_chunks, n_max=num_chunks)

            t2 = time.time()

        mu = m.new_from_template(map_data=flex.double(mean_map_vals.flatten()),
                                                             sparse=m.is_sparse())

        return mu

    def fit_sigma_uncertainty(self, analysis_maps, map_data_size, masked_idxs=None, mask_name=None, q_cut=1.5, cpus=1):
        """Calculate the uncertainty in each of the different maps"""

        print("\t### Fitting sigma_uncertainty!")

        if masked_idxs is None:
            masked_idxs = flex.size_t(range(0, map_data_size))
        else:
            assert max(masked_idxs) < map_data_size, 'masked_idxs out of range of map'
            masked_idxs = flex.size_t(masked_idxs)

        # Extract masked map values from the average map... and sort them
        comp_vals = self.mu.data.select(masked_idxs)

        arg_list = []

        # for i_m, m in enumerate(self.dataset_maps.mask(mask_name=mask_name)):
        for i_m, m in enumerate(analysis_maps):

            if m.meta.map_uncertainty is not None:
                arg_list.append(None)
                continue

            u = UncertaintyCalculator(query_values=m.data.select(masked_idxs),
                                      ref_values=comp_vals)
            arg_list.append(u)

        t1 = time.time()
        num_to_process = len(arg_list) - arg_list.count(None)
        print('1' + ''.join(['{:<5}'.format(i) for i in range(0, num_to_process + 5, 5)])[2:])
        print(' ' * num_to_process + '|\r', end='')
        sys.stdout.flush()
        # TODO: use joblib instead
        # map_uncertainties = easy_mp.pool_map(func=wrapper_run, args=arg_list, processes=cpus, chunksize=1)
        map_uncertainties = jl.Parallel(n_jobs=self.cpus,
                                        verbose=5)(jl.delayed(wrapper_run)(arg)
                                                    for arg
                                                    in arg_list)
        # with worker_client() as client:
        #     map_uncertainties_futures = client.map(wrapper_run, arg_list)
        #     map_uncertainties = client.gather(map_uncertainties_futures)
        print('|')

        for i_m, m in enumerate(analysis_maps):

            map_unc = map_uncertainties[i_m]

            if m.meta.map_uncertainty is not None:
                assert map_unc is None
            else:
                # TODO: remove this print
                # print("Adding map uncertainty for {}".format(m.meta.tag))
                assert map_unc is not None
                m.meta.map_uncertainty = map_unc
                # TODO: Not sure why this is breaking - probably to do with futures print

        # return [m.meta.map_uncertainty for m in self.dataset_maps.mask(mask_name=mask_name)]
        return {m.meta.tag: m.meta.map_uncertainty for m in analysis_maps}

    def fit_sigma_adjusted(self, analysis_maps, uncertainties, map_data_size, cpus=1):

        print("\t### Fitting sigma_adjusted!")

        uncertainties_ordered = [uncertainties[edmap.meta.tag]
                                 for edmap
                                 in analysis_maps]

        self.calculate_statistical_maps(analysis_maps, uncertainties_ordered, map_data_size, cpus=self.cpus)

        return self.statistical_maps.sadj_map

    def calculate_statistical_maps(self, dataset_maps, uncertainties, map_data_size, mask_name=None, ignore_warnings=True, cpus=1):
        """Take the sampled maps and calculate statistics for each grid point across the datasets"""

        # Extract the maps to be used for averaging

        if len(dataset_maps) == 1:

            self._set_statistical_maps_from_array(template_map=self.mu,
                                                  map_array=numpy.zeros((map_data_size, 5)))

            return self.statistical_maps
        else:

            # Create statistics objects for each grid point
            if ignore_warnings:
                warnings.simplefilter('ignore', category=RuntimeWarning)

            # Extract the map uncertainties
            # uncertainties = [m.meta.map_uncertainty for m in dataset_maps]
            assert uncertainties.count(None) == 0, 'some maps have not got associated uncertainties'

            # Chunk the points into groups - Compromise between cpu time and memory usage - 1000 per cpu at 50 datasets
            chunk_size = iceil(1000.0 * cpus * 50.0 / len(dataset_maps))
            chunk_idxs = [i for i in range(0, map_data_size, chunk_size)]
            num_chunks = len(chunk_idxs)

            # Second level of iteration - split the first chunk level between the cpus
            chunk_size_2 = iceil(1.0 * chunk_size / cpus)
            chunk_idxs_2 = [i for i in range(0, chunk_size, chunk_size_2)]
            num_chunks_2 = len(chunk_idxs_2)

            t1 = time.time()

            # Output array of the 5 statistics for each map point
            point_statistics = numpy.zeros((map_data_size, 5))

            tot = 0
            for i_chunk, chunk_start in enumerate(chunk_idxs):
                status_bar_2(n=i_chunk, n_max=num_chunks)

                # Argument list for multiprocessing
                arg_list = []

                # Loop through the secondary chunks and send for multi-core processing
                for i_chunk_2, chunk_start_2 in enumerate(chunk_idxs_2):

                    # Lower limit - always the beginning of the chunk
                    l1 = chunk_start + chunk_start_2
                    # Upper limit - full chunk size, limited by the larger chunk size, or by map size
                    l2 = min(chunk_start + chunk_start_2 + chunk_size_2, chunk_start + chunk_size,
                             map_data_size)

                    if l1 >= l2:
                        continue

                    # Extract map values from the maps
                    map_vals = [m.data[l1:l2] for m in dataset_maps]
                    # Want to iterate over grid points not datasets
                    map_vals = numpy.transpose(map_vals)
                    assert map_vals.shape[1] == len(dataset_maps)

                    # Create DensityStatistics object for analysis of the density variation
                    arg_list.append(DensityStatistics(observations_array=map_vals, uncertainties=uncertainties))

                if not arg_list: continue

                # Calculate the statistis of the grid points
                # TODO: use joblib instead
                # tmp_point_statistics = easy_mp.pool_map(func=wrapper_run, args=arg_list, processes=cpus)
                tmp_point_statistics = jl.Parallel(n_jobs=self.cpus)(jl.delayed(wrapper_run)(arg)
                                                                     for arg
                                                                     in arg_list)
                # with worker_client() as client:
                #     tmp_point_statistics_futures = client.map(wrapper_run, arg_list)
                #     tmp_point_statistics = client.gather(tmp_point_statistics_futures)

                # Put values into the output array
                offset = 0
                for point_vals in tmp_point_statistics:
                    assert point_vals.shape[1] == 5
                    l1 = chunk_start + offset
                    l2 = l1 + point_vals.shape[0]
                    if not (point_statistics[l1:l2, :] == 0.0).all():
                        print('Overwriting data?!')
                        print(point_statistics[l1 - 10:l2 + 10, :])
                        assert point_statistics[l1:l2, :].shape == point_vals.shape, '{} != {}'.format(
                            point_statistics[l1:l2, :].shape, point_vals.shape)
                    point_statistics[l1:l2, :] = point_vals
                    offset += point_vals.shape[0]
                tot += offset

            status_bar_2(n=num_chunks, n_max=num_chunks)

            # Check that we've calculated the right number of things
            assert tot == map_data_size, 'tot {}, map size {}'.format(tot, map_data_size)

            t2 = time.time()

        self._set_statistical_maps_from_array(template_map=self.mu,
                                              map_array=point_statistics,
                                              map_data_size=map_data_size)

    def _set_statistical_maps_from_array(self, template_map, map_array, map_data_size):
        """Set the five non-average-based statistical maps from an array"""

        assert map_array.shape == (map_data_size, 5)

        # Create the other statistical maps
        self.statistical_maps.stds_map = template_map.new_from_template(
            map_data=flex.double(map_array[:, 0].tolist()), sparse=template_map.is_sparse())
        self.statistical_maps.sadj_map = template_map.new_from_template(
            map_data=flex.double(map_array[:, 1].tolist()), sparse=template_map.is_sparse())
        self.statistical_maps.skew_map = template_map.new_from_template(
            map_data=flex.double(map_array[:, 2].tolist()), sparse=template_map.is_sparse())
        self.statistical_maps.kurt_map = template_map.new_from_template(
            map_data=flex.double(map_array[:, 3].tolist()), sparse=template_map.is_sparse())
        self.statistical_maps.bimo_map = template_map.new_from_template(
            map_data=flex.double(map_array[:, 4].tolist()), sparse=template_map.is_sparse())

    def evaluate(self, m):
        """Calculate the z-map relative to the mean and std map"""

        assert self.method in ['none', 'adjusted', 'uncertainty', 'adjusted+uncertainty']
        uncertainty = self.sigma_uncertainty[m.meta.tag]

        # Check that a value has been found/supplied
        if 'uncertainty' in self.method:
            assert uncertainty is not None

        # Extract maps in the right sparseness
        is_sparse = m.is_sparse()
        # Extract mean values (for subtraction)
        comp_vals = self.mu.get_map_data(sparse=is_sparse)

        # Extract the normalisation values (for division)
        if self.method == 'none':
            norm_vals = 1.0
        #        elif method == 'naive':
        #            norm_vals = self.statistical_maps.stds_map.get_map_data(sparse=is_sparse)
        elif self.method == 'adjusted':
            norm_vals = self.sigma_adjusted.get_map_data(sparse=is_sparse)
        elif self.method == 'uncertainty':
            norm_vals = uncertainty
        elif self.method == 'adjusted+uncertainty':
            norm_vals = flex.sqrt(
                self.sigma_adjusted.get_map_data(sparse=is_sparse) ** 2 + uncertainty ** 2)
        else:
            raise Exception('method not found: {!s}'.format(self.method))

        return (m - comp_vals) * (1.0 / norm_vals)
