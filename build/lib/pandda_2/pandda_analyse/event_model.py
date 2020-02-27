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
import scipy.cluster
from hdbscan import HDBSCAN

import joblib as jl

from libtbx.math_utils import ifloor, iceil
from scitbx.array_family import flex

from bamboo.common.status import status_bar, status_bar_2
from bamboo.stats.cluster import find_connected_groups, generate_group_idxs


from giant.grid.utils import idx_to_grid
from giant.structure.select import protein, non_water, find_nearest_atoms
from giant.xray.symmetry import find_symmetry_equivalent_groups
from giant.grid.masks import GridMask
from giant.xray.maps.bdc import calculate_varying_bdc_correlations, calculate_maximum_series_discrepancy, calculate_bdc_subtracted_map

from pandda.analyse.classes import PanddaStatMapList, MapHolderList
from pandda.analyse.functions import DatasetAligner, MapLoader, DensityStatistics, UncertaintyCalculator, wrapper_run

from pandda.analyse.events import PointCluster, Event

from dask.distributed import worker_client


def analyse_dataset_events(dataset,
                           dataset_map,
                           ref_map,
                           events,
                           grid,
                           max_bdc=1.0,
                           min_bdc=0.0,
                           increment=0.01,
                           output_multiplier=1.0,
                           ):

    # ============================================================================>
    # Extract the map data in non-sparse format
    # ============================================================================>
    dset_map_data = dataset_map.get_map_data(sparse=False)
    ref_map_data = ref_map.get_map_data(sparse=False)
    # ============================================================================>
    # Unpack cluster
    # ============================================================================>
    event_stats = OrderedDict()
    for event in events[2]:
        # ============================================================================>
        # Estimate the background correction of the detected feature
        # ============================================================================>
        # Extract sites for this cluster and estimate the background correction for the event

        # Generate custom grid mask for this dataset
        event_mask = GridMask(parent=grid,
                              sites_cart=grid.grid2cart(event.cluster.points,
                                                             origin_shift=True),
                              max_dist=2.0,
                              min_dist=0.0,
                              )

        # Select masks to define regions for bdc calculation
        exp_event_idxs = flex.size_t(event_mask.outer_mask_indices())
        reference_idxs = flex.size_t(grid.global_mask().inner_mask_indices())
        # ============================================================================>
        # Generate BDC-estimation curve and estimate BDC
        # ============================================================================>
        event_remains, event_corrs, global_corrs = calculate_varying_bdc_correlations(
            ref_map_data=ref_map_data,
            query_map_data=dset_map_data,
            feature_idxs=exp_event_idxs,
            reference_idxs=reference_idxs,
            min_remain=1.0 - max_bdc,
            max_remain=1.0 - min_bdc,
            bdc_increment=increment,
            verbose=True)
        event_remain_est = calculate_maximum_series_discrepancy(
            labels=event_remains,
            series_1=global_corrs,
            series_2=event_corrs)

        event_remain_est = min(event_remain_est * output_multiplier,
                               1.0 - min_bdc)
        # ============================================================================>
        # Calculate the map correlations at the selected BDC
        # ============================================================================>
        event_map_data = calculate_bdc_subtracted_map(
            ref_map_data=ref_map_data,
            query_map_data=dset_map_data,
            bdc=1.0 - event_remain_est)
        global_corr = \
            numpy.corrcoef(event_map_data.select(reference_idxs), ref_map_data.select(reference_idxs))[
                0, 1]
        local_corr = numpy.corrcoef(event_map_data.select(exp_event_idxs), ref_map_data.select(exp_event_idxs))[
            0, 1]
        # ============================================================================>
        # Update event parameters
        # ============================================================================>
        event.info.estimated_pseudo_occupancy = event_remain_est
        event.info.estimated_bdc = 1.0 - event_remain_est
        event.info.global_correlation = global_corr
        event.info.local_correlation = local_corr

        # ============================================================================>
        # Find the nearest atom to the event
        # ============================================================================>
        # TODO: restore this?
        atm = find_nearest_atoms(atoms=list(protein(dataset.model.hierarchy).atoms_with_labels()),
                                 query=dataset.model.alignment.ref2nat(
                                     grid.grid2cart(sites_grid=[map(int, event.cluster.centroid)],
                                                         origin_shift=True)))[0]

        event_stats[event.id] = OrderedDict()
        event_stats[event.id]["estimated_pseudo_occupancy"] = event_remain_est
        event_stats[event.id]["estimated_bdc"] = 1.0 - event_remain_est
        event_stats[event.id]["global_corr"] = global_corr
        event_stats[event.id]["local_corr"] = global_corr

    return event_stats


class PanDDAEventModel:

    def __init__(self, statistical_model, clusterer, event_finder, bdc_calculator=None, statistics=[],
                 map_maker=None, event_table_maker=None, cpus=14, tree=None, name=None):

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

    def __call__(self, dataset, reference, name=None):
        self.name = name
        model = PanDDAEventModel(self.statistical_model,
                                 self.clusterer,
                                 self.event_finder,
                                 bdc_calculator=self.bdc_calculator,
                                 statistics=self.statistics,
                                 map_maker=self.map_maker,
                                 event_table_maker=self.event_table_maker,
                                 cpus=self.cpus,
                                 tree=self.tree,
                                 name=name)

        model.tree.update({str(model.name): {"dummy": None}})

        model.dataset = dataset
        model.reference = reference

        return model

    def process(self):
        with self as model:
            # Fit model
            print("Fiting model")
            model.fit()

            # Drop samples before multiprocessing
            model.drop_samples()

            # Set up individual sample loaders
            dtags = self.dataset.partition_samples("test").keys()
            truncated_datasets = self.dataset.sample_loader.truncated_datasets
            sample_loaders = {dtag: lambda d: self.dataset.sample_loader.get_sample(d)
                              for dtag
                              in dtags}

            # Evaluate model for eahc of the datasets
            print("Evaluating events in parallel")
            results = jl.Parallel(n_jobs=self.cpus,
                                  verbose=5)(jl.delayed(self.evaluate_single)(sample_loaders[dtag],
                                                                              truncated_datasets[dtag])
                                             for dtag
                                             in dtags)

            self.statistical_maps = {dtag: res[0]
                                     for dtag, res
                                     in zip(dtags, results)}
            self.clusters = {dtag: res[1]
                             for dtag, res
                             in zip(dtags, results)}
            self.events = {dtag: res[2]
                           for dtag, res
                           in zip(dtags, results)}
            self.bdcs = {dtag: res[3]
                         for dtag, res
                         in zip(dtags, results)}

            # Drop stat maps for memory
            del self.statistical_maps

            # Criticise each dataset
            print("Criticising events in parallel")
            jl.Parallel(n_jobs=self.cpus,
                        verbose=5)(jl.delayed(self.criticise_single)(sample_loaders[dtag],
                                                                     truncated_datasets[dtag],
                                                                     self.events[dtag],
                                                                     self.bdcs[dtag],
                                                                     self.tree)
                                   for dtag
                                   in dtags)

            # Criticise whole model

    def __enter__(self):
        # TODO: move grid getter to its own class
        # Fit the grid the sample loader will use
        self.grid = self.dataset.sample_loader.get_grid(reference=self.reference)

        # Attach grid to clusterer and event finder
        self.clusterer.grid = self.dataset.sample_loader.grid
        self.event_finder.grid = self.dataset.sample_loader.grid
        self.bdc_calculator.grid = self.dataset.sample_loader.grid

        # Get cur resolution for dataset
        test_datasets = self.dataset.partition_datasets("test")
        train_datasets = self.dataset.partition_datasets("train")
        test_datasets.update(train_datasets)
        resolutions = [d.data.summary.high_res for dtag, d in test_datasets.items()]
        max_res = max(resolutions)
        print("Max res of dataset is {}".format(max_res))
        print("Min res of dataset is {}".format(min(resolutions)))

        # Get truncated datasets
        self.dataset.sample_loader.truncate_datasets(self.dataset.datasets)

        # Get ref map
        self.dataset.sample_loader.get_reference(max_res)

        # Attach grid to
        self.map_maker.attach_grid(self.grid)

        # Attach grid to event table maker
        self.event_table_maker.grid = self.grid

        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.samples = None

    def fit(self):
        resolutions_test = max([d.data.summary.high_res for dtag, d
                                in self.dataset.partition_datasets("test").items()])
        resolutions_train = max([d.data.summary.high_res for dtag, d
                                 in self.dataset.partition_datasets("train").items()])
        max_res = max(resolutions_test, resolutions_train)

        # Load train samples for context
        self.dataset.load_samples_from_partition(max_res, "train")

        samples_train = self.dataset.partition_samples("train")
        samples_test = self.dataset.partition_samples("test")

        self.parameters = self.statistical_model.fit(samples_train, samples_test)

        self.dataset.dump_samples()

    def evaluate_single(self, sample_loader, dataset, ref_map):

        print("Getting sample")
        sample = sample_loader(dataset)

        print("Making statistical map")
        statistical_map = self.statistical_model.evaluate(sample)

        print("Clustering")
        clusters = self.clusterer(dataset, statistical_map)

        print("evaluating")
        events = self.event_finder(dataset, clusters[0], clusters[1])

        print("Finding bdcs")
        bdcs = self.bdc_calculator(dataset, sample, ref_map, events)

        return statistical_map, clusters, events, bdcs

    def evaluate_all(self):
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

    # def evaluate_parallel(self):
    #
    #     dtags = self.dataset.partition_samples("test").keys()
    #     truncated_datasets = self.dataset.sample_loader.truncated_datasets
    #     sample_loaders = {dtag: lambda d: self.dataset.sample_loader.get_sample(d)
    #                       for dtag
    #                       in dtags}
    #
    #     results = jl.Parallel(n_jobs=self.cpus,
    #                           verbose=5)(jl.delayed(self.evaluate)(sample_loaders[dtag],
    #                                                                truncated_datasets[dtag])
    #                                      for dtag
    #                                      in dtags)
    #
    #     self.statistical_maps = {dtag: res[0]
    #                              for dtag, res
    #                              in zip(dtags, results)}
    #     self.clusters = {dtag: res[1]
    #                      for dtag, res
    #                      in zip(dtags, results)}
    #     self.events = {dtag: res[2]
    #                    for dtag, res
    #                    in zip(dtags, results)}

    # def evaluate_seriel(self):
    #
    #     samples_test = self.dataset.partition_samples("test")
    #
    #     statistical_maps = {dtag: self.statistical_model.evaluate(sample)
    #                         for dtag, sample
    #                         in samples_test.items()}
    #
    #     print("Clustering")
    #     clusters = {dtag: self.clusterer(self.dataset.datasets[dtag], statistical_map)
    #                 for dtag, statistical_map
    #                 in statistical_maps.items()}
    #
    #     print("evaluating")
    #     events = {self.event_finder(self.dataset.datasets[dtag], cluster[0], cluster[1])
    #               for dtag, cluster
    #               in clusters.items()}

    def criticise_single(self, sample_loader, truncated_dataset, ref_map, events, bdcs, tree):

        dataset_path = p.Path(tree(("processed_datasets", truncated_dataset.tag))[0])

        self.map_maker.process_single(sample_loader, truncated_dataset, ref_map, events, bdcs, dataset_path)

    def criticise_all(self, tree):
        dtags = set(self.dataset.partition_datasets("test").keys()
                    + self.dataset.partition_datasets("train").keys()
                    )

        res = max([d.data.summary.high_res for dtag, d in self.dataset.datasets.items()])
        sample_loaders = {dtag: lambda d: self.dataset.sample_loader.get_sample(res, d)
                          for dtag
                          in dtags}

        self.map_maker.statistical_model = self.statistical_model

        gc.collect()
        # jl.Parallel(n_jobs=int(self.cpus),
        #             verbose=10)(jl.delayed(self.criticise_single)(sample_loaders[dtag],
        #                                                           self.dataset.sample_loader.truncated_datasets[dtag],
        #                                                           self.dataset.sample_loader.ref_map,
        #                                                           self.events[dtag],
        #                                                           self.bdcs[dtag],
        #                                                           self.tree)
        #                         for dtag
        #                         in dtags)

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

    def __init__(self, method="adjusted+uncertainty", cpus=1, multiprocess="dask"):

        self.method = method
        self.cpus = cpus

        self.mu = 0
        self.sigma_uncertainty = {}
        self.sigma_adjusted = 0
        self.multiprocess = multiprocess

        # TODO: legacy requirement
        self.statistical_maps = PanddaStatMapList()

    def fit(self, samples_train, samples_test):

        samples_train = {sample.meta.tag: sample for sample in samples_train}
        samples_test = {sample.meta.tag: sample for sample in samples_test}

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

        return self

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
        if self.multiprocess == "dask":
            # dsk = {}
            # for i, arg in enumerate(numpy.array_split(numpy.array(arg_list), 10)):
            #     dsk["arg_{}".format(i)] = arg
            #     dsk["uncertainty_{}".format(i)] = (map,
            #                                        wrapper_run,
            #                                        "arg_{}".format(i))
            #
            # dsk["map_uncertainties"] = (list,
            #                             ["uncertainty_{}".format(i)
            #                              for i, arg
            #                              in enumerate(numpy.array_split(numpy.array(arg_list), 10))])
            #
            # print("Computing!")
            # with worker_client() as client:
            #     map_uncertainties = client.get(dsk, "map_uncertainties")
            #
            # map_uncertainties = sum(map_uncertainties, [])

            from dask import delayed
            map_uncertainties = map(delayed(wrapper_run), [arg for arg in arg_list])
            # print(map_uncertainties)
            map_uncertainties = delayed(list)(map_uncertainties).compute()
            # print(map_uncertainties)

            # dsk = {}
            # for i, arg in enumerate(arg_list):
            #     dsk["arg_{}".format(i)] = arg
            #     dsk["uncertainty_{}".format(i)] = (wrapper_run, "arg_{}".format(i))
            #
            # dsk["map_uncertainties"] = (list,
            #                             ["uncertainty_{}".format(i) for i, arg in enumerate(arg_list)])
            #
            # with worker_client() as client:
            #     map_uncertainties = client.get(dsk, "map_uncertainties")

            # map_uncertainties = [wrapper_run(arg) for arg in arg_list]

        # elif self.multiprocess == "dask_parallel":
        #     with worker_client() as client:
        #         map_uncertainties_futures = client.map(wrapper_run, arg_list)
        #         map_uncertainties = client.gather(map_uncertainties_futures)
        #
        # else:
        #     map_uncertainties = jl.Parallel(n_jobs=self.cpus,
        #                                     verbose=5)(jl.delayed(wrapper_run)(arg)
        #                                                 for arg
        #                                                 in arg_list)
        map_uncertainties = jl.Parallel(n_jobs=self.cpus,
                                        verbose=5)(jl.delayed(wrapper_run)(arg)
                                                   for arg
                                                   in arg_list)

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

############################################################################
            tot = 0
            chunk_list = []
            for i_chunk, chunk_start in enumerate(chunk_idxs):
                # status_bar_2(n=i_chunk, n_max=num_chunks)

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

                chunk_list.append(arg_list)
                #
                # # Calculate the statistis of the grid points
                # # TODO: use joblib instead
                # # tmp_point_statistics = easy_mp.pool_map(func=wrapper_run, args=arg_list, processes=cpus)
                # if self.multiprocess == "dask":
                #     # tmp_point_statistics = map(wrapper_run, arg_list)
                #     dsk = {}
                #     for i, arg in enumerate(numpy.array_split(numpy.array(arg_list), 10)):
                #         dsk["arg_{}".format(i)] = arg
                #         dsk["tmp_point_statistics_{}".format(i)] = (map,
                #                                                     wrapper_run,
                #                                                     "arg_{}".format(i))
                #
                #     dsk["tmp_point_statistics"] = (list,
                #                                    ["tmp_point_statistics_{}".format(i)
                #                                     for i, arg
                #                                     in enumerate(numpy.array_split(numpy.array(arg_list), 10))])
                #
                #     with worker_client() as client:
                #         tmp_point_statistics = client.get(dsk, "tmp_point_statistics")
                #
                #     tmp_point_statistics = sum(tmp_point_statistics, [])
                #
                #     # dsk = {}
                #     # for i, arg in enumerate(arg_list):
                #     #     dsk["arg_{}".format(i)] = arg
                #     #     dsk["tmp_point_statistics_{}".format(i)] = (wrapper_run,
                #     #                                                 "arg_{}".format(i))
                #     #
                #     # dsk["tmp_point_statistics"] = (list,
                #     #                             ["tmp_point_statistics_{}".format(i) for i, arg in enumerate(arg_list)])
                #     #
                #     # with worker_client() as client:
                #     #     tmp_point_statistics = client.get(dsk, "tmp_point_statistics")
                #
                #     # tmp_point_statistics = [wrapper_run(arg) for arg in arg_list]
                #
                # elif self.multiprocess == "dask_parallel":
                #     with worker_client() as client:
                #         tmp_point_statistics_futures = client.map(wrapper_run, arg_list)
                #         tmp_point_statistics = client.gather(tmp_point_statistics_futures)
                # else:
                #     tmp_point_statistics = jl.Parallel(n_jobs=self.cpus)(jl.delayed(wrapper_run)(arg)
                #                                                          for arg
                #                                                          in arg_list)

            def chunk(l, n):
                for j in xrange(0, len(l), n):
                    yield l[j:j+n]

            def process_chunk(chunk):
                tps = map(wrapper_run, chunk)
                return tps

            def process_chunks(chunks_list):

                results = []
                for chunk in chunk_list:
                    chunk_results = map(wrapper_run, chunk)
                    results.append(chunk_results)

                results = sum(results, [])

                return results

            # dsk = {}
            #
            # print("### Getting tmp point statistics")
            #
            chunks_blocks = list(chunk(chunk_list, 10))
            #
            # print(len(chunks_blocks))
            # print([len(x) for x in chunks_blocks])
            #
            # for i, arg in enumerate(chunks_blocks):
            #     dsk["arg_{}".format(i)] = arg
            #     dsk["tmp_point_statistics_{}".format(i)] = (process_chunks,
            #                                                 "arg_{}".format(i))
            #
            # dsk["tmp_point_statistics"] = (list,
            #                                ["tmp_point_statistics_{}".format(i)
            #                                 for i, arg
            #                                 in enumerate(chunks_blocks)])
            #
            # with worker_client() as client:
            #     tmp_point_statistics = client.get(dsk, "tmp_point_statistics")
            #
            # tmp_point_statistics = sum(tmp_point_statistics, [])

            from dask import delayed
            tmp_point_statistics = map(delayed(process_chunk), [chunk_block for chunk_block in chunk_list])
            # print(tmp_point_statistics)
            tmp_point_statistics = delayed(list)(tmp_point_statistics).compute()
            # print(tmp_point_statistics)

            for i_chunk, chunk_start in enumerate(chunk_idxs):

                # Put values into the output array
                offset = 0
                for point_vals in tmp_point_statistics[i_chunk]:
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


class PanDDADefaultCluster:

    def __init__(self, grid_clustering_cutoff, negative_values, cluster_method, contour_level, outer_mask,
                 inner_mask_symmetry, grid=None, dataset_mask=None):
        self.grid_clustering_cutoff = grid_clustering_cutoff
        self.negative_values = negative_values
        self.cluster_method = cluster_method
        self.contour_level = contour_level

        self.outer_mask = outer_mask
        self.inner_mask_symmetry = inner_mask_symmetry

    def __call__(self, dataset, z_map, grid):

        z_map_data = z_map.get_map_data(sparse=False)

        # ============================================================================>
        # Extract the global mask object from the grid
        # ============================================================================>
        dset_total_temp = grid.global_mask().total_mask_binary().copy()

        # ============================================================================>
        # Generate symmetry masks for this dataset
        # ============================================================================>
        # Generate symmetry contacts for this dataset and align to reference frame
        dataset_sym_copies = dataset.model.crystal_contacts(distance_cutoff=self.outer_mask + 5,
                                                            combine_copies=True)
        dataset_sym_copies.atoms().set_xyz(
            dataset.model.alignment.nat2ref(dataset_sym_copies.atoms().extract_xyz()))

        # Extract protein atoms from the symmetry copies
        dataset_sym_sites_cart = non_water(dataset_sym_copies).atoms().extract_xyz()
        # Generate symmetry contacts grid mask
        dataset_mask = GridMask(parent=grid,
                                sites_cart=dataset_sym_sites_cart,
                                max_dist=self.outer_mask,
                                min_dist=self.inner_mask_symmetry)

        # Combine with the total mask to generate custom mask for this dataset
        dset_total_temp.put(dataset_mask.inner_mask_indices(), 0)
        point_mask_idx = numpy.where(dset_total_temp)[0]

        # Truncate map
        above_val, above_idx, above_gps = self.truncate_map(z_map_data,
                                                            point_mask_idx,
                                                            self.negative_values,
                                                            self.contour_level)

        # Check valid
        # No Cluster points found
        if len(above_val) == 0:
            return 0, []
        # One Cluster point found
        elif len(above_val) == 1:
            return 1, [(above_gps, above_val)]

        # Cluster
        cluter_time_start = time.time()
        cluster_ids = self.cluster_hierarcical(above_gps, self.grid_clustering_cutoff)
        cluter_time_finish = time.time()

        # Get number of clusters
        num_clusters = self.get_num_clusters(cluster_ids)

        # Group clusters
        z_clusters = self.group_clusters(cluster_ids,
                                         num_clusters,
                                         above_gps,
                                         above_val)

        return num_clusters, z_clusters

    def mask_map(self, z_map_data, point_mask_idx):
        # Select these values from the map
        point_mask_idx = flex.size_t(point_mask_idx)
        point_mask_val = z_map_data.select(point_mask_idx)
        return point_mask_val

    def truncate_map(self, z_map_data, point_mask_idx, negative_values=False, contour_level=0):

        # Mask map
        z_map_masked = self.mask_map(z_map_data, point_mask_idx)

        # truncate map
        z_truncated_map_truncated = self.get_truncated_map_mask(z_map_masked,
                                                                negative_values,
                                                                contour_level)

        # Extract these points
        above_val, above_idx, above_gps = self.extract_points(z_map_masked,
                                                              point_mask_idx,
                                                              z_truncated_map_truncated,
                                                              z_map_data)

        return above_val, above_idx, above_gps

    def get_truncated_map_mask(self, point_mask_val, negative_values=False, contour_level=0):
        # Find values above cutoff
        if negative_values:
            above_idx = (point_mask_val >= contour_level).iselection()
            below_idx = (point_mask_val <= -1.0 * contour_level).iselection()
            sel_idx = above_idx.concatenate(below_idx)
        else:
            sel_idx = (point_mask_val >= contour_level).iselection()

        return sel_idx

    def extract_points(self, point_mask_val, point_mask_idx, sel_idx, z_map_data):
        point_mask_idx = flex.size_t(point_mask_idx)
        # Extract values and grid points for these sites
        above_val = point_mask_val.select(sel_idx)
        above_idx = point_mask_idx.select(sel_idx)
        above_gps = flex.vec3_double([idx_to_grid(i, grid_size=z_map_data.all()) for i in above_idx])
        # above_len = len(above_val)

        return above_val, above_idx, above_gps

    def cluster_hdbscan(self, above_gps):
        sample_by_feature = above_gps.to_numpy()
        print(sample_by_feature.shape)

        clusterer = HDBSCAN()
        clusterer.fit(sample_by_feature)
        cluster_ids = list(clusterer.labels_)

        return cluster_ids

    def cluster_hierarcical(self, above_gps, grid_clustering_cutoff):
        cluster_ids = scipy.cluster.hierarchy.fclusterdata(X=above_gps,
                                                           t=grid_clustering_cutoff,
                                                           criterion='distance',
                                                           metric='euclidean',
                                                           method='single')
        cluster_ids = list(cluster_ids)
        return cluster_ids

    def get_num_clusters(self, cluster_ids):
        # Get the number of clusters
        cluster_id_array = numpy.array(cluster_ids)
        unique_cluster_ids = numpy.unique(cluster_id_array)
        num_clusters = len(unique_cluster_ids)

        return num_clusters

    def group_clusters(self, cluster_ids, num_clusters, above_gps, above_val):
        # Group the values by cluster id
        z_clusters = []
        for c_id, c_idxs in generate_group_idxs(cluster_ids):
            c_idxs = flex.size_t(c_idxs)
            c_gps = above_gps.select(c_idxs)
            c_val = above_val.select(c_idxs)
            z_clusters.append((c_gps, c_val))
        assert num_clusters == len(z_clusters)

        return z_clusters


class PanDDADefaultEventFinder:

    def __init__(self, min_blob_z_peak=0, grid_minimum_volume=0, grid_spacing=1, grid=None):

        self.min_blob_z_peak = min_blob_z_peak
        self.grid_minimum_volume = int(float(grid_minimum_volume) / float(grid_spacing)**3)
        self.grid_spacing = grid_spacing

        self.grid = None

    def __call__(self, dataset, num_clusters, z_clusters, grid):

        # TODO: Not very functional...
        self.grid = grid

        print(self.grid_minimum_volume)

        log_strs = []
        print(num_clusters)
        # ============================================================================>
        # Too many points to cluster -- probably a bad dataset
        # ============================================================================>
        if num_clusters == -1:
            # TODO: changed this hard, make sure upstream changes still consistent
            return None

        # ============================================================================>
        #####
        # FILTER/SELECT CLUSTERS OF Z-SCORES
        #####
        # ============================================================================>
        # Filter the clusters by size and peak height
        # ============================================================================>
        if num_clusters > 0:
            num_clusters, z_clusters = self.filter_z_clusters_1(z_clusters=z_clusters)
            self.validate_clusters(z_clusters)
            if num_clusters == 0: log_strs.append('===> Minimum cluster peak/size not reached.')
            print(num_clusters)
        # ============================================================================>
        # Filter the clusters by distance from protein
        # ============================================================================>
        if num_clusters > 0:
            num_clusters, z_clusters = self.filter_z_clusters_2(z_clusters=z_clusters, dataset=dataset)
            self.validate_clusters(z_clusters)
            if num_clusters == 0: log_strs.append('===> Clusters too far from protein.')
            print(num_clusters)
        # ============================================================================>
        # Group Nearby Clusters Together
        # ============================================================================>
        if num_clusters > 0:
            num_clusters, z_clusters = self.group_clusters(z_clusters=z_clusters)
            self.validate_clusters(z_clusters)
            print(num_clusters)
        # ============================================================================>
        # Filter the clusters by symmetry equivalence
        # ============================================================================>
        if num_clusters > 0:
            num_clusters, z_clusters = self.filter_z_clusters_3(z_clusters=z_clusters, dataset=dataset)
            self.validate_clusters(z_clusters)
            print(num_clusters)

        events = self.make_clusters(dataset, z_clusters)
        print(num_clusters)
        return num_clusters, z_clusters, events

    def validate_clusters(self, z_clusters):
        for i, (gps, vals) in enumerate(z_clusters):
            assert len(gps) == len(vals)

    def filter_z_clusters_1(self, z_clusters):
        """Filter the z-clusters on a variety of criteria (size, peak value)"""
        filt_z_clusters = z_clusters

        # Filter out small clusters - get numbers of clusters satisfying the minimum cluster size
        large_clusters = (flex.int([x[1].size() for x in filt_z_clusters]) >= int(self.grid_minimum_volume)).iselection()
        if large_clusters.size() == 0:  return 0, []
        filt_z_clusters = [filt_z_clusters[i] for i in large_clusters]
        # Filter out weak clusters - get numbers of clusters satisfying the minimum z_peak value
        strong_clusters = (flex.double(
            [x[1].min_max_mean().max for x in filt_z_clusters]) >= self.min_blob_z_peak).iselection()
        if strong_clusters.size() == 0: return 0, []
        filt_z_clusters = [filt_z_clusters[i] for i in strong_clusters]

        return len(filt_z_clusters), filt_z_clusters

    def filter_z_clusters_2(self, z_clusters, dataset, min_contact_dist=6):
        """Find and remove clusters more than a minimum distance from the protein"""

        # min_contact_dist - blobs are rejected if they are more than this distance from the protein

        # Extract the protein sites in the reference frame
        ref_sites_cart = dataset.model.alignment.nat2ref(protein(dataset.model.hierarchy).atoms().extract_xyz())
        # Save time - calculate the square of the contact distance
        min_contact_dist_sq = min_contact_dist ** 2

        # Remove any clusters that are more than min_contact_dist from the protein
        filtered_c_idxs = []
        for c_idx, (c_gps, c_val) in enumerate(z_clusters):
            # Extract points in cluster
            cluster_points_cart = self.grid.grid2cart(c_gps)
            # Calculate minimum distance to protein
            for r_site_cart in ref_sites_cart:
                diff_vecs_cart = cluster_points_cart - r_site_cart
                # Keep cluster if minimum distance is less than min_contact_dist
                if min(diff_vecs_cart.dot()) < min_contact_dist_sq:
                    filtered_c_idxs.append(c_idx)
                    break

        filt_z_clusters = [z_clusters[i] for i in filtered_c_idxs]

        return len(filt_z_clusters), filt_z_clusters

    def filter_z_clusters_3(self, z_clusters, dataset, max_contact_dist=8):
        """Find and remove symmetry equivalent clusters"""

        if len(z_clusters) == 1:
            return 1, z_clusters

        # Extract the protein sites in the reference frame
        d_sites_cart = protein(dataset.model.hierarchy).atoms().extract_xyz()
        d_unit_cell = dataset.model.unit_cell
        d_sym_ops = dataset.model.crystal_contact_operators()

        # Cartesianise and fractionalise the points in each of the clusters (in the crystallographic frame)
        points_cart = [None] * len(z_clusters)
        points_frac = [None] * len(z_clusters)
        for c_idx, (c_gps, c_val) in enumerate(z_clusters):
            # Extract points in cluster
            points_cart[c_idx] = dataset.model.alignment.ref2nat(self.grid.grid2cart(c_gps))
            # Fractionalise them to the unit cell of the dataset
            points_frac[c_idx] = d_unit_cell.fractionalize(points_cart[c_idx])
        # Find the sets of clusters that are symmetry related
        sym_equiv_groups = find_symmetry_equivalent_groups(points_frac=points_frac,
                                                           sym_ops=d_sym_ops,
                                                           unit_cell=d_unit_cell,
                                                           cutoff_cart=1.05 * 1.7321 * self.grid_spacing)
        # max_contact_dist - a point contacts an atom if the atoms is within this distance of it
        # Save time - calculate the square of the contact distance
        max_contact_dist_sq = max_contact_dist ** 2
        # Iterate through and chose one from each group to keep
        filt_z_clusters = []
        for g_id, g_idxs in generate_group_idxs(sym_equiv_groups):
            # Count the number of contact for each cluster in the group
            c_contacts = []
            # Iterate through cluster in the group
            for c_idx in g_idxs:
                # Initialise contact counter
                contacts = 0
                # Get the cartesian points for the cluster
                c_points_cart = points_cart[c_idx]
                # Again, use the brute force all-v-all method
                for rp in d_sites_cart:
                    diffs_cart = c_points_cart - rp
                    # Check to see if site closer to cluster than minimum
                    if min(diffs_cart.dot()) < max_contact_dist_sq:
                        contacts += 1
                # Record the number of contacts (over size of cluster)
                c_contacts.append(1.0 * contacts / len(c_points_cart))
            #                if self.log.verbose:
            #                    print('CLUSTER:', c_idx, ', CONTACTS PER POINT:', round(c_contacts[-1],3))

            # Find the cluster with the most contacts
            max_contacts = max(c_contacts)
            if max_contacts == 0:
                raise Exception('MAX CONTACTS IS 0!')
            else:
                cluster_to_keep = g_idxs[c_contacts.index(max_contacts)]
                filt_z_clusters.append(z_clusters[cluster_to_keep])
        #                if self.log.verbose:
        #                    print('KEEPING CLUSTER', cluster_to_keep)
        assert len(filt_z_clusters) == max(
            sym_equiv_groups), 'NUMBER OF UNIQUE GROUPS AND GROUPS TO BE RETURNED NOT THE SAME'

        return len(filt_z_clusters), filt_z_clusters

    def group_clusters(self, z_clusters, separation_cutoff=5):
        """Join clusters that are separated by less than max_separation"""

        if len(z_clusters) == 1:
            return 1, z_clusters

        # Minimum distance between grid points to be joined (squared)
        grid_cutoff_sq = (separation_cutoff / self.grid_spacing) ** 2

        # Record which clusters are to be joined
        connect_array = numpy.zeros((len(z_clusters), len(z_clusters)), dtype=int)
        for i_clust_1, (c_gps_1, c_val_1) in enumerate(z_clusters):
            for i_clust_2, (c_gps_2, c_val_2) in enumerate(z_clusters):
                # Skip if this is the same blob
                if i_clust_1 == i_clust_2:
                    connect_array[(i_clust_1, i_clust_2)] = 1
                    continue
                # Extract the minimum separation of the grid points
                min_dist_sq = min([min((c_gps_2 - gp).dot()) for gp in c_gps_1])
                # Check to see if they should be joined
                if min_dist_sq < grid_cutoff_sq:
                    connect_array[(i_clust_1, i_clust_2)] = 1
        # Cluster the connection array
        cluster_groupings = find_connected_groups(connection_matrix=connect_array)
        # Concatenate smaller clusters into larger clusters
        grouped_clusters = []
        for g_id, g_idxs in generate_group_idxs(cluster_groupings):
            g_gps = []
            [g_gps.extend(z_clusters[i][0]) for i in g_idxs]
            g_gps = flex.vec3_double(g_gps)
            g_val = []
            [g_val.extend(z_clusters[i][1]) for i in g_idxs]
            g_val = flex.double(g_val)
            grouped_clusters.append((g_gps, g_val))

        assert len(grouped_clusters) == max(cluster_groupings)

        return len(grouped_clusters), grouped_clusters

    def make_clusters(self, dataset, z_clusters):

        events = []

        # ============================================================================>
        # Process the identified features
        # ============================================================================>
        for event_idx, (event_points, event_values) in enumerate(z_clusters):
            # Number events from 1
            event_num = event_idx + 1
            # Create a unique identifier for this event
            event_key = (dataset.tag, event_num)
            # ============================================================================>
            # Create a point cluster object
            # ============================================================================>
            point_cluster = PointCluster(id=event_key, points=event_points, values=event_values)
            # ============================================================================>
            # Create an event object
            # ============================================================================>
            event_obj = Event(id=point_cluster.id, cluster=point_cluster)
            event_obj.info.estimated_pseudo_occupancy = 0
            event_obj.info.estimated_bdc = 1.0 - 0
            event_obj.info.global_correlation = 0
            event_obj.info.local_correlation = 0
            # ============================================================================>
            # Append to dataset handler
            # ============================================================================>
            dataset.events.append(event_obj)
            events.append(event_obj)

        return events


class PanDDADefaultBDCCalculator:

    def __init__(self, max_bdc, min_bdc, increment, output_multiplier):
        self.grid = None

        self.max_bdc = max_bdc
        self.min_bdc = min_bdc
        self.increment = increment
        self.output_multiplier = output_multiplier

        self.bdcs = OrderedDict()

    def __call__(self, dataset, dataset_map, ref_map, events, grid):

        self.grid = grid

        # ============================================================================>
        # Extract the map data in non-sparse format
        # ============================================================================>
        dset_map_data = dataset_map.get_map_data(sparse=False)
        ref_map_data = ref_map.get_map_data(sparse=False)
        # ============================================================================>
        # Unpack cluster
        # ============================================================================>
        for event in events[2]:

            # ============================================================================>
            # Estimate the background correction of the detected feature
            # ============================================================================>
            # Extract sites for this cluster and estimate the background correction for the event

            # Generate custom grid mask for this dataset
            event_mask = GridMask(parent=self.grid,
                                  sites_cart=self.grid.grid2cart(event.cluster.points,
                                                                 origin_shift=True),
                                  max_dist=2.0, min_dist=0.0)

            # Select masks to define regions for bdc calculation
            exp_event_idxs = flex.size_t(event_mask.outer_mask_indices())
            reference_idxs = flex.size_t(self.grid.global_mask().inner_mask_indices())
            # ============================================================================>
            # Generate BDC-estimation curve and estimate BDC
            # ============================================================================>
            event_remains, event_corrs, global_corrs = calculate_varying_bdc_correlations(
                ref_map_data=ref_map_data,
                query_map_data=dset_map_data,
                feature_idxs=exp_event_idxs,
                reference_idxs=reference_idxs,
                min_remain=1.0 - self.max_bdc,
                max_remain=1.0 - self.min_bdc,
                bdc_increment=self.increment,
                verbose=True)
            event_remain_est = calculate_maximum_series_discrepancy(
                labels=event_remains,
                series_1=global_corrs,
                series_2=event_corrs)

            event_remain_est = min(event_remain_est * self.output_multiplier,
                                   1.0 - self.min_bdc)
            # ============================================================================>
            # Calculate the map correlations at the selected BDC
            # ============================================================================>
            event_map_data = calculate_bdc_subtracted_map(
                ref_map_data=ref_map_data,
                query_map_data=dset_map_data,
                bdc=1.0 - event_remain_est)
            global_corr = \
            numpy.corrcoef(event_map_data.select(reference_idxs), ref_map_data.select(reference_idxs))[
                0, 1]
            local_corr = numpy.corrcoef(event_map_data.select(exp_event_idxs), ref_map_data.select(exp_event_idxs))[
                0, 1]
            # ============================================================================>
            # Update event parameters
            # ============================================================================>
            event.info.estimated_pseudo_occupancy = event_remain_est
            event.info.estimated_bdc = 1.0-event_remain_est
            event.info.global_correlation = global_corr
            event.info.local_correlation = local_corr

            # ============================================================================>
            # Find the nearest atom to the event
            # ============================================================================>
            # TODO: restore this?
            atm = find_nearest_atoms(atoms=list(protein(dataset.model.hierarchy).atoms_with_labels()),
                                     query=dataset.model.alignment.ref2nat(
                                         self.grid.grid2cart(sites_grid=[map(int, event.cluster.centroid)],
                                                        origin_shift=True)))[0]

            self.bdcs[event.id] = 1.0-event_remain_est

        return self.bdcs