from collections import OrderedDict

import pandas as pd

from pandda_2 import pandda_logging


class ProcessShell:
    def __init__(self,
                 diffraction_data_truncator,
                 reference_map_getter,
                 get_reference_map,
                 map_loader,
                 statistical_model,
                 fit_model,
                 evaluate_model,
                 cluster_outliers,
                 filter_clusters,
                 event_analyser,
                 make_event_map,
                 make_z_map,
                 make_mean_map,
                 make_event_table,
                 process,
                 ):
        # Data
        self.shell_dataset = None
        self.reference = None
        self.grid = None
        self.tree = None
        self.shell_num = None

        # Functions
        self.diffraction_data_truncator = diffraction_data_truncator
        self.reference_map_getter = reference_map_getter
        self.get_reference_map = get_reference_map

        self.map_loader = map_loader
        self.statistical_model = statistical_model
        self.fit = fit_model
        self.evaluate_model = evaluate_model
        self.cluster_outliers = cluster_outliers
        self.filter_clusters = filter_clusters
        self.event_analyser = event_analyser

        self.make_event_map = make_event_map
        self.make_z_map = make_z_map
        self.make_mean_map = make_mean_map
        self.make_event_table = make_event_table

        self.process = process

    def __call__(self,
                 shell_dataset,
                 reference,
                 grid,
                 tree,
                 shell_num,
                 ):

        self.shell_dataset = shell_dataset
        self.reference = reference
        self.grid = grid
        self.tree = tree
        self.shell_num = shell_num

        try:
            print("Checking for event table at: {}".format(tree["shells"][self.shell_num]["event_table"]()))
            event_table = pd.read_csv(tree["shells"][self.shell_num]["event_table"]())
            return event_table
        except:
            print("Shell {} not yet processed!".format(self.shell_num))

        # ###############################################
        # Get resolution
        # ###############################################
        resolutions_test = max([dts.data.summary.high_res
                                for dtag, dts
                                in self.shell_dataset.partition_datasets("test").items()
                                ]
                               )
        resolutions_train = max([dts.data.summary.high_res
                                 for dtag, dts
                                 in self.shell_dataset.partition_datasets("train").items()
                                 ]
                                )
        max_res = max(resolutions_test,
                      resolutions_train,
                      )

        # ###############################################
        # Instantiate sheel variable names
        # ###############################################
        dtags = set(self.shell_dataset.partition_datasets("test").keys()
                    + self.shell_dataset.partition_datasets("train").keys()
                    )

        train_dtags = [dtag
                       for dtag
                       in dtags
                       if (dtag in self.shell_dataset.partition_datasets("train").keys())
                       ]
        test_dtags = [dtag
                      for dtag
                      in dtags
                      if (dtag in self.shell_dataset.partition_datasets("test").keys())
                      ]
        # print("\tNumber of train dtags: {}".format(len(train_dtags)))
        # print("\tNumber of test dtags: {}".format(len(test_dtags)))

        # ###############################################
        # Truncate datasets
        # ###############################################
        print("Truncaitng Diffractions")
        truncated_reference, truncated_datasets = self.diffraction_data_truncator(self.shell_dataset.datasets,
                                                                                  self.reference,
                                                                                  )

        shell_max_res = max_res

        # ###############################################
        # Generate maps
        # ###############################################
        print("Getting Shell refeence map")
        shell_ref_map = self.get_reference_map(self.reference_map_getter,
                                               self.reference,
                                               shell_max_res,
                                               self.grid,
                                               )

        print("Getting XMaps")
        xmap_funcs = {}
        for dtag in truncated_datasets.keys():
            task = TaskWrapper(self.map_loader,
                               truncated_datasets[dtag],
                               self.grid,
                               shell_ref_map,
                               shell_max_res,
                               )
            xmap_funcs[dtag] = task

        # print("\tNum xmap funcs: {}".format(len(xmap_funcs)))
        # print("\tNumber of truncated datasets is: {}".format(len(truncated_datasets)))
        # print("\tNum of keys in xmaps_funcs is:{}".format(xmap_funcs.keys()))
        xmaps = self.process(xmap_funcs)

        # ###############################################
        # Fit statistical model to trianing sets
        # ###############################################
        print("Fitting model")
        # print("\tNumber of xmaps: {}".format(len(xmaps)))

        shell_fit_model = self.fit(self.statistical_model,
                                   [xmaps[dtag] for dtag in train_dtags],
                                   [xmaps[dtag] for dtag in test_dtags],
                                   )

        shell_fit_model_scattered = shell_fit_model

        # ###############################################
        # Find events
        # ###############################################
        print("Evaluating model")
        zmap_funcs = {}
        for dtag in test_dtags:
            task = TaskWrapper(self.evaluate_model,
                               shell_fit_model_scattered,
                               xmaps[dtag],
                               )
            zmap_funcs[dtag] = task
        zmaps = self.process(zmap_funcs)

        print("Getting clusters")
        clusters_funcs = {}
        for dtag in test_dtags:
            task = TaskWrapper(self.cluster_outliers,
                               truncated_datasets[dtag],
                               zmaps[dtag],
                               self.grid,
                               )
            clusters_funcs[dtag] = task
        clusters = self.process(clusters_funcs)

        print("Getting events")
        filter_funcs = {}
        for dtag in test_dtags:
            task = TaskWrapper(self.filter_clusters,
                               truncated_datasets[dtag],
                               clusters[dtag],
                               self.grid,
                               )
            filter_funcs[dtag] = task
        events = self.process(filter_funcs)

        print("Analysign events")
        analyse_funcs = {}
        for dtag in test_dtags:
            task = TaskWrapper(self.event_analyser,
                               truncated_datasets[dtag],
                               xmaps[dtag],
                               shell_ref_map,
                               events[dtag],
                               self.grid,
                               )
            analyse_funcs[dtag] = task
        events_analysed = self.process(analyse_funcs)

        events_analysed_computed = events_analysed

        events_computed = events

        # Criticise each indidual dataset (generate statistics, event map and event table)
        z_maps_ccp4 = OrderedDict()

        event_map_args = OrderedDict()
        z_map_args = OrderedDict()
        event_maps = OrderedDict()

        z_map_funcs = OrderedDict()
        event_map_funcs = OrderedDict()

        for dtag in test_dtags:
            event_maps[dtag] = OrderedDict()

            if len(events_computed[dtag][2]) != 0:
                # print("Making z map args")
                z_map_funcs[dtag] = TaskWrapper(self.make_z_map,
                                                xmaps[dtag],
                                                truncated_datasets[dtag],
                                                self.tree["processed_datasets"][dtag]["z_map"](),
                                                shell_fit_model_scattered,
                                                self.grid,
                                                )

            for event_id, events_analysis_computed in events_analysed_computed[dtag].items():
                # print("Print making event maps args")
                pandda_logging.log_map_making(xmaps[dtag],
                                              shell_ref_map,
                                              events_analysed_computed[dtag][event_id]["estimated_bdc"],
                                              )

                map_bdc = events_analysed_computed[dtag][event_id]["estimated_bdc"]
                fractional_mean_map = shell_ref_map.new_from_template(map_data=shell_ref_map.data * map_bdc,
                                                                      sparse=shell_ref_map.is_sparse(),
                                                                      )

                event_map = xmaps[dtag] - fractional_mean_map

                print(pandda_logging.log_summarise_map(event_map))

                event_map_funcs[(dtag, event_id)] = TaskWrapper(self.make_event_map,
                                                                xmaps[dtag],
                                                                truncated_datasets[dtag],
                                                                shell_ref_map,
                                                                events_computed[dtag][2][
                                                                    int(event_id[1]) - 1],
                                                                events_analysed_computed[dtag][
                                                                    event_id]["estimated_bdc"],
                                                                self.tree["processed_datasets"][dtag]["event_map"]([dtag,
                                                                                 event_id[1],
                                                                                 round(1-events_analysed_computed[dtag][event_id]["estimated_bdc"], 2),
                                                                                                                    ]
                                                                                 ),
                                                                self.grid,
                                                                )

        print("Print making event maps")
        event_maps = self.process(event_map_funcs)

        print("Maziing z map")
        z_maps_ccp4 = self.process(z_map_funcs)

        print("Making mean maps")
        shell_maps = self.make_mean_map(self.reference,
                                        shell_ref_map,
                                        self.grid,
                                        self.tree["shells"][self.shell_num]["mean_map"](),
                                        )

        print("Making event map")
        event_table = self.make_event_table(self.tree["shells"][self.shell_num]["event_table"](),
                                            self.shell_dataset,
                                            events_computed,
                                            self.grid,
                                            events_analysed_computed,
                                            analysed_resolution=max_res,
                                            )

        return event_table

    def repr(self):
        repr = OrderedDict()
        repr["diffraction_data_truncator"] = self.diffraction_data_truncator.repr()
        repr["reference_map_getter"] = self.reference_map_getter.repr()
        repr["get_reference_map"] = self.get_reference_map.repr()

        repr["map_loader"] = self.map_loader.repr()
        repr["statistical_model"] = self.statistical_model.repr()
        repr["fit"] = self.fit.repr()
        repr["evaluate_model"] = self.evaluate_model.repr()
        repr["cluster_outliers"] = self.cluster_outliers.repr()
        repr["filter_clusters"] = self.filter_clusters.repr()
        repr["event_analyser"] = self.event_analyser.repr()

        repr["make_event_map"] = self.make_event_map.repr()
        repr["make_z_map"] = self.make_z_map.repr()
        repr["make_mean_map"] = self.make_mean_map.repr()
        repr["make_event_table"] = self.make_event_table.repr()

        repr["process"] = self.process.repr()

        return repr


class TaskWrapper:
    def __init__(self, func, *args, **kwargs):
        self.func = func
        self.args = args
        self.kwargs = kwargs

    def __call__(self):
        return self.func(*self.args, **self.kwargs)
