from collections import OrderedDict

import pandas as pd


class MakeEventTable:
    def __init__(self):
        pass

    def __call__(self,
                 event_table_path,
                 dataset,
                 events,
                 grid,
                 events_analysed,
                 analysed_resolution,
                 ):
        # Produce the event table
        # dir_path = p.Path(tree([str(name)])[0])
        # event_table_path = dir_path / "event_table.csv"

        event_table = PanDDADefaultEventTableShell()(dataset.partition_datasets("test"),
                                                     events,
                                                     event_table_path,
                                                     grid,
                                                     events_analysed,
                                                     analysed_resolution,
                                                     )

        return event_table

    def repr(self):
        repr = OrderedDict()
        return repr


class PanDDADefaultEventTableShell:

    def __init__(self, order_by=None):
        self.order_by = None
        self.event_info = None

    def __call__(self,
                 datasets,
                 events,
                 table_path,
                 grid,
                 events_analysed,
                 analysed_resolution,
                 ):

        event_table = self.make_event_table(datasets,
                                            events,
                                            grid,
                                            events_analysed,
                                            analysed_resolution,
                                            )
        self.write_event_table(table_path,
                               event_table,
                               )

        return event_table

    def make_event_table(self,
                         datasets,
                         events,
                         grid,
                         events_analysed,
                         analysed_resolution,
                         ):
        # for dtag, dataset in datasets.items():
        #     print(dtag)
        #     # ============================================================================>
        #     # Mark as interesting and add events to the event table
        #     # ============================================================================>
        #     # TODO: dataset.all_masks() will need to be changed
        #     records = []
        #     if dtag in events:
        #         print(events[dtag])
        #         for e in events[dtag][2]:
        #             print(e)
        #             event_record = self.get_event_record(dataset=dataset,
        #                                                  event=e)
        #             records.append(event_record)

        records = []
        for dtag, events in events.items():
            for event in events[2]:
                event_record = self.get_event_record(dataset=datasets[dtag],
                                                     event=event,
                                                     grid=grid,
                                                     analysis=events_analysed[dtag][event.id],
                                                     analysed_resolution=analysed_resolution,
                                                     )
                records.append(event_record)

        return pd.DataFrame(records)

    def get_event_record(self,
                         dataset,
                         event,
                         grid,
                         analysis,
                         analysed_resolution,
                         ):
        """Add event entries to the event table"""

        # Default to site_idx of 0 if no site given
        if event.parent:
            site_idx = event.parent.id
        else:
            site_idx = 0

        # Get native frame loc
        xyz = list(dataset.model.alignment.ref2nat(
            coordinates=grid.grid2cart([event.cluster.peak],
                                       origin_shift=True,
                                       ))[0]
                   )

        # Generate record
        record = {"dtag": dataset.tag,
                  "event_idx": event.id[1],
                  "site_idx": site_idx,
                  "1-BDC": round(analysis["estimated_pseudo_occupancy"], 2),
                  "z_peak": round(event.cluster.max, 2),
                  "z_mean": round(event.cluster.mean, 2),
                  "cluster_size": event.cluster.size,
                  "x": xyz[0],
                  "y": xyz[1],
                  "z": xyz[2],
                  "global_correlation_to_average_map": analysis["global_corr"],
                  "local_correlation_to_average_map": analysis["local_corr"],
                  "analysed_resolution": analysed_resolution,
                  "map_uncertainty": 0,
                  "refx": 0,
                  "refy": 0,
                  "refz": 0,
                  "r_work": dataset.model.input.get_r_rfree_sigma().r_work,
                  "r_free": dataset.model.input.get_r_rfree_sigma().r_free,
                  "rejected - total": False,
                  "noisy zmap": False,
                  "analysed": False,
                  "interesting": False,
                  "exclude_from_zmap_analysis": False,
                  "exclude_from_characterisation": False,
                  }

        return record

    def add_event_to_event_table(self, dataset, event, grid):
        """Add event entries to the event table"""

        # Check event has not been added previously
        self.event_info.loc[event.id, :] = None
        # Default to site_idx of 0 if no site given
        if event.parent:
            site_idx = event.parent.id
        else:
            site_idx = 0
        self.event_info.set_value(event.id, 'site_idx', site_idx)
        # Event and cluster information
        self.event_info.set_value(event.id, '1-BDC', round(1.0 - event.info.estimated_bdc, 2))
        self.event_info.set_value(event.id, 'z_peak', round(event.cluster.max, 2))
        self.event_info.set_value(event.id, 'z_mean', round(event.cluster.mean, 2))
        self.event_info.set_value(event.id, 'cluster_size', event.cluster.size)
        #        self.tables.event_info.set_value(event.id, ['refx','refy','refz'], list(self.grid.grid2cart([event.cluster.peak],origin_shift=False)[0]))
        self.event_info.set_value(event.id, ['x', 'y', 'z'], list(dataset.model.alignment.ref2nat(
            coordinates=grid.grid2cart([event.cluster.peak], origin_shift=True))[0]))
        self.event_info.set_value(event.id, 'global_correlation_to_average_map',
                                  event.info.global_correlation)
        self.event_info.set_value(event.id, 'local_correlation_to_average_map', event.info.local_correlation)

    def write_event_table(self, table_path, event_table):
        # Write the event data only once events have been recorded
        if len(event_table):
            # Sort the event data by z-peak and write out
            sort_eve = event_table.sort_values(by=['dtag', "event_idx"],
                                               ascending=[1, 1],
                                               )
            print(sort_eve)

            # TODO: not sure if removing comb_tab dependancy is reasonable
            # sort_eve = sort_eve.join(comb_tab, how='inner')

            sort_eve.to_csv(path_or_buf=str(table_path))

        else:
            standin_event_table = pd.DataFrame(columns=["dtag",
                                                        "event_idx",
                                                        "site_idx", "1-BDC",
                                                        "z_peak",
                                                        "z_mean",
                                                        "cluster_size",
                                                        "x",
                                                        "y",
                                                        "z",
                                                        "global_correlation_to_average_map",
                                                        "local_correlation_to_average_map",
                                                        "analysed_resolution",
                                                        "map_uncertainty",
                                                        "refx",
                                                        "refy",
                                                        "refz",
                                                        "r_work",
                                                        "r_free",
                                                        "rejected - total",
                                                        "noisy zmap",
                                                        "analysed",
                                                        "interesting",
                                                        "exclude_from_zmap_analysis",
                                                        "exclude_from_characterisation",
                                                        ],
                                               )

            standin_event_table.to_csv(path_or_buf=str(table_path))

