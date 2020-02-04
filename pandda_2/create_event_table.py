from collections import OrderedDict

import pandas as pd


class CreateEventTable:
    def __init__(self):
        pass

    def __call__(self,
                 tree,
                 num_shells,
                 ):
        event_tables_computed = {i: pd.read_csv(tree["shells"][i]["event_table"]())
                                 for i
                                 in range(num_shells)
                                 }

        joined_event_table = merge_event_tables(event_tables_computed)

        return joined_event_table

    def __repr__(self):
        repr = OrderedDict()
        return repr


def merge_event_tables(event_tables):

    print([len(event_table) for shell_num, event_table in event_tables.items()])

    event_table = pd.concat([shell_table
                             for shell_num, shell_table
                             in event_tables.items()
                             if (len(shell_table) > 0)
                             ]
                            )

    event_table = event_table.fillna(0)

    # event_table = event_table.set_index(["dtag", "event_idx"])
    return event_table
