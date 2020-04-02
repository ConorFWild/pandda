class OutputEventTable:
    def __init__(self):
        pass

    def __call__(self,
                 event_table,
                 event_table_path,
                 ):

        event_table_sorted = event_table.sort_values("cluster_size", ascending=False)

        event_table_sorted.to_csv(str(event_table_path))