class OutputEventTable:
    def __init__(self):
        pass

    def __call__(self,
                 event_table,
                 event_table_path,
                 ):
        event_table.to_csv(str(event_table_path))