class OutputSitesTable:
    def __init__(self):
        pass

    def __call__(self,
                 sites_table,
                 sites_table_path,
                 ):
        sites_table.to_csv(sites_table_path)