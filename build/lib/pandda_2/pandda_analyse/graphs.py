class PanDDARunGraphs:

    def __init__(pandda_statistics):  # Produces a set of graphs of statistics

        self.trace = {}

        for statistic in pandda_statistics:
            try:
                statistic.output()
            except Exception as e:
                self.trace[statistic.name] = "{}".format(e)

    def log:

        ...
        "{statistic}: {graphed}".format()