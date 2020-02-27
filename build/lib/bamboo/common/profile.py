import cProfile, pstats, StringIO

class profile_code(object):
    def __init__(self):
        self.start()

    def start():
        self.profiler = cProfile.Profile()
        self.profiler.enable()

    def stop(self, print_stats=True):
        self.profiler.disable()
        s = StringIO.StringIO()
        sortby = 'cumulative'
        ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
        if print_stats:
            ps.print_stats()
            print s.getvalue()
        return ps
