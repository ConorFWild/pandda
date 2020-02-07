import os

from libtbx import easy_pickle

from bamboo.common.logs import Log
from bamboo.common.file import FileManager

class Program(object):
    """Class meant to provide basic functionality for programs and pipelines"""

    _NAME = None
    _TEXT = None
    _VERSION = None

    _allowed_statuses = ['running','done','errored']

    log = Log()

    file_manager = None

    def write_running_parameters_to_log(self, params):
        self.log.heading('Processed parameters')
        self.log(self.master_phil.format(python_object=params).as_str())
        self.log.heading('Parameters different to the defaults')
        self.log(self.master_phil.fetch_diff(source=self.master_phil.format(python_object=params)).as_str())

    def check_for_matplotlib(self, backend=None, interactive=False):
        """Check to see whether we can load matplotlib"""
        self.log('Checking for matplotlib:')
        try:
            import matplotlib
            matplotlib.interactive(interactive)
            from matplotlib import pyplot
            if backend:
                pyplot.switch_backend(backend)
                current_backend = pyplot.get_backend()
                assert current_backend == backend, 'Backend loaded ({}) is not the one requested ({})'.format(current_backend, backend)
            assert pyplot.isinteractive() is interactive, 'Interactive setting is incorrect ({} is not {})'.format(pyplot.isinteractive(), interactive)
            pyplot.style.use('ggplot')
            self.log('pyplot loaded successfully. Using backend "{!s}"'.format(current_backend))
            return True
        except:
            self.log('===================================>>>')
            self.log('>> COULD NOT IMPORT MATPLOTLIB. WILL NOT BE ABLE TO GENERATE GRAPHS.')
            self.log('===================================>>>')
            return False

    def initialise_file_manager(self, rootdir):
        self.file_manager = FileManager(rootdir=rootdir)
        return self.file_manager

    def update_status(self, status):
        """Set log files to indicate the status of the program"""

        assert status in self._allowed_statuses
        # Delete any that may exist
        existing_files = [self.file_manager.get_file('status').format(f) for f in self._allowed_statuses]
        [os.remove(f) for f in existing_files if os.path.exists(f)]
        # Create the new  status file
        with open(self.file_manager.get_file('status').format(status), 'w') as fh: fh.write('')

    def pickle(self, pickle_file, pickle_object, overwrite=True):
        """Takes an object and pickles it"""
        if os.path.exists(pickle_file) and not overwrite:
            self.log('NOT PICKLING: {!s}'.format(os.path.relpath(pickle_file, start=self.out_dir)))
        else:
            self.log('Pickling Object: {!s}'.format(os.path.relpath(pickle_file, start=self.out_dir)))
            easy_pickle.dump(pickle_file, pickle_object)

    def unpickle(self, pickle_file):
        """Takes an object and unpickles it"""
        self.log('Unpickling File: {!s}'.format(os.path.relpath(pickle_file, start=self.out_dir)))
        return easy_pickle.load(pickle_file)

