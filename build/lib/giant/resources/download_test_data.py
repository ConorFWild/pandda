
import os, sys, time

from bamboo.common.path import easy_directory
from bamboo.common.logs import Log

# Adapted from https://blog.shichao.io/2012/10/04/progress_speed_indicator_for_urlretrieve_in_python.html
class _ReportHook(object):
    def __init__(self):
        self.start_time = time.time()
    def __call__(self, count, block_size, total_size):
        duration = time.time() - self.start_time
        progress_size = int(count * block_size)
        speed = int(progress_size / (1024 * duration))
        percent = int(count * block_size * 100 / total_size)
        sys.stdout.write("\r...%d%%, %d MB, %d KB/s, %d seconds passed" %
                        (percent, progress_size / (1024 * 1024), speed, duration))
        sys.stdout.flush()
        if count * block_size >= total_size:
            sys.stdout.write("\n")

class Dataset(object):
    def __init__(self, id, type):
        self.id=id
        self.type=type

class ZenodoDataset(Dataset):
    def __init__(self, id, type, zenodo_id, log=None):
        if log is None: log=Log()
        self.log = log
        super(ZenodoDataset, self).__init__(id=id, type=type)
        self.zenodo_id = zenodo_id
        self.base_url = "https://zenodo.org/record/{}/".format(zenodo_id)
        self.data_url = self.base_url + "files/data.zip"
        self.data_dir = None

    def download(self, output_directory):
        out_dir = easy_directory(output_directory)
        self.data_dir = os.path.join(out_dir, 'data')
        if os.path.exists(self.data_dir):
            self.log('Output data directory already exists -- skipping')
            return
        source_file = self.data_url
        target_file = os.path.join(out_dir, 'data.zip')
        self.log('Downloading zipfile: {} -> {}'.format(source_file, target_file))
        import urllib
        urllib.urlretrieve(source_file, target_file, reporthook=_ReportHook())
        self.log('Unzipping: {} -> {}/'.format(target_file, out_dir))
        import zipfile
        zf = zipfile.ZipFile(target_file)
        zf.extractall(path=out_dir)
        assert os.path.join(self.data_dir), 'Data directory should have been created by this process!'
        os.remove(target_file)

        return self.data_dir

    def show_summary(self, log=None):
        if log is None: log = self.log
        log('Dataset ID: {}'.format(self.id))
        log('Dataset type: {}'.format(self.type))
        log('Dataset Url: {}'.format(self.base_url))
        log('Download location: {}'.format(self.data_dir))

class DatasetList(object):
    def __init__(self, datasets=None):
        self.datasets = []
        self._hash = {}
        if datasets is not None:
            self.add(datasets)

    def add(self, new_datasets):
        if isinstance(new_datasets, DatasetList):
            new_datasets = new_datasets.datasets
        elif not isinstance(new_datasets, list):
            new_datasets = [new_datasets]
        for d in new_datasets:
            assert not self._hash.has_key(d.id)
            self._hash[d.id] = len(self.datasets)
            self.datasets.append(d)
            assert self.get(d.id) is d

    def get(self, id):
        assert self._hash.has_key(id)
        return self.datasets[self._hash.get(id)]

    def show_summary(self, log=None):
        if log is None: log = Log()
        log.subheading('Available datasets')
        for d in self.datasets:
            log.bar()
            d.show_summary(log=log)
        log.bar()

zenodo_dataset_information = [('BAZ2B', 'fragment screen', 48768),
                              ('BRD1A', 'fragment screen', 48769),
                              ('JMJD2', 'fragment screen', 48770),
                              ('SP100', 'fragment screen', 48771)]
# Create class for fragment screens
fragment_screens = DatasetList([ZenodoDataset(*v) for v in zenodo_dataset_information])
# Create class for all test datasets
all_datasets = DatasetList()
all_datasets.add(fragment_screens)


