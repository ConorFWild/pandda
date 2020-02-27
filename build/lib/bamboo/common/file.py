import os, datetime


class FileManager(object):


    def __init__(self, rootdir):
        """Creates and stores output files to be neatly created and retrieved"""

        self.topdir = rootdir
        self.output_dirs = {'root':rootdir}
        self.dir_parents = {'root':None}
        self.output_files = {}

    def add_dir(self, dir_name, dir_tag, top_dir_tag=None, create=True, exists=True):
        """Store a directory `dir_name` under `dir_tag` in directory under `top_dir_tag`"""
        # Check that it hasn't been created already
        assert dir_tag not in self.output_dirs.keys(), 'Directory already added: {}'.format(dir_tag)
        # If no directory given, put in roots
        if top_dir_tag is None: top_dir_tag = 'root'
        # Record the directory's parent
        self.dir_parents[dir_tag] = top_dir_tag
        # Get the path of the parent directory
        top_dir_path = self.output_dirs[top_dir_tag]
        # Create dirname and store
        self.output_dirs[dir_tag] = os.path.join(top_dir_path, dir_name)
        # Create if it doesn't exist
        if create and (not os.path.exists(self.output_dirs[dir_tag])):
            self._make_directory_if_necessary(dir_tag)
        if exists:
            assert os.path.exists(self.output_dirs[dir_tag])

    def get_dir(self, dir_tag):
        """Retrieve a dirname by it's dir_tag"""
        assert dir_tag in self.output_dirs.keys(), 'Directory has not been added: {}'.format(dir_tag)
        return self.output_dirs[dir_tag]

    def add_file(self, file_name, file_tag, dir_tag=None):
        """Store a filename `file_name` under `file_tag` in directory under `dir_tag`"""
        # Check that it hasn't been created already
        assert file_tag not in self.output_files.keys(), 'File already added: {}'.format(file_tag)
        # Check the directory that it's beeing added to exists
        if dir_tag is None: dir_tag = 'root'
        dir_name = self.output_dirs[dir_tag]
        # Create filename and store
        self.output_files[file_tag] = FileObj(file=os.path.join(dir_name, file_name), tag=file_tag)

    def get_file(self, file_tag):
        """Retrieve a filename by it's file_tag"""
        return self.get_file_object(file_tag=file_tag).path

    def get_file_object(self, file_tag):
        """Retrieve a file object by it's file_tag"""
        assert file_tag in self.output_files.keys(), 'File has not been added: {}'.format(file_tag)
        return self.output_files[file_tag]

    def check_and_create_directories(self):
        """Check that all directories exist and create them where necessary"""
        for dir_tag in self.output_dirs.keys():
            self._make_directory_if_necessary(dir_tag=dir_tag)

    def _make_directory_if_necessary(self, dir_tag):
        if (dir_tag != 'root') and (not os.path.exists(self.output_dirs[self.dir_parents[dir_tag]])):
            self._make_directory_if_necessary(self.dir_parents[dir_tag])
        if not os.path.exists(self.output_dirs[dir_tag]):
            os.mkdir(self.output_dirs[dir_tag])


class FileObj(object):


    def __init__(self, file, tag=None):
        """File information object"""

        self.input = file

        self.path = os.path.abspath(file)
        self.tag = tag
        self.dir, self.name = os.path.split(self.path)
        self.base, self.ext = os.path.splitext(self.name)

    def __str__(self):
        return self.path

    def __call__(self):
        return self.path

    def get_datestamp(self, date='created'):
        stats = os.lstat(self.path)
        if date=='created':  return datetime.datetime.fromtimestamp(stats.st_ctime)
        if date=='modified': return datetime.datetime.fromtimestamp(stats.st_mtime)
        if date=='accessed': return datetime.datetime.fromtimestamp(stats.st_atime)
        raise Exception("Invalid option: date='{}'".format(date))

