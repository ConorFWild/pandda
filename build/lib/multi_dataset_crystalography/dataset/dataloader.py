import re

from collections import OrderedDict

# Standard lib
import os
import glob

# Custom

from dataset import PanddaDataset

from multi_dataset_crystalography.constants import *


def load_dataset(pdb=None, mtz=None, dtag=None, num=None):
    # ==============================>
    # Load datasets in parallel
    # ==============================>
    dataset = PanddaDataset.from_file(model_filename=pdb,
                                      data_filename=mtz).label(num=num, tag=dtag)

    # ==============================>
    # Initialise loaded datasets
    # ==============================>

    # ==============================>
    # Catch errors and print at end
    # ==============================>

    # ==============================>
    # Check that dataset has a CRYST line and unit cell information
    # ==============================>
    if dataset.model.crystal_symmetry is None:
        # errors.append((dataset,
        #                'Could not load crystal symmetry information - check pdb CRYST line is present and valid.'))
        raise Exception("Broken dataset")
    if dataset.model.unit_cell is None:
        # errors.append(
        #     (dataset, 'Could not load unit cell information - check pdb CRYST line is present and valid.'))
        raise Exception("Broken dataset")

    # ==============================>
    # Intialise the meta for the dataset
    # ==============================>
    dataset.meta.analysed = False
    dataset.meta.dataset_info = None
    dataset.meta.dataset_map_info = None

    return dataset


def parse_pandda_dir(data_dirs,
                     pdb_style=None,
                     mtz_style=None,
                     pdb_regex=None,
                     mtz_regex=None,
                     dir_regex=None,
                     only_datasets=None,
                     ignore_datasets=None,
                     dataset_prefix=None,
                     ):
    # ==============================>
    # Extract input styles from parameter object
    # ==============================>
    dir_style = data_dirs.strip('/')
    pdb_style = pdb_style.strip('/')
    mtz_style = mtz_style.strip('/')
    # ==============================>
    # Find datasets in the input directories
    # ==============================>
    new_files = []
    empty_directories = []
    for dir in sorted(glob.glob(data_dirs)):
        pdb_files = [f for f in glob.glob(os.path.join(dir, pdb_style)) if os.path.exists(f)]
        mtz_files = [f for f in glob.glob(os.path.join(dir, mtz_style)) if os.path.exists(f)]
        if not (pdb_files and mtz_files):
            # print('EMPTY DIRECTORY: {!s}'.format(dir))
            empty_directories.append(dir)
        elif not pdb_files:
            # print('NO PDB IN DIRECTORY: {!s}'.format(dir))
            empty_directories.append(dir)
        elif not mtz_files:
            # print('NO MTZ IN DIRECTORY: {!s}'.format(dir))
            empty_directories.append(dir)
        else:
            assert len(pdb_files) == 1, 'More than one matching PDB file found: {!s}'.format(
                os.path.join(dir, pdb_style))
            assert len(mtz_files) == 1, 'More than one matching MTZ file found: {!s}'.format(
                os.path.join(dir, mtz_style))
            # ==============================>
            # Found PDB anf MTZ file in directory
            # ==============================>
            new_pdb = pdb_files[0]
            new_mtz = mtz_files[0]
            dataset_tag = [None]
            # ==============================>
            # Regex Matching - PDB file
            # ==============================>
            if '*' in pdb_style:
                pdb_base = os.path.basename(new_pdb)
                if pdb_regex:
                    pdb_regex = pdb_regex
                else:
                    pdb_regex = pdb_style.replace('*', '(.*)')
                pdb_tag = re.findall(pdb_regex, pdb_base)
                assert pdb_tag, 'NO PDB TAG FOUND: {!s} -> {!s}'.format(pdb_regex, pdb_base)
                if isinstance(pdb_tag[0], tuple):
                    pdb_tag = list(pdb_tag[0])[0:1]

            else:
                pdb_regex = pdb_tag = None
            # ==============================>
            # Regex Matching - MTZ file
            # ==============================>
            if '*' in mtz_style:
                mtz_base = os.path.basename(new_mtz)
                if mtz_regex:
                    mtz_regex = mtz_regex
                else:
                    mtz_regex = mtz_style.replace('*', '(.*)')
                mtz_tag = re.findall(mtz_regex, mtz_base)
                assert mtz_tag, 'NO MTZ TAG FOUND: {!s} -> {!s}'.format(mtz_regex, mtz_base)
                if isinstance(mtz_tag[0], tuple):
                    mtz_tag = list(mtz_tag[0])[0:1]
            else:
                mtz_regex = mtz_tag = None
            # ==============================>
            # Regex Matching - Directory
            # ==============================>
            if '*' in dir_style:
                dir_base = os.path.dirname(pdb_files[0])
                if dir_regex:
                    dir_regex = dir_regex
                else:
                    dir_regex = dir_style.replace('*', '(.*)')
                dir_tag = re.findall(dir_regex, dir_base)
                assert dir_tag, 'NO DIR TAG FOUND: {!s} -> {!s}'.format(dir_regex, dir_base)
                if isinstance(dir_tag[0], tuple):
                    dir_tag = list(dir_tag[0])[0:1]
            else:
                dir_regex = dir_tag = None
            # ==============================>
            # Check consistency
            # ==============================>
            if pdb_tag and mtz_tag: assert pdb_tag == mtz_tag, 'PDB-MTZ TAGS ARE NOT IDENTICAL: {} != {}'.format(
                pdb_tag, mtz_tag)
            if dir_tag and pdb_tag: assert dir_tag == pdb_tag, 'DIR-PDB TAGS ARE NOT IDENTICAL: {} != {}'.format(
                dir_tag, pdb_tag)
            if dir_tag and mtz_tag: assert dir_tag == mtz_tag, 'DIR-MTZ TAGS ARE NOT IDENTICAL: {} != {}'.format(
                dir_tag, mtz_tag)
            # ==============================>
            # Extract tag
            # ==============================>
            if dir_tag:
                dataset_tag = dir_tag
            elif pdb_tag:
                dataset_tag = pdb_tag
            elif mtz_tag:
                dataset_tag = mtz_tag
            # ==============================>
            # Add prefix - slightly obsoleted
            # ==============================>
            dataset_tag[0] = str(dataset_tag[0])
            if isinstance(dataset_tag[0], str):
                if dataset_prefix is None:
                    dataset_tag = [dataset_tag[0]]
                else:
                    dataset_tag = [dataset_prefix + dataset_tag[0]]
            else:
                assert dataset_tag[0] is None

            pdb_files = [str(p_file) for p_file in pdb_files]
            mtz_files = [str(p_file) for p_file in mtz_files]

            new_files.append(pdb_files + mtz_files + dataset_tag)

    filtered_new_files = new_files
    # new_files = new_files

    # ==============================>
    # Select only those requested
    # ==============================>
    if only_datasets:
        only_tags = only_datasets.split(',')

        re_filtered_new_files = []
        for i, (pdb, mtz, tag) in enumerate(filtered_new_files):
            if tag in only_tags:
                re_filtered_new_files.append(filtered_new_files[i])

        filtered_new_files = re_filtered_new_files

    # ==============================>
    # Filter out manually labelled datasets to ignore
    # ==============================>
    if ignore_datasets:
        ignore_tags = ignore_datasets.split(',')

        re_filtered_new_files = []
        for i, (pdb, mtz, tag) in enumerate(filtered_new_files):
            if not (tag in ignore_tags):
                re_filtered_new_files.append(filtered_new_files[i])

        filtered_new_files = re_filtered_new_files

    return filtered_new_files


class DefaultDataloader:

    def __init__(self,
                 data_dirs, pdb_style, mtz_style, pdb_regex, mtz_regex, dir_regex, only_datasets, ignore_datasets,
                 dataset_prefix, out_dir, lig_style):
        self.trace = OrderedDict()

        self.name = "DefaultDataloader"

        # self.data_dirs = config["args"]["data_dirs"]
        # self.pdb_style = config["args"]["pdb_style"]
        # self.mtz_style = config["args"]["mtz_style"]
        # self.pdb_regex = config["args"]["pdb_regex"]
        # self.mtz_regex = config["args"]["mtz_regex"]
        # self.dir_regex = config["args"]["dir_regex"]
        # self.only_datasets = config["args"]["only_datasets"]
        # self.ignore_datasets = config["args"]["ignore_datasets"]
        # self.dataset_prefix = config["args"]["dataset_prefix"]
        self.data_dirs = data_dirs
        self.pdb_style = pdb_style
        self.mtz_style = mtz_style
        self.pdb_regex = pdb_regex
        self.mtz_regex = mtz_regex
        self.dir_regex = dir_regex
        self.only_datasets = only_datasets
        self.ignore_datasets = ignore_datasets
        self.dataset_prefix = dataset_prefix

        # self.out_dir = config["args"]["out_dir"]
        # self.lig_style = config["args"]["lig_style"]
        self.out_dir = out_dir
        self.lig_style = lig_style

    def __call__(self):
        pandda_input_directory_parser = PanddaParseInputDirectory(self.data_dirs,
                                                                  self.pdb_style,
                                                                  self.mtz_style,
                                                                  self.pdb_regex,
                                                                  self.mtz_regex,
                                                                  self.dir_regex,
                                                                  self.only_datasets,
                                                                  self.ignore_datasets,
                                                                  self.dataset_prefix)
        filtered_new_files = pandda_input_directory_parser()
        # self.trace[pandda_input_directory_parser.name] = pandda_input_directory_parser.log()

        pandda_dataset_loader = PanddaDatasetLoader(self.out_dir, self.lig_style)
        pandda_datasets = pandda_dataset_loader(filtered_new_files)
        # self.trace[pandda_dataset_loader.name] = pandda_dataset_loader.log()

        return pandda_datasets

    def __repr__(self):
        repr = {"data_dirs": self.data_dirs,
                "pdb_style": self.pdb_style,
                "mtz_style": self.mtz_style,
                }
        return repr

    # def log(self):
    #
    #     return self.trace


class PanddaParseInputDirectory:

    def __init__(self, data_dirs, pdb_style, mtz_style, pdb_regex, mtz_regex, dir_regex, only_datasets, ignore_datasets,
                 dataset_prefix):

        self.name = "PanddaParseInputDirectory"

        self.data_dirs = data_dirs
        self.pdb_style = pdb_style
        self.mtz_style = mtz_style
        self.pdb_regex = pdb_regex
        self.mtz_regex = mtz_regex
        self.dir_regex = dir_regex
        self.only_datasets = only_datasets
        self.ignore_datasets = ignore_datasets
        self.dataset_prefix = dataset_prefix

    def __call__(self):
        """Builds a list of input files from the command line arguments passed"""

        # ==============================>
        # Extract input styles from parameter object
        # ==============================>
        dir_style = self.data_dirs.strip('/')
        pdb_style = self.pdb_style.strip('/')
        mtz_style = self.mtz_style.strip('/')
        # ==============================>
        # Find datasets in the input directories
        # ==============================>
        new_files = []
        empty_directories = []
        for dir in sorted(glob.glob(self.data_dirs)):
            pdb_files = [f for f in glob.glob(os.path.join(dir, pdb_style)) if os.path.exists(f)]
            mtz_files = [f for f in glob.glob(os.path.join(dir, mtz_style)) if os.path.exists(f)]
            if not (pdb_files and mtz_files):
                # print('EMPTY DIRECTORY: {!s}'.format(dir))
                empty_directories.append(dir)
            elif not pdb_files:
                # print('NO PDB IN DIRECTORY: {!s}'.format(dir))
                empty_directories.append(dir)
            elif not mtz_files:
                # print('NO MTZ IN DIRECTORY: {!s}'.format(dir))
                empty_directories.append(dir)
            else:
                assert len(pdb_files) == 1, 'More than one matching PDB file found: {!s}'.format(
                    os.path.join(dir, pdb_style))
                assert len(mtz_files) == 1, 'More than one matching MTZ file found: {!s}'.format(
                    os.path.join(dir, mtz_style))
                # ==============================>
                # Found PDB anf MTZ file in directory
                # ==============================>
                new_pdb = pdb_files[0]
                new_mtz = mtz_files[0]
                dataset_tag = [None]
                # ==============================>
                # Regex Matching - PDB file
                # ==============================>
                if '*' in pdb_style:
                    pdb_base = os.path.basename(new_pdb)
                    if self.pdb_regex:
                        pdb_regex = self.pdb_regex
                    else:
                        pdb_regex = pdb_style.replace('*', '(.*)')
                    pdb_tag = re.findall(pdb_regex, pdb_base)
                    assert pdb_tag, 'NO PDB TAG FOUND: {!s} -> {!s}'.format(pdb_regex, pdb_base)
                    if isinstance(pdb_tag[0], tuple):
                        pdb_tag = list(pdb_tag[0])[0:1]

                else:
                    pdb_regex = pdb_tag = None
                # ==============================>
                # Regex Matching - MTZ file
                # ==============================>
                if '*' in mtz_style:
                    mtz_base = os.path.basename(new_mtz)
                    if self.mtz_regex:
                        mtz_regex = self.mtz_regex
                    else:
                        mtz_regex = mtz_style.replace('*', '(.*)')
                    mtz_tag = re.findall(mtz_regex, mtz_base)
                    assert mtz_tag, 'NO MTZ TAG FOUND: {!s} -> {!s}'.format(mtz_regex, mtz_base)
                    if isinstance(mtz_tag[0], tuple):
                        mtz_tag = list(mtz_tag[0])[0:1]
                else:
                    mtz_regex = mtz_tag = None
                # ==============================>
                # Regex Matching - Directory
                # ==============================>
                if '*' in dir_style:
                    dir_base = os.path.dirname(pdb_files[0])
                    if self.dir_regex:
                        dir_regex = self.dir_regex
                    else:
                        dir_regex = dir_style.replace('*', '(.*)')
                    dir_tag = re.findall(dir_regex, dir_base)
                    assert dir_tag, 'NO DIR TAG FOUND: {!s} -> {!s}'.format(dir_regex, dir_base)
                    if isinstance(dir_tag[0], tuple):
                        dir_tag = list(dir_tag[0])[0:1]
                else:
                    dir_regex = dir_tag = None
                # ==============================>
                # Check consistency
                # ==============================>
                if pdb_tag and mtz_tag: assert pdb_tag == mtz_tag, 'PDB-MTZ TAGS ARE NOT IDENTICAL: {} != {}'.format(
                    pdb_tag, mtz_tag)
                if dir_tag and pdb_tag: assert dir_tag == pdb_tag, 'DIR-PDB TAGS ARE NOT IDENTICAL: {} != {}'.format(
                    dir_tag, pdb_tag)
                if dir_tag and mtz_tag: assert dir_tag == mtz_tag, 'DIR-MTZ TAGS ARE NOT IDENTICAL: {} != {}'.format(
                    dir_tag, mtz_tag)
                # ==============================>
                # Extract tag
                # ==============================>
                if dir_tag:
                    dataset_tag = dir_tag
                elif pdb_tag:
                    dataset_tag = pdb_tag
                elif mtz_tag:
                    dataset_tag = mtz_tag
                # ==============================>
                # Add prefix - slightly obsoleted
                # ==============================>
                dataset_tag[0] = str(dataset_tag[0])
                if isinstance(dataset_tag[0], str):
                    if self.dataset_prefix is None:
                        dataset_tag = [dataset_tag[0]]
                    else:
                        dataset_tag = [self.dataset_prefix + dataset_tag[0]]
                else:
                    assert dataset_tag[0] is None

                pdb_files = [str(p_file) for p_file in pdb_files]
                mtz_files = [str(p_file) for p_file in mtz_files]

                new_files.append(pdb_files + mtz_files + dataset_tag)

        filtered_new_files = new_files
        self.new_files = new_files

        # ==============================>
        # Select only those requested
        # ==============================>
        if self.only_datasets:
            only_tags = self.only_datasets.split(',')

            re_filtered_new_files = []
            for i, (pdb, mtz, tag) in enumerate(filtered_new_files):
                if tag in only_tags:
                    re_filtered_new_files.append(filtered_new_files[i])

            filtered_new_files = re_filtered_new_files

        # ==============================>
        # Filter out manually labelled datasets to ignore
        # ==============================>
        if self.ignore_datasets:
            ignore_tags = self.ignore_datasets.split(',')

            re_filtered_new_files = []
            for i, (pdb, mtz, tag) in enumerate(filtered_new_files):
                if not (tag in ignore_tags):
                    re_filtered_new_files.append(filtered_new_files[i])

            filtered_new_files = re_filtered_new_files

        return filtered_new_files

    # def log(self):
    #     log = OrderedDict()
    #     log["data_dirs"] = "Seaching for datasets in: \n{}\n".format(str(self.data_dirs))
    #     log["found_datasets"] = "Found datasets: \n{}\n".format("\n".join([x[2]
    #                                                                        for x
    #                                                                        in self.new_files]))
    #     log["only_datasets"] = "Restricting to only analyse: \n{}\n".format("\n".join(self.only_datasets))
    #     log["ignored_datasets"] = "Manually ignoring datasets: \n{}\n".format("\n".join(self.ignore_datasets))
    #
    #     return log


class PanddaDatasetLoader:

    def __init__(self, out_dir, lig_style):

        self.name = "PanddaDatasetLoader"

        self.out_dir = out_dir
        self.lig_style = lig_style
        self.datasets = {}

        # Variables fpor trace
        self.errors = []
        self.loaded_dataset_tags = []

    def __call__(self, new_files):
        """Read in maps for the input datasets"""

        # ==============================>
        # Load datasets in parallel
        # ==============================>
        loaded_datasets = {dtag: PanddaDataset.from_file(model_filename=pdb,
                                                         data_filename=mtz).label(num=num, tag=dtag)
                           for num, (pdb, mtz, dtag)
                           in enumerate(new_files)}

        # ==============================>
        # Initialise loaded datasets
        # ==============================>

        # ==============================>
        # Catch errors and print at end
        # ==============================>
        errors = []
        for dtag, dataset in loaded_datasets.items():
            # ==============================>
            # Check that dataset has a CRYST line and unit cell information
            # ==============================>
            if dataset.model.crystal_symmetry is None:
                errors.append((dataset,
                               'Could not load crystal symmetry information - check pdb CRYST line is present and valid.'))
                continue
            if dataset.model.unit_cell is None:
                errors.append(
                    (dataset, 'Could not load unit cell information - check pdb CRYST line is present and valid.'))
                continue
            # ==============================>
            # Intialise the meta for the dataset
            # ==============================>
            dataset.meta.analysed = False
            dataset.meta.dataset_info = None
            dataset.meta.dataset_map_info = None

            # ==============================>
            # Record succesful laoding
            # ==============================>
            self.loaded_dataset_tags.append(dataset.tag)

        return loaded_datasets

    # def log(self):
    #
    #     log = OrderedDict()
    #     log["Loaded Datasets"] = "Loaded datasets : \n{}\n".format("\n".join(self.loaded_dataset_tags))
    #     log["Errors"] = "Could not load datases: \n{}\n".format("\n".join(self.errors))
    #
    #     return log
