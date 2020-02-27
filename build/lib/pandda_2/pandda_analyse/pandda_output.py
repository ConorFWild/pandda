from collections import OrderedDict

import re
from itertools import chain
import pathlib as p


from collections import OrderedDict

# Standard lib
import os
import glob
import shutil
import pathlib as p
# cctbx
# Custom
from bamboo.common.path import easy_directory, rel_symlink


def symlink(source_path,
            target_path,
            ):
    rel_symlink(str(source_path),
                str(target_path),
                )


def get_pdb_path(dataset):
    return p.Path(dataset.model.file)


def get_mtz_path(dataset):
    return p.Path(dataset.data.file)


def get_dataset_tag(dataset):
    return dataset.tag


def make_target_path(dataset_tag, target_dir_path, extension):
    return target_dir_path / (str(dataset_tag) + "-pandda-input" + extension)


def copy_dataset_files(dataset, target_dir_path):
    source_pdb_path = get_pdb_path(dataset)
    source_mtz_path = get_mtz_path(dataset)
    dataset_tag = get_dataset_tag(dataset)

    target_pdb_path = make_target_path(dataset_tag,
                                       target_dir_path,
                                       ".pdb",
                                       )
    target_mtz_path = make_target_path(dataset_tag,
                                       target_dir_path,
                                       ".mtz",
                                       )

    symlink(source_pdb_path,
            target_pdb_path,
            )

    symlink(source_mtz_path,
            target_mtz_path,
            )


class PanddaOutputSetup:

    def __init__(self, out_dir, lig_style):

        self.name = "PanddaDatasetLoader"

        self.out_dir = out_dir
        self.lig_style = lig_style
        self.datasets = {}

        self.tree = None

        # Variables fpor trace
        self.trace = OrderedDict()

    def __call__(self, mcd):

        # ==============================>
        # Clean previous PanDDA
        # ==============================>
        try:
            shutil.rmtree(self.out_dir,
                          ignore_errors=True)
        except Exception as e:
            print(e)

        # ==============================>
        # Make output dir
        # ==============================>
        os.mkdir(str(self.out_dir))

        # ==============================>
        # Get path objects
        # ==============================>
        dataset_template = {"event_map.ccp4": None,
                            "ligand": {"dummy": None}
                            }

        processed_datasets = {dtag: dataset_template for dtag, d in mcd.datasets.items()}

        analyses = {"pandda_analyse_events.csv": None}

        pandda = {"processed_datasets": processed_datasets,
                  "analyses": analyses}

        # ==============================>
        # Initialise pandda output
        # ==============================>
        self.tree = Tree(str(self.out_dir), pandda)

        # ==============================>
        # Initialise dataset output
        # ==============================>
        errors = []
        for dtag, dataset in mcd.datasets.items():

            # ==============================>
            # Get path
            # ==============================>
            dataset_path = p.Path(self.tree(("processed_datasets", dtag))[0])

            # ==============================>
            # Create links to input files
            # ==============================>
            # Links for the dataset input files
            # TODO: seems inelegant
            link_pdb = str(dataset_path / "pandda_input.pdb")
            link_mtz = str(dataset_path / "pandda_input.mtz")

            # Link the input files to the output folder
            if not os.path.exists(link_pdb): rel_symlink(orig=dataset.model.filename, link=link_pdb)
            if not os.path.exists(link_mtz): rel_symlink(orig=dataset.data.filename, link=link_mtz)
            # ==============================>
            # Search for ligand files and copy them to the output ligands folder
            # ==============================>
            lig_files = glob.glob(os.path.join(os.path.dirname(dataset.model.filename), self.lig_style))
            for lig_file in lig_files:
                # Find all files with the same basename but allowing for different extensions. Then link to output folder.
                lig_base = os.path.splitext(lig_file)[0] + '.*'
                lig_matches = glob.glob(lig_base)
                for lig in lig_matches:
                    out_path = os.path.join(str(dataset_path / 'ligand'), os.path.basename(lig))
                    if os.path.exists(lig) and (not os.path.exists(out_path)):
                        try:
                            shutil.copy(lig, out_path)
                        except:
                            pass
            # # ==============================>
            # # Lastly: Update the pointer to the new path (relative to the pandda directory)
            # # ==============================>
            # dataset.model.filename = os.path.relpath(link_pdb, start=self.out_dir)
            # dataset.data.filename = os.path.relpath(link_mtz, start=self.out_dir)

        return self.tree

    def log(self):

        log = OrderedDict()
        log["pandda_created"] = "Output directory created: \n{}\n".format(p.Path(self.out_dir).exists())
        # log["dataset_dirs"] = "Output directory created: \n{}\n".format("\n".join(["{}".format(p.Path(data_dir).exists())
        #                                                                            for data_dir
        #                                                                            in self.tree(("processed_datasets", ".*"))]))

        return log


class Tree:

    def __init__(self, root, structure):
        self.root = root
        self.structure = structure

        self.build_dirs_recursive(root, structure)

    def __call__(self, res):

        paths = self.get_paths_recursive(self.structure, res)

        root_path = p.Path(self.root)

        # append root
        paths = [str(root_path / path) for path in paths]

        return paths

    def get_paths_recursive(self, struc, regexes):

        if len(regexes) == 1:
            return [key for key in struc if re.match(regexes[0], key)]
        else:
            # get all paths that match the rest of the re
            list_of_lists_of_paths = [[p.Path(key) / x for x in self.get_paths_recursive(struc[key], regexes[1:])]
                                      for key
                                      in struc
                                      if re.match(regexes[0], key)]
            # print(list_of_lists_of_paths)

            # unpack
            list_of_paths = list(chain.from_iterable(list_of_lists_of_paths))

            # return
            return list_of_paths

    def build_dirs_recursive(self, root, structure):

        root_path = p.Path(root)

        for key in structure:
            # Get paths to items
            path = root_path / key

            # recurse if not a file
            if structure[key] is not None:
                try:
                    os.mkdir(str(path))
                    self.build_dirs_recursive(path, structure[key])
                except Exception as e:
                    print(e)

    def update(self, new):

        # Update structure
        self.structure.update(new)

        # Build new
        self.build_dirs_recursive(self.root, new)


class Trace:
    def __init__(self,
                 path,
                 dirct,
                 ):

        self.path = path
        self.dirct = dirct
        # print("got trace with path: {}".format(self.path))

    def __getitem__(self,
                    item,
                    ):
        return Trace(self.path / self.dirct.children[item].name,
                     self.dirct.children[item],
                     )

    def __call__(self, *args):
        return self.dirct(*args,
                          path=self.path
                          )

    def make(self, overwrite=True):
        self.dirct.make(path=self.path,
                        overwrite=overwrite,
                        )


class Dir:

    def __init__(self,
                 name=None,
                 children=None,
                 root=None
                 ):
        self.name = name

        if root:
            self.path = p.Path(root) / name
        else:
            self.path = p.Path(name)

        self.children = children

    def __getitem__(self, item):
        return Trace(self.path / self.children[item].name,
                     self.children[item],
                     )

    def __call__(self, path=None):
        if not path:
            path = self.path
        return path

    def make(self, path=None, overwrite=True):
        if not path:
            path = self.path

        # print(str(path))
        if overwrite:
            shutil.rmtree(str(path),
                          ignore_errors=True,
                          )

        os.mkdir(str(path))

        for child_name, child in self.children.items():
            Trace(path / child.name,
                  child,
                  ).make(overwrite=overwrite)


class File:

    def __init__(self,
                 name=None,
                 root=None
                 ):
        self.name = name

        if root:
            self.path = p.Path(root) / self.name
        else:
            self.path = p.Path(self.name)

    def __call__(self, path=None):
        if not path:
            path = self.path
        return path

    def make(self, path=None, overwrite=True):
        return


class IndexedFile:
    def __init__(self,
                 name=None,
                 root=None
                 ):
        self.name = name

        if root:
            self.path = p.Path(root) / self.name
        else:
            self.path = p.Path(self.name)

    def __call__(self, index, path=None):
        if not path:
            path = self.path
        return p.Path(str(path).format(index))

    def make(self, path=None, overwrite=True):
        return

#
# def make(file_obj):
#
#     try:
#     except Exception as e:
#
#
#     try:
#
#
#
#         os.mkdir(str(file_obj.path))

#
#
# def output_tree(out_dir):
#
#     root = StaticDir(out_dir, [])
#
#     datasets = KeyDir()
#     shells = KeyDir(name,
#                   root=hierachy,
#                     children={StaticDir(dtag, children={"zmap": File("zmap"),
#                                                         }
#                                         )
#                               }
#                     )
#
#     hierarchy = StaticDir(root=out_dir,
#                           children={"datasets": datasets,
#                                     "shells": shells,
#                                     }
#                           )
#
#     make(hierarchy)
#
#     return hierarchy