from collections import OrderedDict

# Standard lib

import pathlib as p
# Custom
from path import easy_directory, rel_symlink


class CopyPDB:
    def __init__(self):
        pass

    def __call__(self, dataset, target_dir_path):
        source_pdb_path = self.get_pdb_path(dataset)

        # print(source_pdb_path)

        dataset_tag = get_dataset_tag(dataset)

        target_pdb_path = make_target_path(dataset_tag,
                                           target_dir_path,
                                           ".pdb",
                                           )

        try:
            symlink(source_pdb_path,
                    target_pdb_path,
                    )

        except:
            pass

    def get_pdb_path(self, dataset):

        return p.Path(dataset.model.filename)


class CopyMTZ:
    def __init__(self):
        pass

    def __call__(self, dataset, target_dir_path):
        source_mtz_path = self.get_mtz_path(dataset)

        # print(source_mtz_path)

        dataset_tag = get_dataset_tag(dataset)

        target_mtz_path = make_target_path(dataset_tag,
                                           target_dir_path,
                                           ".mtz",
                                           )

        try:
            symlink(source_mtz_path,
                    target_mtz_path,
                    )
        except:
            pass

    def get_mtz_path(self, dataset):

        return p.Path(dataset.data.filename)


class CopyLigands:
    def __init__(self):
        pass

    def __call__(self, dataset, target_dir_path):
        try:

            source_ligand_path = self.get_ligand_path(dataset)

            # print(source_ligand_path)

            rel_symlink(str(source_ligand_path),
                        str(target_dir_path),
                        )
        except:
            pass

        try:
            initial_ligand_path = self.get_local_ligand_path(dataset)

            # print(initial_ligand_path)

            rel_symlink(str(initial_ligand_path),
                        str(target_dir_path),
                        )
        except:
            pass

    def get_ligand_path(self,
                        dataset,
                        ligand_dir_name="compound",
                        ):
        return p.Path(dataset.data.filename).parent / ligand_dir_name

    def get_local_ligand_path(self,
                              dataset,
                              ligand_dir_name="compound",
                              ):
        return p.Path(dataset.data.filename).parent.glob("*.cif").next()


def symlink(source_path,
            target_path,
            ):
    rel_symlink(str(source_path),
                str(target_path),
                )


def get_dataset_tag(dataset):
    return dataset.tag


def make_target_path(dataset_tag, target_dir_path, extension):
    return target_dir_path / (str(dataset_tag) + "-pandda-input" + extension)


class DatasetFileCopier:
    def __init__(self,
                 copy_pdb=CopyPDB(),
                 copy_mtz=CopyMTZ(),
                 copy_ligand=CopyLigands(),
                 ):
        self.copy_pdb = copy_pdb
        self.copy_mtz = copy_mtz
        self.copy_ligand = copy_ligand

    def __call__(self,
                 dataset,
                 tree,
                 ):
        for dtag, d in dataset.datasets.items():
            target_dir = tree["datasets"][dtag]()

            self.copy_pdb(d,
                          target_dir,
                          )
            self.copy_mtz(d,
                          target_dir,
                          )
            self.copy_ligand(d,
                             target_dir,
                             )

    def repr(self):
        repr = OrderedDict()
        return repr
