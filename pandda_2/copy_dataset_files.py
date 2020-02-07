from collections import OrderedDict

# Standard lib
import subprocess
import pathlib as p
# Custom
from path import easy_directory, rel_symlink

import iotbx.cif
from iotbx.file_reader import any_file

from pandda_2 import processor


class CopyPDB:
    def __init__(self):
        pass

    def __call__(self, dataset, target_dir_path, output_path):
        source_pdb_path = self.get_pdb_path(dataset)

        # print(source_pdb_path)

        dataset_tag = get_dataset_tag(dataset)

        target_pdb_path = make_target_path(dataset_tag,
                                           target_dir_path,
                                           ".pdb",
                                           )

        try:
            symlink(source_pdb_path,
                    output_path,
                    )


        except Exception as e:

            print(e)

    def get_pdb_path(self, dataset):

        return p.Path(dataset.model.filename)


class CopyMTZ:
    def __init__(self):
        pass

    def __call__(self, dataset, target_dir_path, output_path):
        source_mtz_path = self.get_mtz_path(dataset)

        # print(source_mtz_path)

        dataset_tag = get_dataset_tag(dataset)

        target_mtz_path = make_target_path(dataset_tag,
                                           target_dir_path,
                                           ".mtz",
                                           )

        try:
            symlink(source_mtz_path,
                    output_path,
                    )
        except Exception as e:
            print(e)

    def get_mtz_path(self, dataset):

        return p.Path(dataset.data.filename)


class CopyLigands:
    def __init__(self,
                 regex="*.cif",
                 ):
        self.regex = regex

    def __call__(self, dataset, output_path, regex):
        try:

            source_ligand_path = self.get_ligand_path(dataset,
                                                      regex,
                                                      )

            # print(source_ligand_path)

            rel_symlink(str(source_ligand_path),
                        str(output_path / source_ligand_path.name),
                        )

        except Exception as e:
            print(e)

        # try:
        #     initial_ligand_path = self.get_local_ligand_path(dataset)
        #
        #     # print(initial_ligand_path)
        #
        #     rel_symlink(str(initial_ligand_path),
        #                 str(output_path),
        #                 )
        # except:
        #     pass

    def get_ligand_path(self,
                        dataset,
                        regex="**/ligand.cif",
                        ):
        # print(p.Path(dataset.data.filename).parent)
        # print(self.regex)
        return next(p.Path(dataset.data.filename).parent.glob(regex))

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


class CifToPDBPhenix:

    def __init__(self,
                 compound_dir,
                 ):
        self.compound_dir = compound_dir

    def __call__(self, *args, **kwargs):
        compound_dir = self.compound_dir

        processes = []
        for cif_path in compound_dir.glob("*.cif"):
            print("\t{}".format(cif_path))

            # pdb_out_path = cif_path.parent / "{}.pdb".format(cif_path.name)

            command_string = "cd {compound_dir}; module load phenix; phenix.elbow {cif_path} --output=\"{cif_stem}\""

            formatted_command = command_string.format(compound_dir=compound_dir,
                                                      cif_path=cif_path,
                                                      cif_stem=cif_path.stem,
                                                      )

            print(formatted_command)

            p = subprocess.Popen(formatted_command,
                                 shell=True,
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 )
            processes.append(p)

        output = [process.communicate()
                  for process
                  in processes
                  ]
        print(output)

        # cif_file = open(str(cif_path), "r")
        # cif_string = cif_file.read()
        # print("\t{}".format(cif_string))
        #
        # structures = iotbx.cif.reader(input_string=cif_string).build_crystal_structures()
        # print(structures)
        #
        # structure_file = any_file(str(cif_path))
        # print(structure_file)
        # print(dir(structure_file))
        #
        # structure = structure_file.file_object
        #
        # print(structure)
        # print(dir(structure))
        #
        # crystal_strucutres = structure.build_crystal_structures()
        #
        # print(crystal_strucutres)
        #
        # builder = structure.builder
        # print(builder)
        #
        # model = structure.model
        # print(model)
        #
        # # print(builder())
        # print(model())
        # print(dir(model()))
        #
        # if len(structures) == 0:
        #     print("no structures: Skipping!")
        #     continue
        #
        # else:
        #     structure = structures[0]
        # print("Converting to pdb")
        # pdb = structure.as_pdb_file(str(pdb_out_path))
        # print(pdb)

        # parser = MMCIFParser()
        # structure = parser.get_structure("test",
        #                                  str(cif_path),
        #                                  )
        #
        # io = PDBIO()
        # io.set_structure(structure)
        # io.save(pdb_out_path)


class DatasetFileCopier:
    def __init__(self,
                 copy_pdb=CopyPDB(),
                 copy_mtz=CopyMTZ(),
                 copy_ligand=CopyLigands(),
                 cif_to_pdb=CifToPDBPhenix,
                 processor=processor.ProcessorJoblib(),
                 autobuild=False
                 ):
        self.copy_pdb = copy_pdb
        self.copy_mtz = copy_mtz
        self.copy_ligand = copy_ligand
        self.cif_to_pdb = cif_to_pdb
        self.processor = processor
        self.autobuild = autobuild
        # TODO: make these settable from main
        self.ligand_cif_regex = "**/ligand.cif"
        self.ligand_pdb_regex = "**/ligand.pdb"

    def __call__(self,
                 dataset,
                 tree,
                 ):
        for dtag, d in dataset.datasets.items():
            target_dir = tree["processed_datasets"][dtag]()
            print("\tCopying pdb for dataset {}".format(dtag))
            self.copy_pdb(d,
                          target_dir,
                          tree["processed_datasets"][dtag]["initial_model"](),
                          )
            print("\tCopying mtz for dataset {}".format(dtag))

            self.copy_mtz(d,
                          target_dir,
                          tree["processed_datasets"][dtag]["initial_data"](),
                          )
            print("\tCopying ligand for dataset {}".format(dtag))

            self.copy_ligand(d,
                             tree["processed_datasets"][dtag]["ligand_files"](),
                             regex=self.ligand_cif_regex,
                             )

            self.copy_ligand(d,
                             tree["processed_datasets"][dtag]["ligand_files"](),
                             regex=self.ligand_pdb_regex,
                             )

            self.copy_ligand(d,
                             tree["processed_datasets"][dtag]["autobuilt"](),
                             self.ligand_cif_regex,
                             )

            print("\tConverting cifs to ligands for dataset {}".format(dtag))

        # Convert any cifs to usable pdbs
        if self.autobuild:
            self.processor([self.cif_to_pdb(compound_dir=tree["processed_datasets"][dtag]["ligand_files"]())
                            for dtag
                            in dataset.datasets.keys()
                            ]
                           )

    def repr(self):
        repr = OrderedDict()
        return repr
