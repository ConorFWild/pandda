import os, subprocess, pathlib

import numpy as np
from biopandas.pdb import PandasPdb


class QFitParameters:
    def __init__(self,
                 density_cutoff,
                 ):
        self.density_cutoff = 0


class AutobuildQFit:
    def __init__(self,
                 ):
        pass

    def __call__(self,
                 protein_model_path,
                 ligand_model_path,
                 event_map_path,
                 output_dir_path,
                 resolution,
                 event,
                 qfit_parameters=QFitParameters(density_cutoff=0),
                 ):
        placed_ligand_model_path = get_placed_ligand_model_path(ligand_model_path)

        place_ligand(ligand_model_path,
                     event,
                     placed_ligand_model_path,
                     )

        command_string = "/dls/science/groups/i04-1/conor_dev/ccp4/base/bin/qfit_ligand {event_map} {resolution} {ligand_model_path} -r {protein_model_path} -d {output_dir_path}"
        command_string.format(event_map=event_map_path,
                              resolution=resolution,
                              ligand_model_path=ligand_model_path,
                              protein_model_path=protein_model_path,
                              output_dir_path=output_dir_path,
                              )

        p = subprocess.Popen(command_string,
                             shell=True,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             )

        stdout, stderr = p.communicate()


def place_ligand(ligand_model_path,
                 event,
                 output_path,
                 ):
    ligand = PandasPdb().read_pdb(str(ligand_model_path))

    ligand_coords = ligand.df["HETATM"][["x_coord", "y_coord", "z_coord"]]

    ligand_com = ligand_coords.mean(axis=1)


    print(event)
    event_com = np.array(event.x, event.y, event.z)

    translation = event_com - ligand_com

    ligand_coords["x_coord", "y_coord", "z_coord"] += translation

    ligand.to_pdb(path=output_path,
                  records=None,
                  gz=False,
                  append_newline=True,
                  )


def get_placed_ligand_model_path(ligand_model_path):
    placed_ligand_model_path = pathlib.Path(ligand_model_path).parent / "placed_ligand.pdb"
    return placed_ligand_model_path
