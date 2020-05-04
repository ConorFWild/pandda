import os
import argparse
import re
import subprocess
import shutil
from pathlib import Path

import pandas as pd

import joblib

# from biopandas.pdb import PandasPdb
import gemmi


class AutobuildingCommand:
    def __init__(self,
                 out_dir_path=None,
                 mtz_path=None,
                 ligand_path=None,
                 receptor_path=None,
                 coord=(0, 0, 0),
                 ):
        env = "module load phenix"
        ligand_fit_command = "phenix.ligandfit"
        ligand_fit_args = "data={mtz} ligand={ligand} model={receptor} search_center=[{x},{y},{z}] search_dist=6"
        ligand_fit_args_formatted = ligand_fit_args.format(mtz=mtz_path,
                                                           ligand=ligand_path,
                                                           receptor=receptor_path,
                                                           x=coord[0],
                                                           y=coord[1],
                                                           z=coord[2],
                                                           )
        self.command = "{env}; cd {out_dir_path}; {ligand_fit_command} {args}".format(env=env,
                                                                                      out_dir_path=out_dir_path,
                                                                                      ligand_fit_command=ligand_fit_command,
                                                                                      args=ligand_fit_args_formatted,
                                                                                      )

    def __repr__(self):
        return self.command


def execute(command):
    # print("\t\t: {}".format(command))
    submit_proc = subprocess.Popen(str(command),
                                   shell=True,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE,
                                   )
    stdout, stderr = submit_proc.communicate()
    return stdout, stderr


def event_map_to_mtz(event_map_path: Path,
                     output_path,
                     resolution,
                     col_f="FWT",
                     col_ph="PHWT",
                     gemmi_path: Path = "/dls/science/groups/i04-1/conor_dev/gemmi/gemmi",
                     ):
    command = "module load gcc/4.9.3; source /dls/science/groups/i04-1/conor_dev/anaconda/bin/activate env_clipper_no_mkl; {gemmi_path} map2sf {event_map_path} {output_path} {col_f} {col_ph} --dmin={resolution}"
    formatted_command = command.format(gemmi_path=gemmi_path,
                                       event_map_path=event_map_path,
                                       output_path=output_path,
                                       col_f=col_f,
                                       col_ph=col_ph,
                                       resolution=resolution,
                                       )

    p = subprocess.Popen(formatted_command,
                         shell=True,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE,
                         )

    stdout, stderr = p.communicate()

    return output_path


def autobuild_event(event):
    event_mtz_path = event.pandda_event_dir / "{}_{}.mtz".format(event.dtag, event.event_idx)

    event_map_to_mtz(event.event_map_path,
                     event_mtz_path,
                     event.analysed_resolution,
                     )

    out_dir_path = event.pandda_event_dir / "autobuild_event_{}".format(event.event_idx)

    try:
        shutil.rmtree(str(out_dir_path))
    except:
        pass

    os.mkdir(str(out_dir_path))

    autobuilding_command = AutobuildingCommand(out_dir_path=out_dir_path,
                                               mtz_path=event_mtz_path,
                                               ligand_path=event.ligand_path,
                                               receptor_path=event.receptor_path,
                                               coord=event.coords,
                                               )

    stdout, stderr = execute(autobuilding_command)

    result = AutobuildingResult(event,
                                stdout,
                                stderr,
                                )

    return result


class Option:
    def __init__(self, short, long, required=False, default=None, help=""):
        self.short = short
        self.long = long
        self.required = required
        self.default = default
        self.help = help

    def __call__(self, args):
        return vars(args)[self.long]


class Config:
    def __init__(self):
        input_pandda_dir = Option("i", "input_pandda_dir", required=True)
        options = [input_pandda_dir]

        args = self.get_args(options)

        self.input_pandda_dir = Path(input_pandda_dir(args))

    def get_args(self, options):
        parser = argparse.ArgumentParser()
        for option in options:
            parser.add_argument("-{}".format(option.short),
                                "--{}".format(option.long),
                                help=option.help,
                                required=option.required,
                                default=option.default,
                                )
        return parser.parse_args()


class PanDDAFilesystemModel:
    def __init__(self, pandda_root_dir):
        self.pandda_root_dir = pandda_root_dir

        self.pandda_analyse_dir = pandda_root_dir / "analyses"
        self.pandda_inspect_events_path = self.pandda_analyse_dir / "pandda_inspect_events.csv"
        self.autobuilding_results_table = self.pandda_analyse_dir / "autobuilding_results.csv"

        self.pandda_processed_datasets_dir = pandda_root_dir / "processed_datasets"
        self.pandda_processed_datasets_dirs = list(self.pandda_analyse_dir.glob("*"))


def map_parallel(f, datasets):
    results = joblib.Parallel(n_jobs=20,
                              verbose=50,
                              )(joblib.delayed(f)(dataset)
                                for dataset
                                in datasets)

    return results


class Event:
    def __init__(self,
                 pandda_event_dir,
                 dtag,
                 event_idx,
                 event_map_path,
                 ligand_path,
                 receptor_path,
                 coords,
                 analysed_resolution,
                 ):
        self.dtag = dtag
        self.event_idx = event_idx
        self.pandda_event_dir = pandda_event_dir
        self.event_map_path = event_map_path
        self.ligand_path = ligand_path
        self.receptor_path = receptor_path
        self.coords = coords
        self.analysed_resolution = analysed_resolution


def get_event_table(path):
    return pd.read_csv(str(path))


def get_event_map_path(pandda_event_dir, dtag, event_idx, occupancy):
    event_map_path = pandda_event_dir / "{}-event_{}_1-BDC_{}_map.native.ccp4".format(dtag,
                                                                                      event_idx,
                                                                                      occupancy,
                                                                                      )
    return event_map_path


def get_ligand_path(pandda_event_dir):
    event_dir = pandda_event_dir

    ligands = list((event_dir / "ligand_files").glob("*.pdb"))
    ligand_strings = [str(ligand_path) for ligand_path in ligands if ligand_path.name != "tmp.pdb"]

    if len(ligand_strings) == 0:
        return None

    ligand_pdb_path: Path = Path(min(ligand_strings,
                                     key=len,
                                     )
                                 )
    return ligand_pdb_path


def get_receptor_path(pandda_event_dir, dtag):
    receptor_path = pandda_event_dir / "{}-pandda-input.pdb".format(dtag)
    return receptor_path


def get_coords(row):
    return (row["x"], row["y"], row["z"])


def get_analyed_resolution(row):
    return row["analysed_resolution"]


def get_events(event_table, fs):
    events = []
    for index, row in event_table.iterrows():
        dtag = row["dtag"]
        event_idx = row["event_idx"]
        occupancy = row["1-BDC"]
        pandda_event_dir = fs.pandda_processed_datasets_dir / "{}".format(dtag)
        event_map_path = get_event_map_path(pandda_event_dir,
                                            dtag, event_idx, occupancy
                                            )
        ligand_path = get_ligand_path(pandda_event_dir)
        if ligand_path is None:
            continue
        receptor_path = get_receptor_path(pandda_event_dir, dtag)
        coords = get_coords(row)
        analysed_resolution = get_analyed_resolution(row)

        event = Event(pandda_event_dir,
                      dtag,
                      event_idx,
                      event_map_path,
                      ligand_path,
                      receptor_path,
                      coords,
                      analysed_resolution,
                      )
        events.append(event)

    return events


class AutobuildingResult:
    def __init__(self, event: Event, stdout, stderr):
        event_autobuilding_dir = event.pandda_event_dir / "autobuild_event_{}".format(event.event_idx)
        event_ligandfit_dir = event_autobuilding_dir / "LigandFit_run_1_"
        autobuilding_results_file = event_ligandfit_dir / "LigandFit_summary.dat"

        with open(str(autobuilding_results_file), "r") as f:
            result_string = f.read()

        rscc_regex = "[\s]+1[\s]+[0-9\.]+[\s]+([0-9\.]+)"
        match = re.findall(rscc_regex, result_string)
        rscc_string = match[0]
        self.dtag = event.dtag
        self.event_idx = event.event_idx
        self.rscc = float(rscc_string)
        self.stdout = stdout
        self.stderr = stderr


class ResultsTable:
    def __init__(self, results):
        records = []
        for result in results:
            record = {}
            record["dtag"] = result.dtag
            record["event_idx"] = result.event_idx
            record["rscc"] = result.rscc

            records.append(record)

        self.table = pd.DataFrame(records)

    def to_csv(self, path):
        self.table.to_csv(str(path))


def get_highest_rscc_events(events,
                            results_table,
                            ):
    unique_dtags = results_table["dtag"].unique()

    events_map = {(event.dtag, event.event_idx): event for event in events}

    max_events = []
    for unique_dtag in unique_dtags:
        events_table = results_table[results_table["dtag"] == unique_dtag]
        max_rscc_label = events_table["rscc"].idxmax()
        event_row = events_table.loc[max_rscc_label]
        dtag = event_row["dtag"]
        event_idx = event_row["event_idx"]
        max_events.append(events_map[(dtag, event_idx)])

    return max_events


# def copy_event_to_processed_models(event: Event, fs):
#     event_autobuilding_dir = event.pandda_event_dir / "autobuild_event_{}".format(event.event_idx)
#     event_ligandfit_dir = event_autobuilding_dir / "LigandFit_run_1_"
#     event_build_path = event_ligandfit_dir / "ligand_fit_1.pdb"
#
#     initial_model_path = event_autobuilding_dir / "{}-pandda-input.pdb".format(event.dtag)
#
#     pandda_inspect_model_dir = event.pandda_event_dir / "modelled_structures"
#     pandda_inspect_model_path = pandda_inspect_model_dir / "{}-pandda-model.pdb".format(event.dtag)
#
#     initial_model = PandasPdb().read_pdb(str(initial_model_path))
#     best_autobuild_model = PandasPdb().read_psb(str(event_build_path))
#
#     initial_model.df["HETATM"] = initial_model.df["HETATM"].append(best_autobuild_model.df["HETATM"])
#
#     renumber(initial_model.df["HETATM"])
#
#     initial_model.to_pdb(str(pandda_inspect_model_path))

def merge_model(event, fs):
    # Define initial model path, event model path and output path
    initial_model_path = event.pandda_event_dir / "{}-pandda-input.pdb".format(event.dtag)
    autobuilds_dir = event.pandda_event_dir / "autobuild_event_{}".format(event.event_idx)
    autobuild_dir = autobuilds_dir / "LigandFit_run_1_"
    autobuild_path = autobuild_dir / "ligand_fit_1.pdb"
    output_path = event.pandda_event_dir / "modelled_structures" / "{}-pandda-model.pdb".format(event.dtag)

    # Read inital model
    initial_structure = gemmi.gemmi.read_structure(str(initial_model_path))
    initial_model = initial_structure[0]

    # Read event model
    event_structure = gemmi.gemmi.read_structure(str(autobuild_path))
    event_model = event_structure[0]
    ligand_chain = event_model[0]

    # Create new chain and append event residue
    initial_model.add_chain(ligand_chain,
                            pos=-1,
                            )

    # Save to pandda modelled structures dir
    initial_structure.write_pdb(str(output_path))


def merge_models(events,
                 results_table,
                 fs,
                 ):
    highest_rscc_events = get_highest_rscc_events(events,
                                                  results_table,
                                                  )

    print("\t\tAfter filetering duplicate events got {} events".format(len(highest_rscc_events)))
    for event in highest_rscc_events:
        merge_model(event, fs)


def main():
    config = Config()

    print("Building I/O model...")
    fs = PanDDAFilesystemModel(config.input_pandda_dir)
    print("\tFound {} dataset dirs".format(len(fs.pandda_processed_datasets_dirs)))

    print("Geting event table...")
    event_table = get_event_table(fs.pandda_inspect_events_path)
    print("\tFound {} PanDDA events".format(len(event_table)))

    print("Getting event models...")
    events = get_events(event_table,
                        fs,
                        )
    print("\tGot models of {} events".format(len(events)))

    print("Autobuilding...")
    autobuilding_results = map_parallel(autobuild_event,
                                        events,
                                        )
    print("\tAutobuilt {} events".format(len(autobuilding_results)))
    for result in autobuilding_results: print("\t{} {} : RSCC: {}".format(result.dtag, result.event_idx, result.rscc))

    print("Making autobuilding results table...")
    results_table = ResultsTable(autobuilding_results)
    print("\tMade autobuilding resutls table")

    print("Merging best models")
    merge_models(events,
                 results_table.table,
                 fs,
                 )

    print("Outputing autobuilding results table...")
    results_table.to_csv(fs.autobuilding_results_table)
    print("\tOutput autobuilding results table")


if __name__ == "__main__":
    main()
