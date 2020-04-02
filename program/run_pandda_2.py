from __future__ import print_function
import logging

# # Imports
import sys, time
from pathlib import Path
import pprint
import json
from pandda_2 import (config,
                      pandda_phil,
                      options,
                      process_shell,
                      checks,
                      pandda_logging,

                      )


class TaskWrapper:
    def __init__(self, func, *args, **kwargs):
        self.func = func
        self.args = args
        self.kwargs = kwargs

    def __call__(self):
        return self.func(*self.args, **self.kwargs)


def dump_config_to_json(pandda_config ,
                        output_path,
                        ):
    import json

    record = {"data_dirs": str(pandda_config .input.data_dirs),
              "out_dir": str(pandda_config .output.out_dir),
              }

    json_string = json.dumps(record)

    with open(str(output_path), "w") as f:
        f.write(json_string)


def main():
    # Get start time
    pandda_start_time = time.time()

    # Parse Config files and command line arguments
    working_phil = config.extract_params_default(master_phil=pandda_phil.pandda_phil,
                                                 args=sys.argv[1:],
                                                 blank_arg_prepend=None,
                                                 home_scope=None
                                                 ).extract()

    # Maps options to code abstraction: Phil -> Config
    pandda_config = config.Config(working_phil)
    dump_config_to_json(pandda_config ,
                        Path(str(pandda_config .output.out_dir)) / "pandda.json",
                        )
    log = pandda_logging.PanDDALog(log_file=Path(str(pandda_config .output.out_dir)) / "log.txt")
    log(pandda_logging.log_startup(working_phil))
    log(pandda_logging.log_config(pandda_config ))
    # try:

    # Options: maps a config to code abstraction
    pandda_options = options.Options(pandda_config )
    checks.check_config(pandda_config )

    # Load the datase
    dataset = pandda_options.load_dataset()
    log(pandda_logging.log_load_datasets(dataset.datasets))

    # Get the reference
    reference = pandda_options.get_reference(dataset.datasets)
    log(pandda_logging.log_get_reference(reference))

    # Transform the dataset: check data, filter by RMSD, filter by wilson, scale diffraction, align
    dataset = pandda_options.transform_dataset(dataset, reference)
    log(pandda_logging.log_transform_dataset(dataset))

    # Partition the dataset
    dataset.partitions = pandda_options.partitioner(dataset.datasets)
    log(pandda_logging.log_partitioning(dataset))

    # Generate the grid
    grid = pandda_options.get_grid(reference)
    log(pandda_logging.log_get_grid(grid))

    # Define the resolution shells to work on
    shells = {shell_num: shell_dataset
              for shell_num, shell_dataset
              in pandda_options.create_shells(dataset)
              }
    log(pandda_logging.log_get_shells(shells))

    # Generate the output tree
    tree = pandda_options.output(dataset,
                          shells,
                          )
    log(pandda_logging.log_output(tree))

    # Define the resolution shell processing tasks
    shell_processors = []
    for shell_num, shell_dataset in shells.items():
        shell_p = TaskWrapper(pandda_options.process_shell,
                              shell_dataset=shell_dataset,
                              reference=reference,
                              grid=grid,
                              tree=tree,
                              shell_num=shell_num,
                              )
        shell_processors.append(shell_p)
    log(pandda_logging.log_shell_tasks(shell_processors))

    # Process the resolution shells
    event_tables = pandda_options.processer(shell_processors,
                                     output_paths=[tree["shells"][shell_num]["event_table"]()
                                                   for shell_num
                                                   in range(len(shell_processors))
                                                   ],
                                     result_loader=None,
                                     shared_tmp_dir=tree["shells"](),
                                     )
    log(pandda_logging.log_process_shells(event_tables))

    # Create the pandda event table
    event_table = pandda_options.create_event_table(tree,
                                             len(shells),
                                             )

    log(pandda_logging.log_make_event_table(event_table))

    # Create the site table
    sites_table, events_table_with_sites = pandda_options.create_sites_table(event_table,
                                                                      grid,
                                                                      reference,
                                                                      )
    log(pandda_logging.log_sites_table(sites_table))

    # Output the site table
    pandda_options.output_sites_table(sites_table,
                               tree["analyses"]["pandda_analyse_sites"](),
                               )
    log(pandda_logging.log_output_sites_table(sites_table,
                                              tree["analyses"]["pandda_analyse_sites"](),
                                              )
        )

    # Output the event table
    pandda_options.output_event_table(events_table_with_sites,
                               tree["analyses"]["pandda_analyse_events"](),
                               )

    pandda_finish_time = time.time()
    log(pandda_logging.log_finish_pandda(pandda_start_time,
                                         pandda_finish_time,
                                         )
        )
    # except Exception as e:
    #     log(str(e))


if __name__ == "__main__":
    main()

# print("Autobuilding")
# autobuilders = {}
# for index, event in events_table_with_sites.iterrows():
#
#     print("index: {}".format(index))
#     print("event: {}".format(event))
#     dtag = index[0]
#     analysed_resolution = event["analysed_resolution"]
#     bdc = event["1-BDC"]
#     event_idx = int(float(index[1]))
#
#
#
#     autobuilders[event_idx] = autobuilder(
#         protein_model_path=tree["processed_datasets"][dtag]["initial_model"](),
#         ligand_model_path=tree["processed_datasets"][dtag]["ligand_files"]().glob("*.pdb").next(),  # TODO: fix
#         event_map_path=tree["processed_datasets"][dtag]["event_map"]([dtag,
#                                                                       event_idx,
#                                                                       bdc,
#                                                                       ],
#                                                                      ),
#         output_dir_path=tree["processed_datasets"][dtag]["autobuilt"](),
#         resolution=analysed_resolution,
#         event=event,
#     )
#     print(tree["processed_datasets"][dtag]["initial_model"](),
#           tree["processed_datasets"][dtag]["ligand_files"]().glob("*.pdb").next(),
#           tree["processed_datasets"][dtag]["autobuilt"](),
#           analysed_resolution,
#           event,
#           )

#     autobuild_event = TaskWrapper(autobuilder,
#                                   protein_model_path=tree["processed_datasets"][dtag]["initial_model"](),
#                                   ligand_model_path=tree["processed_datasets"][dtag]["ligand_files"]().glob("*.pdb").next(),  # TODO: fix
#                                   event_map_path=tree["processed_datasets"][dtag]["event_map"]([dtag,
#                                                                                                event_idx ,
#                                                                                                bdc,
#                                                                                                 ],
#                                                                                                ),
#                                   output_dir_path=tree["processed_datasets"][dtag]["autobuilt"](),
#                                   resolution=analysed_resolution,
#                                   event=event,
#                                   )
#     autobuilders.append(autobuild_event)
# autobuilt = processer(autobuilders,
#                          output_paths=[tree["processed_datasets"][dtag]["autobuilt"]()
#                                        for dtag
#                                        in events_table_with_sites["dtag"]
#                                        ],
#                          result_loader=None,
#                          shared_tmp_dir=tree["processed_datasets"](),
#                          )

# autobuilt = process_in_shell(autobuilders)
