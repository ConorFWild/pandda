from __future__ import print_function

# # Imports
import sys, time
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


if __name__ == "__main__":

    # Get start time
    pandda_start_time = time.time()

    # Parse Config files and command line arguments
    working_phil = config.extract_params_default(master_phil=pandda_phil.pandda_phil,
                                                 args=sys.argv[1:],
                                                 blank_arg_prepend=None,
                                                 home_scope=None
                                                 ).extract()

    # Maps options to code abstraction: Phil -> Config
    config = config.Config(working_phil)

    # Options: maps a config to options
    options = options.Options(config)
    checks.check_config(config)

    # Define program
    dataset = options.load_dataset()
    pandda_logging.log_load_datasets(dataset.datasets)

    reference = options.get_reference(dataset.datasets)
    pandda_logging.log_get_reference(reference)

    dataset = options.transform_dataset(dataset, reference)
    pandda_logging.log_transform_dataset(dataset)

    dataset.partitions = options.partitioner(dataset.datasets)
    pandda_logging.log_partitioning(dataset)

    grid = options.get_grid(reference)
    pandda_logging.log_get_grid(grid)

    shells = {shell_num: shell_dataset
              for shell_num, shell_dataset
              in options.create_shells(dataset)
              }
    pandda_logging.log_get_shells(shells)

    tree = options.output(dataset,
                          shells,
                          )
    pandda_logging.log_output(tree)

    print("Processing shells")
    shell_processors = []
    for shell_num, shell_dataset in shells.items():

        shell_p = TaskWrapper(options.process_shell,
                              shell_dataset=shell_dataset,
                              reference=reference,
                              grid=grid,
                              tree=tree,
                              shell_num=shell_num,
                              )
        shell_processors.append(shell_p)
    event_tables = options.processer(shell_processors,
                                     output_paths=[tree["shells"][shell_num]["event_table"]()
                                                   for shell_num
                                                   in range(len(shell_processors))
                                                   ],
                                     result_loader=None,
                                     shared_tmp_dir=tree["shells"](),
                                     )

    print("Creating event table")
    event_table = options.create_event_table(tree,
                                             len(shells),
                                             )

    print(event_table)

    print("Creating sites table")
    sites_table, events_table_with_sites = options.create_sites_table(event_table,
                                                                      grid,
                                                                      reference,
                                                                      )

    print(events_table_with_sites)

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

    # exit()

    print("Outputting event table")
    options.output_sites_table(sites_table,
                               tree["analyses"]["pandda_analyse_sites"](),
                               )

    print("Outputting event table")
    options.output_event_table(events_table_with_sites,
                               tree["analyses"]["pandda_analyse_events"](),
                               )

    pandda_finish_time = time.time()
    print("PanDDA ran in: {}".format(pandda_finish_time - pandda_start_time))
