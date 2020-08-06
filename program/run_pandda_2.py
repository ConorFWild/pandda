from __future__ import print_function

# Imports
import sys
import os
import time
from pathlib import Path
from pandda_2 import (config,
                      pandda_phil,
                      options,
                      checks,
                      pandda_logging,
                      task_wrapper,
                      store,
                      )


def main():
    ####################################################################################################################
    # Program Setup
    ####################################################################################################################
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
    try:
        os.mkdir(str(pandda_config.output.out_dir))
    except:
        pass
    store.dump_config_to_json(pandda_config,
                              Path(str(pandda_config.output.out_dir)) / "pandda.json",
                              )

    # Setup logging
    log = pandda_logging.PanDDALog(log_file=Path(str(pandda_config.output.out_dir)) / "log.txt")
    log(pandda_logging.log_startup(working_phil))
    log(pandda_logging.log_config(pandda_config))

    ####################################################################################################################
    # Program
    ####################################################################################################################

    try:
        # Options: maps a config to code abstraction
        pandda_options = options.Options(pandda_config)
        checks.check_config(pandda_config)

        # Load the datasets
        dataset = pandda_options.load_dataset()
        log(pandda_logging.log_load_datasets(dataset.datasets))

        # Partition the dataset
        dataset.partitions = pandda_options.partitioner(dataset.datasets)
        log(pandda_logging.log_partitioning(dataset))

        # Get the reference
        reference = pandda_options.get_reference(dataset.partition_datasets("train"))
        log(pandda_logging.log_get_reference(reference))

        # Transform the dataset: check data, filter by RMSD, filter by wilson, scale diffraction, align
        dataset = pandda_options.transform_dataset(dataset, reference)
        log(pandda_logging.log_transform_dataset(dataset))

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
            shell_p = task_wrapper.TaskWrapper(pandda_options.process_shell,
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
        log(pandda_logging.add_shell_logs([tree["shells"][shell_num]() / "shell_log.txt"
                                           for shell_num
                                           in range(len(shell_processors))
                                           ]
                                          )
            )

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

    ####################################################################################################################
    # Program Wide error catching
    ####################################################################################################################
    except Exception as e:
        # Catch and log any error
        log(pandda_logging.error(e))


if __name__ == "__main__":
    main()
