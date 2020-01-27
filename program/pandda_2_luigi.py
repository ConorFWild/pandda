from __future__ import print_function

# # Imports
import sys, time
import pprint
import json
from pandda_2 import (config,
                      pandda_phil,
                      options,
                      load_dataset,
                      transform_dataset,
                      get_dataset,
                      get_reference,
                      get_grid,
                      define_tree,
                      make_tree,
                      copy_dataset_files,
                      output,
                      create_shells,
                      map_loader,
                      get_reference_map,
                      statistical_model,
                      fit_model,
                      evaluate_model,
                      cluster_outliers,
                      filter_clusters,
                      event_analyser,
                      make_event_map,
                      make_z_map,
                      make_mean_map,
                      make_event_table,
                      process_shell,
                      process_shells,
                      processor,
                      create_sites_table,
                      output_sites_table,
                      create_event_table,
                      output_event_table,
                      standard_pandda,
                      autobuild,
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

    # Get dataset loader
    load_dataset = load_dataset.LoadDataset(dataloader=options.dataloader,
                                            sample_loader=options.sample_loader,
                                            )

    # Get reference loader
    get_reference = options.get_reference

    # Get dataset transformer
    transform_dataset = transform_dataset.TransformDataset(transform_data_check=options.transform_data_check,
                                                           transform_scale_diffraction=options.transform_scale_diffraction,
                                                           transform_filter_structure=options.transform_filter_structure,
                                                           transform_filter_wilson=options.transform_filter_wilson,
                                                           transform_align=options.transform_align,
                                                           )

    # Get grid loader
    get_grid = options.grid_loader

    # Get partitioner
    partitioner = options.partitioner

    # Get output handler
    define_tree = define_tree.DefineTree(output_dir=config.output.out_dir)
    make_tree = make_tree.MakeTree(overwrite=config.output.overwrite)
    copy_dataset_files = copy_dataset_files.DatasetFileCopier()

    output = output.Output(define_tree=define_tree,
                           make_tree=make_tree,
                           copy_dataset_files=copy_dataset_files,
                           )

    # Get resolution shell scheme
    create_shells = create_shells.CreateShells(min_train_datasets=60,
                                               max_test_datasets=60,
                                               cutoff=0.1,
                                               )

    # Get Resolution shell processor
    diffraction_data_truncator = options.diffraction_data_truncator
    reference_map_getter = options.reference_map_getter
    get_reference_map = get_reference_map.GetReferenceMap()

    map_loader = options.map_loader
    statistical_model = options.statistical_model
    fit_model = fit_model.FitModel()
    evaluate_model = evaluate_model.EvaluateModel()
    cluster_outliers = options.clusterer
    filter_clusters = options.event_finder
    event_analyser = options.event_analyser

    make_event_map = make_event_map.MakeEventMap()
    make_z_map = make_z_map.MakeZMap()
    make_mean_map = make_mean_map.MakeMeanMap()
    make_event_table = make_event_table.MakeEventTable()

    process_in_shell = processor.ProcessorDictJoblib()

    process_shell = process_shell.ProcessShell(diffraction_data_truncator=diffraction_data_truncator,
                                               reference_map_getter=reference_map_getter,
                                               get_reference_map=get_reference_map,
                                               map_loader=map_loader,
                                               statistical_model=statistical_model,
                                               fit_model=fit_model,
                                               evaluate_model=evaluate_model,
                                               cluster_outliers=cluster_outliers,
                                               filter_clusters=filter_clusters,
                                               event_analyser=event_analyser,
                                               make_event_map=make_event_map,
                                               make_z_map=make_z_map,
                                               make_mean_map=make_mean_map,
                                               make_event_table=make_event_table,
                                               process=process_in_shell,
                                               )

    processer = processor.ProcessorLuigi(jobs=10,
                                         parallel_env="smp",
                                         n_cpu=12,
                                         run_locally=False,
                                         )

    # Get site table creator
    create_sites_table = create_sites_table.CreateSitesTable()

    # Get event table outputter
    output_sites_table = output_sites_table.OutputSitesTable()

    # Get event table processor
    create_event_table = create_event_table.CreateEventTable()

    # Get event table outputter
    output_event_table = output_event_table.OutputEventTable()

    # Autobuilders
    autobuilder = autobuild.AutobuildQFit()

    # Define program
    print("Loading dataset")
    dataset = load_dataset()

    print("Loading reference")
    reference = get_reference(dataset.datasets)

    print("Transforming dataset")
    dataset = transform_dataset(dataset, reference)

    print("Partitioning")
    dataset.partitions = partitioner(dataset.datasets)

    print("Getting grid")
    grid = get_grid(reference)

    print("Partitioning shells")
    shells = {shell_num: shell_dataset
              for shell_num, shell_dataset
              in create_shells(dataset)
              }

    for shell_num, shell_dataset in shells.items():
        print(shell_dataset.get_partition("test"))
        print(shell_dataset.get_partition("train"))

    print("Output")
    tree = output(dataset,
                  shells,
                  )

    print("Processing shells")
    shell_processors = []
    for shell_num, shell_dataset in shells.items():
        # if shell_num == 5:
        #     process_shell(shell_dataset = shell_dataset,
        #                   reference = reference,
        #                   grid = grid,
        #                   tree = tree,
        #                   shell_num = shell_num,
        #                   )


        shell_p = TaskWrapper(process_shell,
                              shell_dataset=shell_dataset,
                              reference=reference,
                              grid=grid,
                              tree=tree,
                              shell_num=shell_num,
                              )
        shell_processors.append(shell_p)
    event_tables = processer(shell_processors,
                             output_paths=[tree["shells"][shell_num]["event_table"]()
                                           for shell_num
                                           in range(len(shell_processors))
                                           ],
                             result_loader=None,
                             shared_tmp_dir=tree["shells"](),
                             )

    print("Creating event table")
    event_table = create_event_table(tree,
                                     len(shells),
                                     )

    print(event_table)

    print("Creating sites table")
    sites_table, events_table_with_sites = create_sites_table(event_table,
                                                              grid,
                                                              reference,
                                                              )

    print(events_table_with_sites)

    # print("Autobuilding")
    # autobuilders = []
    # for index, event in events_table_with_sites.iterrows():
    #     dtag = event["dtag"]
    #     analysed_resolution = event["analysed_resolution"]
    #     bdc = event["1-BDC"]
    #     event_idx = event["event_idx"]
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

    print("Outputting event table")
    output_sites_table(sites_table,
                       tree["analyses"]["pandda_analyse_sites"](),
                       )

    print("Outputting event table")
    output_event_table(events_table_with_sites,
                       tree["analyses"]["pandda_analyse_events"](),
                       )

    pandda_finish_time = time.time()
    print("PanDDA ran in: {}".format(pandda_finish_time - pandda_start_time))
