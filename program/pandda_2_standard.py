from __future__ import print_function

# # Imports
import sys
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
                      create_event_table,
                      output_event_table,
                      standard_pandda,
                      )

if __name__ == "__main__":
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

    output = output.Output(define_tree=define_tree,
                           make_tree=make_tree,
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

    processer = processor.Processor()

    # Get event table processor
    create_event_table = create_event_table.CreateEventTable()

    # Get event table outputter
    output_event_table = output_event_table.OutputEventTable()

    # Define program
    pandda = standard_pandda.StandardPandda(load_dataset=load_dataset,
                                            get_reference=get_reference,
                                            transform_dataset=transform_dataset,
                                            get_grid=get_grid,
                                            partitioner=partitioner,
                                            output=output,
                                            create_shells=create_shells,
                                            process_shell=process_shell,
                                            processor=processer,
                                            create_event_table=create_event_table,
                                            output_event_table=output_event_table,
                                            )

    ####################################################

    # Summarise the PanDDA to be performed
    print("pandda summary")
    print(json.dumps(pandda.repr(),
                     indent=8,
                     )
          )

    # run the pandda
    pandda()
