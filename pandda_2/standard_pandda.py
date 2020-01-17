from collections import OrderedDict


class StandardPandda:
    def __init__(self,
                 load_dataset,
                 get_reference,
                 transform_dataset,
                 get_grid,
                 partitioner,
                 output,
                 create_shells,
                 process_shell,
                 processor,
                 create_event_table,
                 output_event_table,
                 ):
        self.load_dataset = load_dataset
        self.get_reference = get_reference
        self.transform_dataset = transform_dataset
        self.get_grid = get_grid
        self.partitioner = partitioner
        self.output = output
        self.create_shells = create_shells
        self.processor = processor
        self.process_shell = process_shell
        self.create_event_table = create_event_table
        self.output_event_table = output_event_table

    def __call__(self):
        print("Loading dataset")
        dataset = self.load_dataset()

        print("Loading reference")
        reference = self.get_reference(dataset.datasets)

        print("Transforming dataset")
        dataset = self.transform_dataset(dataset, reference)

        print("Partitioning")
        dataset.partitions = self.partitioner(dataset.datasets)

        print("Getting grid")
        grid = self.get_grid(reference)

        print("Partitioning shells")
        shells = {shell_num: shell_dataset
                  for shell_num, shell_dataset
                  in self.create_shells(dataset)
                  }

        print("Output")
        tree = self.output(dataset,
                           shells,
                           )

        print("Processing shells")
        shell_processors = []
        for shell_num, shell_dataset in shells.items():
            shell_p = TaskWrapper(self.process_shell,
                                  shell_dataset=shell_dataset,
                                  reference=reference,
                                  grid=grid,
                                  tree=tree,
                                  shell_num=shell_num,
                                  )
            shell_processors.append(shell_p)
        event_tables = self.processor(shell_processors)

        print("Creating event table")
        event_table = self.create_event_table(tree,
                                              len(shells),
                                              )

        print("Outputting event table")
        self.output_event_table(event_table,
                                tree["analyses"]["pandda_analyse_events"](),
                                )

    def repr(self):
        repr = OrderedDict()
        repr["load_dataset"] = self.load_dataset.__repr__()
        repr["transform_dataset"] = self.transform_dataset.__repr__()
        repr["get_reference"] = self.get_reference.__repr__()
        repr["get_grid"] = self.get_grid.__repr__()
        repr["create_shells"] = self.create_shells.__repr__()
        repr["output"] = self.output.__repr__()
        repr["process_shells"] = self.process_shell.repr()
        repr["processor"] = self.processor.repr()
        repr["create_event_table"] = self.create_event_table.__repr__()

        return repr


class TaskWrapper:
    def __init__(self, func, *args, **kwargs):
        self.func = func
        self.args = args
        self.kwargs = kwargs

    def __call__(self):
        return self.func(*self.args, **self.kwargs)
