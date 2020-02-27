import pathlib as p


class Criticiser:

    def __init__(self):
        pass

    def __call__(self, tree, map_maker, map_loader, truncated_dataset, ref_map, events, bdcs, grid):
        dataset_path = p.Path(tree(("processed_datasets", truncated_dataset.tag))[0])

        map_maker.process_single(map_loader, truncated_dataset, ref_map, events, bdcs, dataset_path, grid)


class CriticiserAll:

    def __init__(self, map_loader, event_table_maker, map_maker):
        self.map_loader = map_loader
        self.event_table_maker = event_table_maker
        self.map_maker = map_maker

    def __call__(self, map_resolution, reference_map, grid, tree, bdcs, events, ref_map,
                 truncated_dataset, statistical_model, name, dataset):

        # Produce maps that are shared by iteration
        dir_path = p.Path(tree([str(name)])[0])
        dir_path_string = str(dir_path)
        self.map_maker.process_shell(self.map_loader,
                                     truncated_dataset,
                                     ref_map,
                                     events,
                                     bdcs,
                                     dir_path_string,
                                     grid,
                                     statistical_model)

        # Produce the event table
        event_table_path = dir_path / "event_table.csv"

        event_table = self.event_table_maker(dataset.partition_datasets("test"),
                                             events,
                                             event_table_path)

        return event_table
