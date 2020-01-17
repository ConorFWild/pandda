class Output:
    def __init__(self,
                 define_tree,
                 make_tree,
                 copy_dataset_files,
                 ):
        self.define_tree = define_tree
        self.make_tree = make_tree
        self.copy_dataset_files = copy_dataset_files

    def __call__(self,
                 dataset,
                 shells,
                 ):
        tree = self.define_tree(dataset,
                                shells,
                                )
        self.make_tree(tree)

        self.copy_dataset_files(dataset,
                                tree,
                                )
        return tree

    def __repr__(self):
        repr = {}
        return repr
