class GetDataset:
    def __init__(self,
                 load_dataset,
                 transform_dataset,
                 ):
        self.load_dataset = load_dataset
        self.transform_dataset = transform_dataset

    def __call__(self,):
        dataset = self.load_dataset()
        dataset = self.transform_dataset(dataset)
        return dataset

    def __repr__(self):

        repr = {"load_dataset": self.load_dataset.__repr__(),
                "transform_dataset": self.transform_dataset.__repr__(),
                }
        return repr