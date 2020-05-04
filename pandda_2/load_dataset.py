import multi_dataset_crystalography


class LoadDataset:
    """
    A class for loading a multicrystal dataset
    dataloader: () -> Dict[Dataset]
    sample_loader: TODO: DEPRECATED
    """
    def __init__(self,
                 dataloader,
                 ):
        self.dataloader = dataloader

    def __call__(self):
        pandda_dataset = multi_dataset_crystalography.dataset.dataset.MultiCrystalDatasetPlain(dataloader=self.dataloader,
                                                                                          )

        return pandda_dataset

    def __repr__(self):
        repr = {"dataloader": self.dataloader.__repr__(),
                "sample_loader": self.sample_loader.__repr__(),
                }
        return repr