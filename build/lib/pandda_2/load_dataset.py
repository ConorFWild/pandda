import multi_dataset_crystalography


class LoadDataset:
    """
    A class for loading a multicrystal dataset
    dataloader: () -> Dict[Dataset]
    sample_loader: TODO: DEPRECATED
    """
    def __init__(self,
                 dataloader,
                 sample_loader,
                 ):
        self.dataloader = dataloader
        self.sample_loader = sample_loader

    def __call__(self):
        pandda_dataset = multi_dataset_crystalography.dataset.dataset.MultiCrystalDataset(dataloader=self.dataloader,
                                                                                          sample_loader=self.sample_loader,
                                                                                          )

        return pandda_dataset

    def __repr__(self):
        repr = {"dataloader": self.dataloader.__repr__(),
                "sample_loader": self.sample_loader.__repr__(),
                }
        return repr