from collections import OrderedDict


class DefaultPanDDAPartitionsGetter:

    def __init__(self, test, train,):
        self.test = test
        self.train = train

    def __call__(self, datasets):
        partitions = OrderedDict()

        dtags = [dtag for dtag in datasets]

        # Get train set
        if bool(self.train):
            # Get train datasets that are in datasets
            partitions["train"] = [dtag for dtag in self.train if (dtag in dtags)]
        else:
            partitions["train"] = [dtag for dtag in dtags]

        # Get test set
        if bool(self.test):
            partitions["test"] = [dtag for dtag in self.test if (dtag in dtags)]
        else:
            partitions["test"] = [dtag for dtag in dtags]

        return partitions
