from collections import OrderedDict


class PrePanDDAPartitions:

    def __init__(self, max_build_datasets, min_build_datasets):

        self.max_build_datasets = max_build_datasets
        self.min_build_datasets = min_build_datasets

        self.partitions = {}

    def __call__(self, partition, samples):
        return {dtag: samples
                for dtag, samples
                in samples.items()
                if dtag in self.partitions[partition]}

    def get_train(self, datasets, high_res_cutoff):

        self.partitions["train"] = self.select_datasets_for_density_characterisation(datasets, high_res_cutoff)
        return self.partitions["train"]

    def get_test(self, datasets, high_res_cutoff):

        self.partitions["test"] = self.select_datasets_for_density_characterisation(datasets, high_res_cutoff)
        return self.partitions["test"]

    def select_datasets_for_density_characterisation(self, datasets, high_res_cutoff):
        """Select all datasets with resolution better than high_res_cutoff"""

        building_datasets = []

        # Counter for the number of datasets to select
        total_build = 0
        # Select from the datasets that haven't been rejected
        for dtag, dataset in datasets.items():
            # Check the resolution of the dataset
            if dataset.data.summary.high_res > high_res_cutoff:
                building_datasets.append(dtag)

                # Check to see if the number of datasets to use in building has been reached
                total_build += 1
                if total_build >= self.max_build_datasets:
                    break

        return building_datasets

    def select_datasets_for_analysis(self, datasets, high_res_large_cutoff, high_res_small_cutoff):
        analysis_datasets = []

        # Select from the datasets that haven't been rejected
        for dtag, dataset in datasets.items():

            # Check the resolution of the dataset (is not too low)
            if dataset.data.summary.high_res > high_res_large_cutoff:
                continue
            # Check the resolution of the dataset (is not too high)
            elif dataset.data.summary.high_res <= high_res_small_cutoff:
                continue
            else:
                analysis_datasets.append(dataset.tag)

        return analysis_datasets


class DefaultPanDDAPartitionsGetter:

    def __init__(self, test, train, not_test, not_train):
        self.test = test
        self.train = train
        self.not_test = not_test
        self.not_train = not_train

    def __call__(self, datasets):
        partitions = OrderedDict()

        dtags = [dtag for dtag in datasets]

        # Get train set
        if self.train:
            # Get train datasets that are in datasets
            partitions["train"] = [dtag for dtag in self.train if (dtag in dtags)]
        else:
            partitions["train"] = dtags

        # Exclude non-training datasets
        if self.not_train:
            partitions["train"] = [dtag for dtag in partitions["train"] if not (dtag in self.not_train)]

        # Get test set
        if self.train:
            # Get train datasets that are in datasets
            partitions["test"] = [dtag for dtag in self.test if (dtag in dtags)]
        else:
            partitions["test"] = dtags

        # Exclude non-training datasets
        if self.not_train:
            partitions["test"] = [dtag for dtag in partitions["test"] if not (dtag in self.not_test)]

        return partitions


#
# class DefaultPanDDAPartitions:
#
#     def __init__(self):
#         self.partitions = OrderedDict()
#
#     def __call__(self, partition, samples):
#         return {dtag: samples
#                 for dtag, samples
#                 in samples.items()
#                 if dtag in self.partitions[partition]}
#
#     def set_partition(self, name, dtags):
#         self.partitions[name] = dtags
#
#     def get_partition(self, name):
#         return self.partitions[name]
#
#     def sample_partition(self, partition, samples):
#         return {dtag: samples
#                 for dtag, samples
#                 in samples.items()
#                 if dtag in self.partitions[partition]}

    # def get_train(self):
    #     return self.partitions["train"]
    #
    # def get_test(self):
    #     return self.partitions["test"]
    #
    # def set_train(self, train_dtags):
    #
    #     self.partitions["train"] = train_dtags
    #     return self.partitions["train"]
    #
    # def set_test(self, test_dtags):
    #
    #     self.partitions["test"] = test_dtags
    #     return self.partitions["test"]
