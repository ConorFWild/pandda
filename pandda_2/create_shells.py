from collections import OrderedDict

class CreateShells:

    def __init__(self, min_train_datasets=60, max_test_datasets=60, cutoff=0.1):
        self.min_train_datasets = min_train_datasets
        self.max_test_datasets = max_test_datasets
        self.cutoff = cutoff

    def __call__(self, pandda_dataset):

        # Get all acceptable training datasets
        all_train_datasets = pandda_dataset.partition_datasets("train")
        print("Got all train datasets")

        # Sort dataset tags by resolution
        all_train_datasets_sorted = sorted(all_train_datasets.items(), key=lambda kv: kv[1].data.summary.high_res)
        print("Sorted datasets")

        # Accumulate datasets until number of characterisation datasets reached
        # The upper bound on this is the lowest res shell possible
        train_datasets = {}
        train_res_max = 0
        while len(train_datasets) < self.min_train_datasets:
            train_dataset = all_train_datasets_sorted.pop(0)
            print(train_dataset)
            train_datasets[train_dataset[0]] = train_dataset[1]
            train_res_max = train_dataset[1].data.summary.high_res
        print("Collected trian datasets")

        # Get all acceptable training datasets
        all_test_datasets = pandda_dataset.partition_datasets("test")
        print("Got all test datasets")

        # Sort dataset tags by resolution
        all_test_datasets_sorted = sorted(all_test_datasets.items(), key=lambda kv: kv[1].data.summary.high_res)
        print("sorted test datasets")

        # Accumulate datasets up to the max - per - iteration and then yield
        test_datasets = []
        test_res_max = 0
        dataset_res_max = train_res_max
        idx = 0
        for test_dataset in all_test_datasets_sorted:

            test_res_current = test_dataset[1].data.summary.high_res

            # If bwould not go outside batch res bounds or make the dataset too large, add, update res and continue
            if ((test_res_current - dataset_res_max) < self.cutoff) and (len(test_datasets) < int(self.max_test_datasets)):
                test_datasets.append(test_dataset)
                test_res_max = test_dataset[1].data.summary.high_res
                continue

            # Else yield the dataset for processing
            print("yielding dataset")

            # Convert test datasets to dict
            test_datasets_dict = {key: value
                                  for key, value
                                  in test_datasets
                                  }

            # Update test datasets with train datasets
            new_dataset = {key: value for key, value in test_datasets}
            new_dataset.update(train_datasets)

            print(new_dataset)

            # Make new mcd with these datasets
            dataset = pandda_dataset.new_from_datasets(datasets=new_dataset)

            dataset.set_partition("test",
                                  list(test_datasets_dict.keys()),
                                  )
            dataset.set_partition("train",
                                  list(train_datasets.keys()),
                                  )

            # yield the dataset for processing
            print("Dataset {} of length {}; res limits ({},{})"
                  .format(idx,
                          len(dataset.datasets),
                          dataset_res_max,
                          test_datasets[-1][1].data.summary.high_res))

            yield idx, dataset

            # Empty old list
            test_datasets = []

            # Continue adding
            test_datasets.append(test_dataset)

            # Update dataset batch upper bound
            dataset_res_max = test_res_current

            # update idx
            idx = idx + 1

    def __repr__(self):
        repr = OrderedDict()
        repr["min_train_datasets"] = self.min_train_datasets
        repr["max_test_datasets"] = self.max_test_datasets
        repr["cutoff"] = self.cutoff
        return repr