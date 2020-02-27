from functools import wraps
import json

import numpy as np


def logging_method(f, method="print"):
    @wraps(f)
    def wrapper(*args, **kwargs):
        result = f(*args, **kwargs)
        if method == "print":
            print(result)
        else:
            raise Exception("Invalid logging options: {}".format(method))
        return

    return wrapper


@logging_method
def dict_to_str(dict):
    logstr = str(json.dumps(dict, indent=4))
    return logstr


@logging_method
def log_config(config):
    logstr = str(config)
    return logstr


@logging_method
def log_load_datasets(datasets_dict):
    logstr = "Loaded datasets"
    # logstr += dict_to_str(datasets_dict)
    logstr += str(datasets_dict)

    return logstr


@logging_method
def log_get_reference(reference):
    logstr = "Selected reference"
    return logstr


@logging_method
def log_transform_dataset(dataset):
    print("Transformed dataset")


@logging_method
def log_partitioning(dataset):
    logstr = "partitioned dataset\n"
    # logstr += dict_to_str(dataset.partitions)
    logstr += str(dataset.partitions)
    return logstr


@logging_method
def log_get_grid(grid):
    logstr = "Got grid"
    return logstr


@logging_method
def log_get_shells(shells):
    logstr = "Got shells"
    for shell_num, shell_dataset in shells.items():
        logstr += str(shell_dataset.get_partition("test"))
        logstr += str(shell_dataset.get_partition("train"))
    return logstr


@logging_method
def log_output(tree):
    logstr = "Created output"
    return logstr


def log_map_making(sample,
                   ref_map,
                   bdc,
                   ):
    import numpy as np

    logstr = "OUTPUTTING EVENT MAP!\n"
    logstr += "{}\n".format(sample)
    logstr += "{}\n".format(ref_map)
    logstr += "{}\n".format(bdc)

    logstr += "{}\n".format(np.mean(sample.data))
    logstr += "{}\n".format(np.mean(ref_map.data))

    logstr += "{}\n".format(np.min(sample.data))
    logstr += "{}\n".format(np.min(ref_map.data))

    logstr += "{}\n".format(np.max(sample.data))
    logstr += "{}\n".format(np.max(ref_map.data))

    print(logstr)

    return logstr


def log_summarise_map(electron_density_map):
    data = electron_density_map.data

    logstr = "summarising map!\n"
    logstr += "{}\n".format(np.mean(data))
    logstr += "{}\n".format(np.std(data))
    logstr += "{}\n".format(np.max(data))
    logstr += "{}\n".format(np.min(data))

    return logstr
