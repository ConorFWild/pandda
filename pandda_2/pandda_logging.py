import logging
from functools import wraps
import json

import numpy as np


class PanDDALog:
    def __init__(self,
                 log_file,
                 name="PanDDA",
                 level=logging.DEBUG,
                 log_format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
                 ):
        self.log_file = log_file
        # self.log = logging.getLogger(name)
        # formtter = logging.Formatter(log_format)
        #
        # file_handler = logging.FileHandler(str(log_file))
        # file_handler.setFormatter(formtter)
        # file_handler.setLevel(level)
        #
        # stream_handler = logging.StreamHandler()
        # stream_handler.setFormatter(formtter)
        # stream_handler.setLevel(level)
        #
        # self.log.addHandler(file_handler)
        # self.log.addHandler(stream_handler)

    def __call__(self, logstr):
        print(logstr)
        # self.log.info(logstr)
        with open(str(self.log_file), "a+") as f:
            f.write(logstr)


def title(string):
    head = "########################################################\n"
    mid = "#\t" + string + "\n"
    tail = "########################################################\n"
    return head + mid + tail


def indent(string, level=1):
    for i in range(level):
        string = "\t" + string
    return string


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


# @logging_method
def dict_to_str(dict):
    logstr = str(json.dumps(dict, indent=4))
    return logstr


def log_startup(working_phil):
    logstr = ""
    logstr += title("Starting PanDDA")
    logstr += str(working_phil)
    return logstr


# @logging_method
def log_config(config):
    logstr = str(config)
    return logstr


# @logging_method
def log_load_datasets(datasets_dict):
    logstr = ""
    logstr += title("Loaded datasets")
    # logstr += dict_to_str(datasets_dict)
    logstr += str(datasets_dict)

    return logstr


# @logging_method
def log_get_reference(reference):
    logstr = title("Selected reference")
    return logstr


# @logging_method
def log_transform_dataset(dataset):
    logstr = title("Transformed dataset")
    return logstr


# @logging_method
def log_partitioning(dataset):
    logstr = title("partitioned dataset")
    # logstr += dict_to_str(dataset.partitions)
    logstr += str(dataset.partitions)
    return logstr


# '@logging_method
def log_get_grid(grid):
    logstr = title("Got grid")
    return logstr


# '@logging_method
def log_get_shells(shells):
    logstr = title("Got shells")
    for shell_num, shell_dataset in shells.items():
        logstr += indent(str(shell_dataset.get_partition("test")))
        logstr += indent(str(shell_dataset.get_partition("train")))
    return logstr


# @logging_method
def log_output(tree):
    logstr = title("Created output")
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


def log_shell_tasks(shell_processors):
    logstr = ""
    logstr += title("Prepared shell processors")
    logstr += indent("Processing {} shells...".format(len(shell_processors)))
    return logstr


def log_process_shells(event_tables):
    logstr = ""
    logstr += title("Completed processing shells!")
    return logstr


def log_make_event_table(event_table):
    logstr = ""
    logstr += title("Made events table")
    logstr += "Discovered {} events".format(len(event_table))
    return logstr


def log_sites_table(sites_table):
    logstr = ""
    logstr += title("Made sites table")
    logstr += "Discovered {} sites".format(len(sites_table))
    return logstr


def log_output_sites_table(sites_table,
                           tree,
                           ):
    logstr = ""
    logstr += title("Output sites table")
    logstr += "Sites table of len: {}".format(len(sites_table))
    return logstr


def log_finish_pandda(pandda_start_time,
                      pandda_finish_time,
                      ):
    logstr = ""
    logstr += title("PanDDA Finished sucessfully!")
    logstr += "PanDDA finished in {}".format(pandda_finish_time - pandda_start_time)
    return logstr
