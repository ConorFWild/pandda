from pathlib import Path


def check_config(config):
    check_out_dir(config.output)

    check_in_dir(config.input)


def check_out_dir(output):
    out_dir = Path(output.out_dir)

    if out_dir.parent.exists():
        pass
    else:
        raise Exception("The output directory path {} does not exists".format(out_dir.parent))



def check_in_dir(input):
    in_dir = Path(input.data_dirs)

    if in_dir.parent.exists():
        pass
    else:
        raise Exception("The input data dirs {} do not exist".format(in_dir.parent))
