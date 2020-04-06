import json

def dump_config_to_json(pandda_config,
                        output_path,
                        ):
    record = {"data_dirs": str(pandda_config.input.data_dirs),
              "out_dir": str(pandda_config.output.out_dir),
              }

    json_string = json.dumps(record)

    with open(str(output_path), "w") as f:
        f.write(json_string)