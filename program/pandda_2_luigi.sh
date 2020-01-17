#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

python $DIR/../experiments/pandda_2_luigi/pandda_2_luigi.py "$@"

