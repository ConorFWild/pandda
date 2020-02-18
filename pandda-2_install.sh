#!/bin/bash
source /ccp4/bin/ccp4.setup-sh
cd /pandda
ccp4-python -m pip install --upgrade .
ccp4-python -m pip --no-cache-dir install --global-option=build_ext --global-option="-I/usr/include/python2.7 -lpython2.7" biopandas==0.2.4
