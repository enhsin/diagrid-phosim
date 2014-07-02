#!/usr/bin/env bash

bin/phosim.py validation/validation_4D_catalog -c validation/validation_4D_00_commands -e 0 "$@"
bin/phosim.py validation/validation_4D_catalog -c validation/validation_4D_01_commands -e 0 "$@"


