#!/usr/bin/env bash

bin/phosim.py validation/validation_4C_catalog -c validation/validation_4C_00_commands -e 0 "$@"
bin/phosim.py validation/validation_4C_catalog -c validation/validation_4C_01_commands -e 0 "$@"

