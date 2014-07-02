#!/usr/bin/env bash

bin/phosim.py validation/validation_2F_00_catalog -c validation/validation_2F_commands -e 0 "$@"
bin/phosim.py validation/validation_2F_01_catalog -c validation/validation_2F_commands -e 0 "$@"
bin/phosim.py validation/validation_2F_02_catalog -c validation/validation_2F_commands -e 0 "$@"
bin/phosim.py validation/validation_2F_03_catalog -c validation/validation_2F_commands -e 0 "$@"
bin/phosim.py validation/validation_2F_04_catalog -c validation/validation_2F_commands -e 0 "$@"
bin/phosim.py validation/validation_2F_05_catalog -c validation/validation_2F_commands -e 0 "$@"

