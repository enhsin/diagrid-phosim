#!/usr/bin/env bash

bin/phosim.py validation/validation_1B_catalog -c validation/validation_1B_00_commands -e 0 "$@"
bin/phosim.py validation/validation_1B_catalog -c validation/validation_1B_01_commands -e 0 "$@"
bin/phosim.py validation/validation_1B_catalog -c validation/validation_1B_02_commands -e 0 "$@"
bin/phosim.py validation/validation_1B_catalog -c validation/validation_1B_03_commands -e 0 "$@"
bin/phosim.py validation/validation_1B_catalog -c validation/validation_1B_04_commands -e 0 "$@"
bin/phosim.py validation/validation_1B_catalog -c validation/validation_1B_05_commands -e 0 "$@"
bin/phosim.py validation/validation_1B_catalog -c validation/validation_1B_06_commands -e 0 "$@"
bin/phosim.py validation/validation_1B_catalog -c validation/validation_1B_07_commands -e 0 "$@"
bin/phosim.py validation/validation_1B_catalog -c validation/validation_1B_08_commands -e 0 "$@"

