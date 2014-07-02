#!/usr/bin/env bash

bin/phosim.py validation/validation_4B_00_catalog -c validation/validation_4B_00_commands -e 0 "$@"
bin/phosim.py validation/validation_4B_01_catalog -c validation/validation_4B_00_commands -e 0 "$@"
bin/phosim.py validation/validation_4B_02_catalog -c validation/validation_4B_00_commands -e 0 "$@"
bin/phosim.py validation/validation_4B_03_catalog -c validation/validation_4B_00_commands -e 0 "$@"
bin/phosim.py validation/validation_4B_04_catalog -c validation/validation_4B_00_commands -e 0 "$@"
bin/phosim.py validation/validation_4B_05_catalog -c validation/validation_4B_00_commands -e 0 "$@"
bin/phosim.py validation/validation_4B_00_catalog -c validation/validation_4B_01_commands -e 0 "$@"
bin/phosim.py validation/validation_4B_01_catalog -c validation/validation_4B_01_commands -e 0 "$@"
bin/phosim.py validation/validation_4B_02_catalog -c validation/validation_4B_01_commands -e 0 "$@"
bin/phosim.py validation/validation_4B_03_catalog -c validation/validation_4B_01_commands -e 0 "$@"
bin/phosim.py validation/validation_4B_04_catalog -c validation/validation_4B_01_commands -e 0 "$@"
bin/phosim.py validation/validation_4B_05_catalog -c validation/validation_4B_01_commands -e 0 "$@"

