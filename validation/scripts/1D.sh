#!/usr/bin/env bash

bin/phosim.py validation/validation_1D_00_catalog -c validation/validation_1D_commands -e 0 "$@"
bin/phosim.py validation/validation_1D_01_catalog -c validation/validation_1D_commands -e 0 "$@"
bin/phosim.py validation/validation_1D_02_catalog -c validation/validation_1D_commands -e 0 "$@"
bin/phosim.py validation/validation_1D_03_catalog -c validation/validation_1D_commands -e 0 "$@"
bin/phosim.py validation/validation_1D_04_catalog -c validation/validation_1D_commands -e 0 "$@"
bin/phosim.py validation/validation_1D_05_catalog -c validation/validation_1D_commands -e 0 "$@"
bin/phosim.py validation/validation_1D_06_catalog -c validation/validation_1D_commands -e 0 "$@"
bin/phosim.py validation/validation_1D_07_catalog -c validation/validation_1D_commands -e 0 "$@"
bin/phosim.py validation/validation_1D_08_catalog -c validation/validation_1D_commands -e 0 "$@"
bin/phosim.py validation/validation_1D_09_catalog -c validation/validation_1D_commands -e 0 "$@"

