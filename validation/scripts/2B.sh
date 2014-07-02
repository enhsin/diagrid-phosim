#!/usr/bin/env bash

bin/phosim.py validation/validation_2B_00_catalog -c validation/validation_2B_0_commands -e 0 "$@"
bin/phosim.py validation/validation_2B_01_catalog -c validation/validation_2B_0_commands -e 0 "$@"
bin/phosim.py validation/validation_2B_02_catalog -c validation/validation_2B_0_commands -e 0 "$@"
bin/phosim.py validation/validation_2B_03_catalog -c validation/validation_2B_0_commands -e 0 "$@"
bin/phosim.py validation/validation_2B_04_catalog -c validation/validation_2B_0_commands -e 0 "$@"
bin/phosim.py validation/validation_2B_05_catalog -c validation/validation_2B_0_commands -e 0 "$@"
bin/phosim.py validation/validation_2B_06_catalog -c validation/validation_2B_0_commands -e 0 "$@"
bin/phosim.py validation/validation_2B_07_catalog -c validation/validation_2B_0_commands -e 0 "$@"
bin/phosim.py validation/validation_2B_08_catalog -c validation/validation_2B_0_commands -e 0 "$@"
bin/phosim.py validation/validation_2B_09_catalog -c validation/validation_2B_0_commands -e 0 "$@"
bin/phosim.py validation/validation_2B_10_catalog -c validation/validation_2B_1_commands -e 0 "$@"
bin/phosim.py validation/validation_2B_11_catalog -c validation/validation_2B_1_commands -e 0 "$@"
bin/phosim.py validation/validation_2B_12_catalog -c validation/validation_2B_2_commands -e 0 "$@"
bin/phosim.py validation/validation_2B_13_catalog -c validation/validation_2B_3_commands -e 0 "$@"
bin/phosim.py validation/validation_2B_14_catalog -c validation/validation_2B_4_commands -e 0 "$@"

