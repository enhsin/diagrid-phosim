#!/usr/bin/env bash

bin/phosim.py validation/validation_2C_00_catalog -c validation/validation_2C_commands -e 0 "$@"
bin/phosim.py validation/validation_2C_01_catalog -c validation/validation_2C_commands -e 0 "$@"
bin/phosim.py validation/validation_2C_02_catalog -c validation/validation_2C_commands -e 0 "$@"
bin/phosim.py validation/validation_2C_03_catalog -c validation/validation_2C_commands -e 0 "$@"
bin/phosim.py validation/validation_2C_04_catalog -c validation/validation_2C_commands -e 0 "$@"
bin/phosim.py validation/validation_2C_05_catalog -c validation/validation_2C_commands -e 0 "$@"
bin/phosim.py validation/validation_2C_10_catalog -c validation/validation_2C_commands -e 0 "$@"
bin/phosim.py validation/validation_2C_11_catalog -c validation/validation_2C_commands -e 0 "$@"
bin/phosim.py validation/validation_2C_12_catalog -c validation/validation_2C_commands -e 0 "$@"
bin/phosim.py validation/validation_2C_13_catalog -c validation/validation_2C_commands -e 0 "$@"
bin/phosim.py validation/validation_2C_14_catalog -c validation/validation_2C_commands -e 0 "$@"
bin/phosim.py validation/validation_2C_15_catalog -c validation/validation_2C_commands -e 0 "$@"
bin/phosim.py validation/validation_2C_20_catalog -c validation/validation_2C_commands -e 0 "$@"
bin/phosim.py validation/validation_2C_21_catalog -c validation/validation_2C_commands -e 0 "$@"
bin/phosim.py validation/validation_2C_22_catalog -c validation/validation_2C_commands -e 0 "$@"
bin/phosim.py validation/validation_2C_23_catalog -c validation/validation_2C_commands -e 0 "$@"
bin/phosim.py validation/validation_2C_24_catalog -c validation/validation_2C_commands -e 0 "$@"
bin/phosim.py validation/validation_2C_25_catalog -c validation/validation_2C_commands -e 0 "$@"
rm validation/*_220?_*.gz
rm validation/*_221?_*.gz
rm validation/*_222?_*.gz

