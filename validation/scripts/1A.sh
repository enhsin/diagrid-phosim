#!/usr/bin/env bash

bin/phosim.py validation/validation_1A_00_catalog -c validation/validation_1A_00_commands -e 0 "$@" ; \
bin/phosim.py validation/validation_1A_02_catalog -c validation/validation_1A_02_commands -e 0 "$@" ; \
bin/phosim.py validation/validation_1A_04_catalog -c validation/validation_1A_04_commands -e 0 "$@" ; \
bin/phosim.py validation/validation_1A_00_catalog -c validation/validation_1A_10_commands -e 0 "$@" ; \
bin/phosim.py validation/validation_1A_02_catalog -c validation/validation_1A_12_commands -e 0 "$@" ; \
bin/phosim.py validation/validation_1A_04_catalog -c validation/validation_1A_14_commands -e 0 "$@" ; \
bin/phosim.py validation/validation_1A_00_catalog -c validation/validation_1A_20_commands -e 0 "$@" ; \
bin/phosim.py validation/validation_1A_02_catalog -c validation/validation_1A_22_commands -e 0 "$@" ; \
bin/phosim.py validation/validation_1A_04_catalog -c validation/validation_1A_24_commands -e 0 "$@" ; \
bin/phosim.py validation/validation_1A_01_catalog -c validation/validation_1A_01_commands -e 0 "$@" ; \
bin/phosim.py validation/validation_1A_03_catalog -c validation/validation_1A_03_commands -e 0 "$@" ; \
bin/phosim.py validation/validation_1A_01_catalog -c validation/validation_1A_11_commands -e 0 "$@" ; \
bin/phosim.py validation/validation_1A_03_catalog -c validation/validation_1A_13_commands -e 0 "$@" ; \
bin/phosim.py validation/validation_1A_01_catalog -c validation/validation_1A_21_commands -e 0 "$@" ; \
bin/phosim.py validation/validation_1A_03_catalog -c validation/validation_1A_23_commands -e 0 "$@" ; \
bin/phosim.py validation/validation_1A_05_catalog -c validation/validation_1A_25_commands -e 0 "$@" ; \
bin/phosim.py validation/validation_1A_05_catalog -c validation/validation_1A_05_commands -e 0 "$@"

