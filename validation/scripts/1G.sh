#!/usr/bin/env bash

bin/phosim.py validation/validation_1G_01_catalog -c validation/validation_1G_commands "$@" -e 0
bin/phosim.py validation/validation_1G_02_catalog -c validation/validation_1G_commands "$@" -e 0
bin/phosim.py validation/validation_1G_03_catalog -c validation/validation_1G_commands "$@" -e 0
bin/phosim.py validation/validation_1G_04_catalog -c validation/validation_1G_commands "$@" -e 0
bin/phosim.py validation/validation_1G_11_catalog -c validation/validation_1G_commands "$@" -e 0
bin/phosim.py validation/validation_1G_12_catalog -c validation/validation_1G_commands "$@" -e 0
bin/phosim.py validation/validation_1G_13_catalog -c validation/validation_1G_commands "$@" -e 0
bin/phosim.py validation/validation_1G_14_catalog -c validation/validation_1G_commands "$@" -e 0


