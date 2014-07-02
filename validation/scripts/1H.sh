#!/usr/bin/env bash

bin/phosim.py validation/validation_1H_01_catalog -c validation/validation_1H_commands "$@" -e 0 --keepscreens=1
bin/phosim.py validation/validation_1H_02_catalog -c validation/validation_1H_commands "$@" -e 0 --keepscreens=1
bin/phosim.py validation/validation_1H_03_catalog -c validation/validation_1H_commands "$@" -e 0 --keepscreens=1
bin/phosim.py validation/validation_1H_04_catalog -c validation/validation_1H_commands "$@" -e 0 --keepscreens=1
bin/phosim.py validation/validation_1H_05_catalog -c validation/validation_1H_commands "$@" -e 0 --keepscreens=1
bin/phosim.py validation/validation_1H_06_catalog -c validation/validation_1H_commands "$@" -e 0 --keepscreens=1
bin/phosim.py validation/validation_1H_07_catalog -c validation/validation_1H_commands "$@" -e 0 --keepscreens=1
bin/phosim.py validation/validation_1H_08_catalog -c validation/validation_1H_commands "$@" -e 0 --keepscreens=1
bin/phosim.py validation/validation_1H_09_catalog -c validation/validation_1H_commands "$@" -e 0 --keepscreens=1
rm validation/lsst_e_199?_*


