#!/usr/bin/env bash

bin/phosim.py validation/validation_3A_0_catalog -c validation/validation_3A_commands "$@" 
bin/phosim.py validation/validation_3A_1_catalog -c validation/validation_3A_commands "$@" 
bin/phosim.py validation/validation_3A_2_catalog -c validation/validation_3A_commands "$@" 
bin/phosim.py validation/validation_3A_3_catalog -c validation/validation_3A_commands "$@" 
bin/phosim.py validation/validation_1A_05_catalog -c validation/validation_1A_15_commands -e 0 "$@"

