#!/usr/bin/env bash

bin/phosim.py validation/validation_1C_00_catalog -c validation/validation_1C_commands -e 0 "$@" 
bin/phosim.py validation/validation_1C_01_catalog -c validation/validation_1C_commands -e 0 "$@" 
bin/phosim.py validation/validation_1C_02_catalog -c validation/validation_1C_commands -e 0 "$@" 
bin/phosim.py validation/validation_1C_03_catalog -c validation/validation_1C_commands -e 0 "$@" 
bin/phosim.py validation/validation_1C_04_catalog -c validation/validation_1C_commands -e 0 "$@" 
bin/phosim.py validation/validation_1E_00_catalog -c validation/validation_1E_commands -e 0 "$@" 
bin/phosim.py validation/validation_1E_01_catalog -c validation/validation_1E_commands -e 0 "$@" 
bin/phosim.py validation/validation_1E_02_catalog -c validation/validation_1E_commands -e 0 "$@" 
bin/phosim.py validation/validation_1E_03_catalog -c validation/validation_1E_commands -e 0 "$@" 
bin/phosim.py validation/validation_1E_04_catalog -c validation/validation_1E_commands -e 0 "$@" 
bin/phosim.py validation/validation_1F_00_catalog -c validation/validation_1F_commands -e 0 "$@" 
bin/phosim.py validation/validation_1F_01_catalog -c validation/validation_1F_commands -e 0 "$@" 
bin/phosim.py validation/validation_1F_02_catalog -c validation/validation_1F_commands -e 0 "$@"
bin/phosim.py validation/validation_1F_03_catalog -c validation/validation_1F_commands -e 0 "$@" 
bin/phosim.py validation/validation_1F_04_catalog -c validation/validation_1F_commands -e 0 "$@" 
date > validation/speed_out_0 
bin/phosim.py validation/validation_4A_00_catalog -c validation/validation_4A_00_commands "$@" 
date > validation/speed_out_1
bin/phosim.py validation/validation_4A_01_catalog -c validation/validation_4A_01_commands "$@" >& validation/speed_out
date > validation/speed_out_2 

