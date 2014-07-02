#!/usr/bin/env bash

bin/phosim.py validation/validation_2A_0_catalog -c validation/validation_2A_commands -e 0 "$@"
mv validation/output.fits validation/output_0.fits
bin/phosim.py validation/validation_2A_1_catalog -c validation/validation_2A_commands -e 0 "$@"
mv validation/output.fits validation/output_1.fits
bin/phosim.py validation/validation_2A_2_catalog -c validation/validation_2A_commands -e 0 "$@"
mv validation/output.fits validation/output_2.fits
bin/phosim.py validation/validation_2A_3_catalog -c validation/validation_2A_commands -e 0 "$@"
mv validation/output.fits validation/output_3.fits
bin/phosim.py validation/validation_2A_4_catalog -c validation/validation_2A_commands -e 0 "$@"
mv validation/output.fits validation/output_4.fits

