CC = ${CCS}
ERASE = rm -f
CFLAGS = -g -Wall -c -I ../ -I ../../
LFLAGS = -g -Wall -lfftw3 -lcfitsio -lm -lz
UT_CFLAGS = -I${BOOST_INC_DIR} 
UT_LFLAGS = -L${BOOST_LIB_DIR} -lboost_unit_test_framework

OBJECTS = test_phosim_parser.o test_fits_utils.o test_file_reader.o \
	test_ancillary.o phosim_parser.o readtext.o rng_mwc.o \
	fits_utils.o file_reader.o 

all:
	make test_phosim_parser
	make test_fits_utils
	make test_readtext
#	make test_ancillary

rng_mwc.o: ../raytrace/rng_mwc.cpp
	$(CC) $(CFLAGS) ../raytrace/rng_mwc.cpp -o rng_mwc.o

test_phosim_parser: test_phosim_parser.o phosim_parser.o
	$(CC) test_phosim_parser.o phosim_parser.o $(LFLAGS) \
	-o test_phosim_parser

test_readtext: test_readtext.o readtext.o
	$(CC) test_readtext.o readtext.o $(LFLAGS) \
	-o test_readtext

test_fits_utils: test_fits_utils.o fits_utils.o rng_mwc.o
	$(CC) test_fits_utils.o fits_utils.o rng_mwc.o $(LFLAGS) \
	-o test_fits_utils

test_file_reader: test_file_reader.o file_reader.h file_reader.o
	$(CC) test_file_reader.o file_reader.o $(LFLAGS) $(UT_LFLAGS) \
	-o test_file_reader

test_ancillary: test_ancillary.o phosim_parser.o fits_utils.o rng_mwc.o
	$(CC) test_ancillary.o phosim_parser.o fits_utils.o rng_mwc.o \
	$(LFLAGS) $(UT_LFLAGS) \
	-o test_ancillary

phosim_parser.o: phosim_parser.h phosim_parser.cpp
	$(CC) $(CFLAGS) phosim_parser.cpp -o phosim_parser.o

readtext.o: readtext.h readtext.cpp
	$(CC) $(CFLAGS) readtext.cpp -o readtext.o

fits_utils.o: fits_utils.h fits_utils.cpp
	$(CC) $(CFLAGS) fits_utils.cpp -o fits_utils.o

file_reader.o: file_reader.h file_reader.cpp
	$(CC) $(CFLAGS) file_reader.cpp -o file_reader.o

test_phosim_parser.o: test_phosim_parser.cpp phosim_parser.h
	$(CC) $(CFLAGS) test_phosim_parser.cpp -o test_phosim_parser.o

test_readtext.o: test_readtext.cpp readtext.h
	$(CC) $(CFLAGS) test_readtext.cpp -o test_readtext.o

test_fits_utils.o: test_fits_utils.cpp fits_utils.h
	$(CC) $(CFLAGS) test_fits_utils.cpp -o test_fits_utils.o

test_file_reader.o: test_file_reader.cpp file_reader.h
	$(CC) $(CFLAGS) $(UT_CFLAGS) test_file_reader.cpp -o test_file_reader.o

test_ancillary.o: test_ancillary.cpp fits_utils.h phosim_parser.h
	$(CC) $(CFLAGS) $(UT_CFLAGS) test_ancillary.cpp \
	-o test_ancillary.o
clean: 
	$(ERASE) $(OBJECTS) test_phosim_parser test_fits_utils \
	test_file_reader test_ancillary test_readtext
