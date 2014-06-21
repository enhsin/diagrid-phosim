CC 	= ${CCS}
ERASE = rm -f
CFLAGS = -g -Wall -c -I ../
LFLAGS = -g -Bstatic -lcfitsio -Wl,-Bdynamic -lfftw3 -lm -lpthread
OBJECT = atmosphere.o turb2d.o cloud.o airglow.o \
	rng_mwc.o fits_utils.o atmosphere_creator.o phosim_parser.o readtext.o 

all: atmosphere

rng_mwc.o: ../raytrace/rng_mwc.cpp
	$(CC) $(CFLAGS) ../raytrace/rng_mwc.cpp -o rng_mwc.o

phosim_parser.o: ../ancillary/phosim_parser.h ../ancillary/phosim_parser.cpp
	$(CC) $(CFLAGS) ../ancillary/phosim_parser.cpp -o phosim_parser.o

readtext.o: ../ancillary/readtext.h ../ancillary/readtext.cpp
	$(CC) $(CFLAGS) ../ancillary/readtext.cpp -o readtext.o

fits_utils.o: ../ancillary/fits_utils.cpp ../ancillary/fits_utils.h
	$(CC) $(CFLAGS) ../ancillary/fits_utils.cpp -o fits_utils.o

atmosphere: $(OBJECT)
	$(CC) $(OBJECT) $(LFLAGS) -o atmosphere

atmosphere.o: atmosphere.cpp airglow.h cloud.h turb2d.h
	$(CC) $(CFLAGS) atmosphere.cpp -o atmosphere.o

atmosphere_creator.o: atmosphere_creator.h atmosphere_creator.cpp
	$(CC) $(CFLAGS) atmosphere_creator.cpp -o atmosphere_creator.o

turb2d.o: turb2d.cpp
	$(CC) $(CFLAGS) turb2d.cpp -o turb2d.o

cloud.o: cloud.cpp cloud.h ../ancillary/fits_utils.h ../ancillary/fftw_utils.h
	$(CC) $(CFLAGS) cloud.cpp -o cloud.o

airglow.o: airglow.cpp airglow.h ../ancillary/fits_utils.h ../ancillary/fftw_utils.h
	$(CC) $(CFLAGS) airglow.cpp -o airglow.o

test_atmosphere.o: test_atmosphere.cpp ../ancillary/phosim_parser.h
	$(CC) $(CFLAGS) test_atmosphere.cpp -o test_atmosphere.o

test_atmosphere: test_atmosphere.o phosim_parser.o
	$(CC) test_atmosphere.o phosim_parser.o $(LFLAGS) -o test_atmosphere

clean: 
	$(ERASE) $(OBJECT) atmosphere test_atmosphere
