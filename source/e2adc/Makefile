CC 	= ${CCP}
ERASE = rm -f
CFLAGS = -g -O3 -Wall -c -I${CFIO_INC_DIR} -I ../
LFLAGS = -g -Wall -I${CFIO_INC_DIR} -L${CFIODIR} -lcfitsio -lm
OBJECTS = e2adc.o rng_mwc.o readtext.o fits_utils.o

all: e2adc

e2adc: $(OBJECTS)
	$(CC) $(OBJECTS) $(LFLAGS) -o e2adc
	cp e2adc ../../bin            

rng_mwc.o: ../raytrace/rng_mwc.cpp
	$(CC) $(CFLAGS) ../raytrace/rng_mwc.cpp -o rng_mwc.o

readtext.o: ../ancillary/readtext.h ../ancillary/readtext.cpp
	$(CC) $(CFLAGS) ../ancillary/readtext.cpp -o readtext.o

fits_utils.o: ../ancillary/fits_utils.cpp ../ancillary/fits_utils.h
	$(CC) $(CFLAGS) ../ancillary/fits_utils.cpp -o fits_utils.o

e2adc.o: e2adc.cpp ../ancillary/fits_utils.h
	$(CC) $(CFLAGS) e2adc.cpp -o e2adc.o

clean: 
	$(ERASE) $(OBJECTS)
