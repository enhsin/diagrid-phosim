CC 	= ${CCS}
ERASE = rm -f
CFLAGS = -g -Wall -c  -I../
LFLAGS = -g -Wall   -lm
OBJECTS = instrument.o rng_mwc.o instrumentfiles.o phosim_parser.o readtext.o

all: instrument

instrument: $(OBJECTS)
	 $(CC) $(OBJECTS) $(LFLAGS) -o instrument

rng_mwc.o: ../raytrace/rng_mwc.cpp
	$(CC) $(CFLAGS) ../raytrace/rng_mwc.cpp -o rng_mwc.o

instrument.o:	instrument.cpp
	$(CC) $(CFLAGS) instrument.cpp -o instrument.o 

instrumentfiles.o:	instrumentfiles.cpp
	$(CC) $(CFLAGS) instrumentfiles.cpp -o instrumentfiles.o

phosim_parser.o: ../ancillary/phosim_parser.h ../ancillary/phosim_parser.cpp
	$(CC) $(CFLAGS) ../ancillary/phosim_parser.cpp -o phosim_parser.o

readtext.o: ../ancillary/readtext.h ../ancillary/readtext.cpp
	$(CC) $(CFLAGS) ../ancillary/readtext.cpp -o readtext.o

clean: 
	$(ERASE) $(OBJECTS)
