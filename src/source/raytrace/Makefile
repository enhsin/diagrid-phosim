CC = ${CCP}
ERASE = rm -f
CFLAGS = -g -O3 -ffast-math -Wall -c -I ../
LFLAGS = -g -O3 -ffast-math -Wl,-Bstatic -lcfitsio -Wl,-Bdynamic -lfftw3 -lm -lz -lpthread
OBJECTS =   main.o observation.o image.o rng_mwc.o raytrace.o galaxy.o dust.o event.o grating.o counter.o silicon.o readtext.o surface.o coating.o perturbation.o obstruction.o screen.o medium.o air.o chip.o contamination.o fea.o

all: raytrace

raytrace: $(OBJECTS)
	 $(CC) $(OBJECTS) $(LFLAGS) -o ./raytrace

main.o: main.cpp
	$(CC) $(CFLAGS) main.cpp -o main.o

observation.o: observation.cpp
	$(CC) $(CFLAGS) observation.cpp -o observation.o

image.o: image.cpp
	$(CC) $(CFLAGS) image.cpp -o image.o

rng_mwc.o: rng_mwc.cpp
	$(CC) $(CFLAGS) rng_mwc.cpp -o rng_mwc.o

raytrace.o: raytrace.cpp
	$(CC) $(CFLAGS) raytrace.cpp -o raytrace.o

galaxy.o: galaxy.cpp
	$(CC) $(CFLAGS) galaxy.cpp -o galaxy.o

dust.o: dust.cpp
	$(CC) $(CFLAGS) dust.cpp -o dust.o

event.o: event.cpp
	$(CC) $(CFLAGS) event.cpp -o event.o

silicon.o: silicon.cpp
	$(CC) $(CFLAGS) silicon.cpp -o silicon.o

surface.o: surface.cpp
	$(CC) $(CFLAGS) surface.cpp -o surface.o

coating.o: coating.cpp
	$(CC) $(CFLAGS) coating.cpp -o coating.o

perturbation.o: perturbation.cpp
	$(CC) $(CFLAGS) perturbation.cpp -o perturbation.o

screen.o: screen.cpp
	$(CC) $(CFLAGS) screen.cpp -o screen.o

medium.o: medium.cpp
	$(CC) $(CFLAGS) medium.cpp -o medium.o

air.o: air.cpp
	$(CC) $(CFLAGS) air.cpp -o air.o

obstruction.o: obstruction.cpp
	$(CC) $(CFLAGS) obstruction.cpp -o obstruction.o

contamination.o: contamination.cpp
	$(CC) $(CFLAGS) contamination.cpp -o contamination.o

chip.o: chip.cpp
	$(CC) $(CFLAGS) chip.cpp -o chip.o

grating.o: grating.cpp
	$(CC) $(CFLAGS) grating.cpp -o grating.o

counter.o: counter.cpp
	$(CC) $(CFLAGS) counter.cpp -o counter.o

readtext.o: ../ancillary/readtext.h ../ancillary/readtext.cpp
	$(CC) $(CFLAGS) ../ancillary/readtext.cpp -o readtext.o

fea.o: fea.cpp
	$(CC) $(CFLAGS) fea.cpp -o fea.o

clean:
	$(ERASE) $(OBJECTS)
