CC 	= ${CCP}
ERASE = rm -f
CFLAGS = -g -Wall -c -I ../
LFLAGS = -g -Wall -lm -lz
OBJECTS = trim.o readtext.o

all: trim

trim: $(OBJECTS)
	 $(CC) $(OBJECTS) $(LFLAGS) -o ./trim
	cp trim ../../bin/.

readtext.o: ../ancillary/readtext.h ../ancillary/readtext.cpp
	$(CC) $(CFLAGS) ../ancillary/readtext.cpp -o readtext.o

trim.o: trim.cpp
	$(CC) $(CFLAGS) trim.cpp -o trim.o

clean:
	$(ERASE) $(OBJECTS)
