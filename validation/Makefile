CC 	= ${CCP}
ERASE = rm -f
CFLAGS = -g -Wall -c
LFLAGS = -g -Wall -lm -lfftw3
OBJECTS = unittest.o ../source/raytrace/raytrace.o

all: unittest

unittest: $(OBJECTS)
	 $(CC) $(OBJECTS) $(LFLAGS) -o unittest

unittest.o: unittest.cpp
	$(CC) $(CFLAGS) unittest.cpp -o unittest.o

clean: 
	$(ERASE) $(OBJECTS)
