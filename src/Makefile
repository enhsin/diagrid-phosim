

all:
	rm -f ../source/raytrace/*.o; rm -f ../source/trim/*.o; rm -f ../source/e2adc/*.o; rm -f ../source/instrument/*.o; rm -f ../source/atmosphere/*.o ; rm -f ../source/ancillary/*.o
	csh -c "source ../bin/setup ; cd ../source/raytrace ; make"
	csh -c "source ../bin/setup ; cd ../source/trim ; make"
	csh -c "source ../bin/setup ; cd ../source/e2adc ; make"
	csh -c "source ../bin/setup ; cd ../source/instrument ; make"
	csh -c "source ../bin/setup ; cd ../source/atmosphere ; make"
	csh -c "source ../bin/setup ; cd ../source/ancillary ; make "
	csh -c "source ../bin/setup ; cd ../validation ; make "

install:
	install -D ../source/raytrace/raytrace ../bin/
	install -D ../source/trim/trim ../bin/
	install -D ../source/e2adc/e2adc ../bin/
	install -D ../source/instrument/instrument ../bin/
	install -D ../source/atmosphere/atmosphere ../bin/

clean:
	rm -f ../source/raytrace/*.o; rm -f ../source/trim/*.o; rm -f ../source/e2adc/*.o; rm -f ../source/instrument/*.o; rm -f ../source/atmosphere/*.o ; rm -f ../source/ancillary/*.o

distclean: clean
	rm -f ../bin/raytrace; rm -f ../bin/trim; rm -f ../bin/e2adc; rm -f ../bin/instrument; rm -f ../bin/atmosphere
