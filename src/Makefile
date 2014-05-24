all:
	make -C source/raytrace CCP=g++ CCS=g++
	make -C source/trim CCP=g++ CCS=g++
	make -C source/e2adc CCP=g++ CCS=g++
	make -C source/instrument CCP=g++ CCS=g++
	make -C source/atmosphere CCP=g++ CCS=g++
	make -C source/ancillary CCP=g++ CCS=g++
	make -C validation CCP=g++ CCS=g++

install: all
	install -D source/raytrace/raytrace ../bin/
	install -D source/trim/trim ../bin/
	install -D source/e2adc/e2adc ../bin/
	install -D source/instrument/instrument ../bin/
	install -D source/atmosphere/atmosphere ../bin/
	install -D validation/unittest ../bin/

clean:
	rm -f source/raytrace/*.o; rm -f source/trim/*.o; rm -f source/e2adc/*.o; rm -f source/instrument/*.o; rm -f source/atmosphere/*.o ; rm -f source/ancillary/*.o; rm -f validation/*.o

distclean: clean
	rm -f ../bin/raytrace; rm -f ../bin/trim; rm -f ../bin/e2adc; rm -f ../bin/instrument; rm -f ../bin/atmosphere; rm -f  ../bin/unittest

