CC = m68k-palmos-gcc
CFLAGS = -g -mdebug-labels -O2 -Wall

all: hchem.prc

hchem.prc: hchem bin.stamp
	build-prc hchem.prc "HandyChem" HCHM hchem *.bin

bin.stamp: hchem.rcp resids.h gnu.pbitm
	pilrc hchem.rcp
	touch bin.stamp

hchem: hchem.o u_struct.o
hchem.o: hchem.c
u_struct.o: u_struct.c

clean:
	-rm -f *.o hchem *.bin *.stamp *.[pg]rc
