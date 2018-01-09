VERSION = 2014
TARGET = Linux-i686
#
PROG = mga-$(VERSION).x
# Choose system 
# Linux: 
#        -malign-double uses 64 bit data alignment. This flag does not 
#                       work in Opteron machines, since this machines 
#                       are already 64 bits.
#   FFLAGS = -O3 -i4 -ident -Wall -static -malign-double
   FFLAGS = 
   CC  = gcc -lm
   FC  = gfortran
   FFC = gfortran
#
START    = $(shell date)
WORKDIR  = $(shell pwd)
DATE     = $(shell date +%Y-%m-%d)
#
# libraries
#
BLAS_LIB = -L/usr/lib/sse -lblas -L/usr/lib/sse/atlas -lblas
#
#
FILES = makefile input.example
#
#
MAIN_OBJ = src/main.o 
#
MGA_OBJS = src/rdinput.o src/alloc.o src/exitg.o src/readtable.o src/genindv.o \
	   src/centerxyz.o src/rotate.o src/placeatomonzaxis.o \
	   src/placeatomonplan.o src/ran0.o src/geninp.o src/genxyz.o \
	   src/runprog.o src/getxyz.o src/order.o src/placeatomonyzplan.o\
	   src/getenergy.o src/mating.o src/ran1.o src/gasdev.o src/expdev.o \
           src/getxyzga.o src/rdzmat.o src/xyz2zmat.o src/genzmat.o \
	   src/zmat2xyz.o src/matsc.o src/matdc.o src/mutation.o \
           src/rtgroup.o src/predator.o src/reoptpop.o src/distmeasure.o \
           src/nga.o src/iga.o src/confirmconv.o src/history.o src/orderfpop.o \
	   src/getxyzobabel.o  src/getenergyobabel.o src/popbackup.o src/readbackup.o \
# 
OBJS = $(MAIN_OBJ) $(MGA_OBJS) 
#
mga:    $(OBJS)
	$(CC) $(FFLAGS) -o $(PROG) $(OBJS) 
	size $(PROG)
	ln -sf $(PROG) mga.x
	chmod a+xr $(PROG) mga.x
#
#
clean:	 
	rm -f $(OBJS)
#
veryclean:  
	make clean
	rm -f $(PROG) src/*.o
#
distclean:
	make veryclean 
	rm -rf  *.* *~ core "*#" "#*" "*#" "#*" tes* \
		src/cab.f src/*~ src/*.log src/*.dat src/*.ptc src/*.spc \
		src/*.xmgr src/*.x src/tes* src/LIXO src/Lixo src/lixo \
		src/core src/a.out doc/*~ doc/core src/"*#" src/"#*" \
		doc/"*#" doc/"#*" 
#