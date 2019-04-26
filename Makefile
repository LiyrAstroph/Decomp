HELL=/bin/bash

CC       = mpicc
OPTIMIZE = -O2 -Wall -finline-functions
#OPTIMIZE += -DDebug

#------------target system---------
SYSTEM="Linux"

ifeq ($(SYSTEM), "Linux")
NCORE      :=$(grep -c ^processor /proc/cpuinfo)
GSL_INCL    = $(shell pkg-config --cflags gsl) 
GSL_LIBS    = $(shell pkg-config --libs gsl) 
LAPACK_INCL = -I/usr/include/lapacke
LAPACK_LIBS = -L/usr/lib64 -llapacke -llapack -lblas
DNEST_INCL  = -I /home/liyropt/Projects/GIT/DNest/
DNEST_LIBS  = -L /home/liyropt/Projects/GIT/DNest -ldnest

MPICHINCL     = $(shell pkg-config --cflags mpich) 
MPICHLIB    = $(shell pkg-config --libs mpich) 
endif


EXEC     = decomp
SRC      = ./src
OBJS     = $(SRC)/main.o $(SRC)/allvars.o $(SRC)/mathfun.o $(SRC)/system.o \
           $(SRC)/read.o $(SRC)/run.o $(SRC)/help.o $(SRC)/recon.o         \
					 $(SRC)/mc_dnest.o
 
INCL     = Makefile  $(SRC)/allvars.h $(SRC)/proto.h $(SRC)/mc_dnest.h

OPTIONS  = $(OPTIMIZE)
CFLAGS   = $(OPTIONS) $(GSL_INCL) $(LAPACK_INCL) $(MPICHINCL) $(DNEST_INCL)
LIBS     = $(GSL_LIBS) $(LAPACK_LIBS) $(MPICHLIB) $(DNEST_LIBS)

$(EXEC):$(OBJS)
	cd $(SRC)
	$(CC) $(OBJS) $(LIBS) -o $@

$(OBJS): $(INCL)

clean:
	rm -f $(OBJS) $(EXEC)