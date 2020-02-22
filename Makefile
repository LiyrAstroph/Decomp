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
DNEST_INCL  = -I /home/liyropt/Projects/GIT/CDNest/
DNEST_LIBS  = -L /home/liyropt/Projects/GIT/CDNest -ldnest

MPICHINCL     = $(shell pkg-config --cflags mpich) 
MPICHLIB    = $(shell pkg-config --libs mpich) 
endif

ifeq ($(SYSTEM), "Tianhe")
GSL_INCL =
GSL_LIBS = -lgsl -lgslcblas
MPICHLIB =
MPIINCL  =
LAPACK_INCL = -I /THL8/software/lapack/3.8.0-icc16-dyn/include
LAPACK_LIBS = #/THL8/software/lapack/3.7.0-icc16/liblapacke.a /THL8/software/lapack/3.7.0-icc16/liblapack.a /THL8/software/lapack/3.7.0-icc16/libblas.a 
FFTW_INCL =
FFTW_LIBS = -lfftw3

DNEST_INCL  = -I /THL8/home/liyanrong3/soft/CDNest
DNEST_LIBS  = -L /THL8/home/liyanrong3/soft/CDNest -ldnest
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
