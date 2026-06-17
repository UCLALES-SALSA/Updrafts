###############################################################
#
# Location of code (in $ROOT) and location where model is to be built $BIN
#
ROOT      :=..
BIN       = $(ROOT)/bin
ARCH      := $(shell uname)
#
# Source directory
SRC_GPF =$(ROOT)/src

OBJ_DIR  =$(ROOT)/obj
ECHO    = /bin/echo
RM      = /bin/rm -f


# Cygwin
# ======
INCL = -I/usr/local/include -I/usr/include
LIBS = '-L/lib -L/usr/local/lib -L/usr/lib -llapack -lblas -lnlopt -lm -lnetcdf -lnetcdff'
F90 = gfortran
FFLAGS = -O2 -fdefault-real-8
# Debug (slower but sometimes useful):
#FFLAGS = -O0 -fdefault-real-8 -fimplicit-none -fbacktrace -g -fbounds-check -ffpe-trap=invalid,zero,overflow,denormal -fcheck=all -Wall -Wno-tabs -Wextra

# Puhti - Intel
# =============
# - Note: no lnopt, so optimization does not work!
#module load intel-oneapi-compilers-classic/2021.6.0
#module load intel-oneapi-mpi/2021.6.0
#module load netcdf-fortran/4.5.4
#module load hdf5/1.12.2-mpi
#LIBS = '-L/appl/spack/v018/install-tree/intel-2021.6.0/hdf5-1.12.2-rzlc7b/lib -lhdf5_hl -lhdf5 -L/appl/spack/v018/install-tree/intel-2021.6.0/netcdf-fortran-4.5.4-nucfze/lib -lnetcdff -L/appl/spack/v018/install-tree/intel-2021.6.0/netlib-scalapack-2.2.0-7vqmei/lib -lscalapack'
#F90 = ifort
#FFLAGS =-O2 -real-size 64

n24 j25 j26 h26:
	cd $(OBJ_DIR); $(MAKE) LES_ARC=$@ \
	FFLAGS='$(FFLAGS) $(INCL)' F90=$(F90) \
	OUT=$(BIN)/gp.$@ LIBS=$(LIBS) \
	SRCGPF=$(SRC_GPF)


.PHONY : clean

#
# cleaning
# --------------------
#
clean:
	$(ECHO) "cleaning"
	$(RM) $(OBJ_DIR)/*.mod $(OBJ_DIR)/*.o $(BIN)/gp.*

FORCE: 
.PRECIOUS: $(LIBS)
