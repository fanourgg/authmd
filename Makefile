EXE=authmd
EXEDIR=../exe/

$(shell makedepf90 ${DEFS} -o ${EXEDIR}${EXE} *.f90  > Makefile_deps)
include Makefile_deps

all: ${EXE}
$(FOBJ): Makefile

#FC=gfortran
#FFLAGS=
FC=ifort
FFLAGS=-fpp
#FFLAGS=-C -debug -traceback -check noarg_temp_created -warn unused -fpp -openmp
#FFLAGS=-C -debug -traceback -fpp #-openmp #-warn unused
FFLAGS=-C -check all -debug -traceback -fpp -warn  notruncated_source -warn unused
#FFLAGS=-fast -ftz -fpp -openmp
FFLAGS=-O2 -ftz -fpp -Bstatic #-openmp
#FFLAGS=-check all -warn all,nodec,interfaces -gen_interfaces -traceback -fpe0 -fpstkchk -fpp
LDFLAGS=${FFLAGS}
#FC=gfortran
#FFLAGS=-O -cpp

DEFS=
INC=
LIBS=

deps:
	makedepf90 ${DEFS} -o ${EXE} *.f90 > Makefile_deps

%.o:%.f90
	$(FC) $(FFLAGS) $(DEFS) $(INC) -c $< -o $@

clean: 
	rm -rf $(FOBJ) *.o *.mod Makefile_deps ${EXE}


