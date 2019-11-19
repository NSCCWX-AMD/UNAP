FLAGS = -fPIC -DLABEL_INT -DSCALAR_DOUBLE -O3 -g -DSWTIMER
F_FLAGS = -cpp
DEBUG_FLAGS  = -DDEBUG
# DEBUG_FLAGS += -Wall -Wextra -Werror

## default
PLA = x86

ifeq (${PLA}, x86)
	## x86 compiler
	CPP       = icpc
	F90       = ifort
	AR        = ar
	FLINKER   = mpiifort -lstdc++
	LINKER    = mpiicc
	MPI_DIR   = /usr/sw-cluster/mpi2
	LIBS_NAME = -lunap
endif

ifeq (${PLA}, sw)
	## sw compiler
	CPP        = swg++453
	F90        = mpif90
	AR         = swar
	LINK_FLAGS = -static
	LINKER     = swld453 ${LINK_FLAGS}
	FLINKER    = swld453-fort ${LINK_FLAGS}
	MPI_DIR    = /usr/sw-mpp/mpi2
	LIBS_NAME  = -lswunap
	FLAGS     += -DSW_SLAVE
endif

MPI_INCS_DIR = -I${MPI_DIR}/include
MPI_LIBS_DIR = -L${MPI_DIR}/lib

ifeq (${PLA}, sw)
	CPP_FLAGS += -mieee
	F_FLAGS   += -OPT:IEEE_arithmetic=2
endif

TEST_ROOT = ${PWD}
LIBS_DIR = -L../build/lib
INCS_DIR = -I../build/lnInclude

ifeq (${PLA}, sw)
	UNAT_INCS_DIR = -I${TEST_ROOT}/../../UNAT/include
	LIBS_DIR  += -L${TEST_ROOT}/../../UNAT/lib
	LIBS_NAME += -lUNAT -lmetis

	#swlu
	INCS_DIR += -I/home/export/online1/swmore/opensource/swlu/include
	LIBS_DIR  += -L/home/export/online1/swmore/opensource/swlu/lib
	LIBS_NAME += -lswlu_mpi
	# LIBS_NAME += -lswlu
endif

EXE: EXE.o
	${LINKER} -o ${EXE} $^ ${LIBS_DIR} ${LIBS_NAME} ${MPI_LIBS_DIR} -lmpi -lmpicxx

EXE.o: ${EXE}.cpp
	${CPP} ${FLAGS} ${DEBUG_FLAGS} ${CPP_FLAGS} -c -o $@ ${INCS_DIR} ${UNAT_INCS_DIR} ${MPI_INCS_DIR} $^

FEXE: FEXE.o
	${FLINKER} -o ${FEXE} $^ ${LIBS_DIR} ${LIBS_NAME}

FEXE.o: ${FEXE}.f90
	${F90} ${FLAGS} ${F_FLAGS} -c -o $@ ${INCS_DIR} ${MPI_INCS_DIR} $^

clean:
	rm -f *.o
