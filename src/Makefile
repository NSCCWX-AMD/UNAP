FLAGS = -fPIC -DLABEL_INT -DSCALAR_DOUBLE -g -O3 -DSWTIMER
DEBUG_FLAGS  = -DDEBUG
# DEBUG_FLAGS += -Wall -Wextra -Werror

## default
PLA = x86

ifeq (${PLA}, x86)
	## x86 compiler
	CPP       = icpc
	C  		  = icc
	F90       = ifort
	AR        = ar
	FLINKER   = mpiifort -lstdc++
	LINKER    = mpiicc
	MPI_DIR   = /usr/sw-cluster/mpi2
	LIB_NAME  = unap
	LIB_NAME_FULL = lib${LIB_NAME}.so
endif

ifeq (${PLA}, sw)
	## sw compiler
	CPP           = swg++453
	C             = sw5cc
	F90           = mpif90
	AR            = swar
	LINK_FLAGS    = -static
	LINKER        = swld453 ${LINK_FLAGS}
	FLINKER       = swld453-fort ${LINK_FLAGS}
	MPI_DIR       = /usr/sw-mpp/mpi2
	LIB_NAME      = swunap
	LIB_NAME_FULL = lib${LIB_NAME}.a
endif

MPI_INCS_DIR = ${MPI_DIR}/include

ifeq (${PLA}, sw)
	CPP_FLAGS += -mieee
	F_FLAGS   += -OPT:IEEE_arithmetic=2
endif

SRC_ROOT = ${PWD}

VPATH = ${SRC_ROOT}/base/:\
		${SRC_ROOT}/base/vectors/:\
		${SRC_ROOT}/base/mpi/:\
		${SRC_ROOT}/matrix/:\
		${SRC_ROOT}/matrix/lduMatrix/:\
		${SRC_ROOT}/matrix/lduMatrix/lduAgglomeration/:\
		${SRC_ROOT}/matrix/csrMatrix/:\
		${SRC_ROOT}/matrix/interfaces/:\
		${SRC_ROOT}/matrix/interfaces/labelPair/:\
		${SRC_ROOT}/matrix/interfaces/patch/:\
		${SRC_ROOT}/matrix/eigenvalue/eigenDiagPCG/:\
		${SRC_ROOT}/solvers/PCG/:\
		${SRC_ROOT}/solvers/MG/:\
		${SRC_ROOT}/solvers/PBiCGStab/:\
		${SRC_ROOT}/preconditioners/lduPreconditioners/lduDiagonal/:\
		${SRC_ROOT}/preconditioners/lduPreconditioners/lduDIC/:\
		${SRC_ROOT}/preconditioners/lduPreconditioners/lduDILU/:\
		${SRC_ROOT}/smoothers/lduSmoothers/lduGaussSeidelSmoother/:\
		${SRC_ROOT}/smoothers/chebySmoother/:\
		${SRC_ROOT}/tools/timers/:\
		${SRC_ROOT}/tools/fortranInterface/:\
		${SRC_ROOT}/tools/matrixConversion/:\
		${SRC_ROOT}/tools/readData/readFromOpenFOAM/:\
		${SRC_ROOT}/tools/readData/readFromHypre/:\
		${SRC_ROOT}/tools/printData/printUNAP/:\
		${SRC_ROOT}/tools/swArrays: ${SRC_ROOT}/tools/swArrays/restrictInterpolation:\
		${SRC_ROOT}/tools/swAthread: \
		${SRC_ROOT}/tools/reorderMatrix: \
		${SRC_ROOT}/tools/rcm

CPPFILES = unap.cpp unapVectorOperations.cpp unapMPI.cpp \
		   unapMatrixSolver.cpp unapMatrixTest.cpp unapMatrixPrecond.cpp unapMatrixAgglomeration.cpp\
		   PCG.cpp PBiCGStab.cpp \
		   MG.cpp scalingFactor.cpp MGSolve.cpp Vcycle.cpp solveCoarsestLevel.cpp \
		   interfaces.cpp labelPair.cpp patch.cpp \
		   eigenDiagPCG.cpp \
		   lduMatrix.cpp lduMatrixSpMV.cpp lduAgglomeration.cpp lduAgglomerate.cpp \
		   lduAgglomerationAddressing.cpp lduAgglomerateMatrix.cpp lduAgglomerationCombineLevels.cpp \
		   lduDiagPrecond.cpp lduDICPrecond.cpp lduDILUPrecond.cpp \
		   lduGaussSeidelSmoother.cpp  chebySmoother.cpp \
		   swTimer.cpp \
		   fortranInterface.cpp \
		   matrixConversion.cpp coo2ldu.cpp \
		   readFromOpenFOAM.cpp readFromHypre.cpp \
		   printUNAP.cpp cooRCM.cpp

CFILES = rcm.c



ROOT_DIR = ${SRC_ROOT}/..
OBJS_DIR = ${ROOT_DIR}/build/objs
INCS_DIR = ${ROOT_DIR}/build/lnInclude
LIBS_DIR = ${ROOT_DIR}/build/lib

DIRS_CREATED = ${OBJS_DIR} ${INCS_DIR} ${LIBS_DIR}
${shell mkdir -p ${DIRS_CREATED}}


CPPOBJS  = ${CPPFILES:%.cpp=${OBJS_DIR}/%.o}
COBJS    = ${CFILES:%.c=${OBJS_DIR}/%.o}
OBJS_ALL = ${CPPOBJS}
OBJS_ALL+= ${COBJS}

.PHONY: all link clean

ifeq (${PLA}, sw)
CHOST_FLAGS   = -host -OPT:IEEE_arith=2
CSLAVE_FLAGS  = -slave -OPT:IEEE_arith=2 -msimd
CPP_FLAGS    += -DSW_SLAVE
UNAT_INCS_DIR = -I${ROOT_DIR}/../UNAT/include

CHOSTFILES  = gSum_host.c residualNormFactor_host.c vectorOps_host.c \
			  vectorCopy_host.c scalingFactor_host.c swAthread.c \
			  interpolateData_host.c restrictData_host.c

CSLAVEFILES = gSum_slave.c residualNormFactor_slave.c vectorOps_slave.c \
			  vectorCopy_slave.c userFunc_SAXPY.c scalingFactor_slave.c \
			  interpolateData_slave.c restrictData_slave.c

CPPFILES += lduAgglomerationReorder.cpp swAgglomeration.cpp

CHOSTOBJS  = ${CHOSTFILES:%.c=${OBJS_DIR}/%.o}
CSLAVEOBJS = ${CSLAVEFILES:%.c=${OBJS_DIR}/%.o}

OBJS_ALL += ${CHOSTOBJS} ${CSLAVEOBJS}

${CHOSTOBJS}: ${OBJS_DIR}/%.o : %.c
	${C} ${CHOST_FLAGS} -c -o $@ -I${INCS_DIR} $^

${CSLAVEOBJS}: ${OBJS_DIR}/%.o : %.c
	${C} ${CSLAVE_FLAGS} -c -o $@ -I${INCS_DIR} $^
endif

all: ${OBJS_ALL}
	${AR} cur ${LIBS_DIR}/${LIB_NAME_FULL} $^

link:
	rm -f ${INCS_DIR}/*
	find ${SRC_ROOT}/ -path ${INCS_DIR} -prune -o -name "*.h*" -type f -exec ln -s {} ${INCS_DIR}/ \;

${CPPOBJS}: ${OBJS_DIR}/%.o : %.cpp
	${CPP} ${FLAGS} ${DEBUG_FLAGS} ${CPP_FLAGS} -c -o $@ -I${INCS_DIR} ${UNAT_INCS_DIR} -I${MPI_INCS_DIR} $^

${COBJS}: ${OBJS_DIR}/%.o : %.c
	${C} ${CHOST_FLAGS} -c -o $@ -I${INCS_DIR} $^

clean:
	rm -rf ${SRC_ROOT}/../build/objs
