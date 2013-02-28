EXE			+= main
SOURCES	+= main.cpp mesh.cpp residue.cpp jacobian.cpp
OBJECTS	+= main.o conditions.o mesh.o residue.o jacobian.o
#PETSC_ARCH=arch-linux2-cxx-debug
PETSC_ARCH =arch-linux2-cxx-opt

CFLAGS   = -Wall -Wextra
FFLAGS   =
CPPFLAGS = -I. -I${FEPIC_DIR} $(FEP_INCLUDE) -Wall -Wextra -fopenmp -m64 -msse2 -L$(FEP_LIBS_DIR) -lfepic $(FEP_LDFLAGS)
FPPFLAGS =

.PHONY: clean

# tests some variables
ifeq "" "$(wildcard ${FEPIC_DIR})"
$(error variable FEPIC_DIR was not defined or is an invalid directory)
endif 

ifeq "" "$(wildcard ${PETSC_DIR})"
$(error variable PETSC_DIR was not defined or is an invalid directory)
endif 

ifeq "" "$(wildcard ${PETSC_DIR}/${PETSC_ARCH})"
$(error variable PETSC_ARCH was not defined or is an invalid directory)
endif 

ifeq "" "${PETSC_ARCH}"
$(error PETSC_ARCH is was properly defined)
endif


#onde é procurado as bibliotecas
#SLINKER		+= -Wl,-rpath,/home/felipe/Slibs/smesh -L/home/felipe/Slibs/smesh

# "variables" must be included twice
include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules
include ${FEPIC_DIR}/conf/variables

PETSC_KSP_LIB += -L$(FEP_LIBS_DIR) -lfepic $(FEP_LDFLAGS)

main: ${OBJECTS} chkopts
	-${CLINKER} -o $(EXE) $(OBJECTS) ${PETSC_KSP_LIB}

conditions.o: conditions.cpp
	g++ -c $(FEP_INCLUDE) conditions.cpp -o conditions.o

clean::
	${RM} *.o


include ${PETSC_DIR}/conf/test


