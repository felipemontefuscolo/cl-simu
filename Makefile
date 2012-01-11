EXE			+= main
SOURCES	+= main.cpp
OBJECTS	+= main.o conditions.o
#PETSC_ARCH=linux-gnu-c++-debug
PETSC_ARCH=linux-gnu-c++

CFLAGS   = -Wall -Wextra
FFLAGS   = 
CPPFLAGS = -I. -I${FEPIC_DIR} $(FEP_INCLUDE) -Wall -Wextra -fopenmp -m64 -msse2 -L$(FEP_LIBS_DIR) -lfepic $(FEP_LDFLAGS)
FPPFLAGS = 

#onde é procurado as bibliotecas
#SLINKER		+= -Wl,-rpath,/home/felipe/Slibs/smesh -L/home/felipe/Slibs/smesh

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules
include ${FEPIC_DIR}/conf/variables

#PETSC_COMPILE  +=  -I. -I${FEPIC_DIR} $(FEP_INCLUDE) -Wall -Wextra -fopenmp -m64 -msse2 
PETSC_KSP_LIB += -L$(FEP_LIBS_DIR) -lfepic $(FEP_LDFLAGS)


main: main.o conditions.o chkopts
	-${CLINKER} -o $(EXE) $(OBJECTS) ${PETSC_KSP_LIB}


conditions.o: conditions.cpp
	g++ -c -Wall -Wextra -I${FEPIC_DIR} conditions.cpp -o conditions.o

dclean:
	${RM} *.o

include ${PETSC_DIR}/conf/test
