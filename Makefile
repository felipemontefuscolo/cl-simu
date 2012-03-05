EXE			+= main
SOURCES	+= main.cpp mesh.cpp residue.cpp jacobian.cpp
OBJECTS	+= main.o conditions.o mesh.o residue.o jacobian.o
#PETSC_ARCH=linux-gnu-c++-debug
PETSC_ARCH=linux-gnu-c++
#PETSC_DIR=/home/felipe/libs/petsc-3.2-p2
PETSC_DIR=/home/felipe/libs/petsc-dev

CFLAGS   = -Wall -Wextra
FFLAGS   = 
CPPFLAGS = -I. -I${FEPIC_DIR} $(FEP_INCLUDE) -Wall -Wextra -fopenmp -m64 -msse2 -L$(FEP_LIBS_DIR) -lfepic $(FEP_LDFLAGS)
FPPFLAGS = 

#onde é procurado as bibliotecas
#SLINKER		+= -Wl,-rpath,/home/felipe/Slibs/smesh -L/home/felipe/Slibs/smesh

#necessário incluir duas vezes
include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules
include ${FEPIC_DIR}/conf/variables

#PETSC_COMPILE  +=  -I. -I${FEPIC_DIR} $(FEP_INCLUDE) -Wall -Wextra -fopenmp -m64 -msse2 
PETSC_KSP_LIB += -L$(FEP_LIBS_DIR) -lfepic $(FEP_LDFLAGS)


#main: main.o conditions.o jacobian.o residue.o chkopts
main: ${OBJECTS} chkopts
	-${CLINKER} -o $(EXE) $(OBJECTS) ${PETSC_KSP_LIB}

conditions.o: conditions.cpp
	g++ -c -I${FEPIC_DIR} conditions.cpp -o conditions.o

#jacobian.o: jacobian.cpp
#	g++ -c -Wall -Wextra -I${FEPIC_DIR} jacobian.cpp -o jacobian.o

dclean:
	${RM} *.o

include ${PETSC_DIR}/conf/test
