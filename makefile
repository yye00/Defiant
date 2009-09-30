include ${PETSC_DIR}/conf/base

all:		Defiant

DEFIANT_OBJS = Defiant.o CreateDestroy.o IMPES2Ph.o RelativePermeability.o Transmissibility.o \
               InterpolateProperties.o Coordinates.o PVT.o IMPES3Ph.o IO.o RockProperties.o \
               DefiantUtilities.o Examples.o NewtonRaphson2Ph.o NewtonRaphson3Ph.o

Defiant:  ${DEFIANT_OBJS} chkopts
	-${CLINKER} -g3 -O0 -o Defiant.exe ${DEFIANT_OBJS} ${PETSC_SNES_LIB}
	${RM} Defiant.o

	
