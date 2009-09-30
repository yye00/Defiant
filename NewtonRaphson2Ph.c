/*
 * NewtonRaphson2Ph.c
 *
 *  Created on: Sep 20, 2009
 *      Author: yye00
 */

#include "Defiant.h"

#undef __FUNCT__
#define __FUNCT__ "DefiantNewtonRaphson2PhGravity"
extern PetscErrorCode DefiantNewtonRaphson2PhGravity(BlackOilReservoirSimulation* MySim)
{

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DefiantNewtonRaphson2PhCapillaryPressure"
extern PetscErrorCode DefiantNewtonRaphson2PhCapillaryPressure(BlackOilReservoirSimulation* MySim)
{

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DefiantNewtonRaphson2PhHandleWells"
extern PetscErrorCode DefiantNewtonRaphson2PhHandleWells(BlackOilReservoirSimulation* MySim)
{

  PetscFunctionReturn(0);
}

/* assemble matrix and RHS */
#undef __FUNCT__
#define __FUNCT__ "DefiantNewtonRaphson2PhAssembleMatrix"
extern PetscErrorCode DefiantNewtonRaphson2PhAssembleMatrix(BlackOilReservoirSimulation* MySim)
{

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DefiantNewtonRaphson2PhAssembleRHS"
extern PetscErrorCode DefiantNewtonRaphson2PhAssembleRHS(BlackOilReservoirSimulation* MySim)
{

  PetscFunctionReturn(0);
}

/* Solver routines, can have many of those */
#undef __FUNCT__
#define __FUNCT__ "DefiantNewtonRaphson2PhDMMGComputeRHS"
extern PetscErrorCode DefiantNewtonRaphson2PhDMMGComputeRHS(DMMG dmmg, Vec b)
{

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DefiantNewtonRaphson2PhDMMGComputeMatrix"
extern PetscErrorCode DefiantNewtonRaphson2PhDMMGComputeMatrix(DMMG dmmg, Mat jac, Mat B)
{

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DefiantNewtonRaphson2PhSolve"
extern PetscErrorCode DefiantNewtonRaphson2PhSolve(BlackOilReservoirSimulation* MySim)
{

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DefiantNewtonRaphson2PhSolveDMMG"
extern PetscErrorCode DefiantNewtonRaphson2PhSolveDMMG(BlackOilReservoirSimulation* MySim)
{

  PetscFunctionReturn(0);
}

/* Iterating routines */
#undef __FUNCT__
#define __FUNCT__ "DefiantNewtonRaphson2PhIterate"
extern PetscErrorCode DefiantNewtonRaphson2PhIterate(BlackOilReservoirSimulation* MySim)
{

  PetscFunctionReturn(0);
}
