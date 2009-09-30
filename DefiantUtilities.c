/*
 * DefiantUtilities.c
 *
 *  Created on: Sep 12, 2009
 *      Author: yye00
 */

#include "Defiant.h"


#undef __FUNCT__
#define __FUNCT__ "DefiantTableXFromY"
extern PetscErrorCode DefiantTableXFromY(DefiantTable *MyTable, PetscReal *X, PetscReal *Y)
{
  PetscErrorCode ierr;
  PetscInt i;

  PetscFunctionBegin;
  *X = -1.0;
  for (i=1; i<MyTable->NumberOfEntries; i++)
  {
    if(*Y > MyTable->Y[i-1] && *Y < MyTable->Y[i] )
      *X = MyTable->X[i-1]+(*Y-MyTable->Y[i-1])*(MyTable->X[i]-MyTable->X[i-1])/(MyTable->Y[i]-MyTable->Y[i-1]);
  }

  if(ABS(*X+1.0)<EPSILON) ierr = -1;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DefiantTableYFromX"
extern PetscErrorCode DefiantTableYFromX(DefiantTable *MyTable, PetscReal *X, PetscReal *Y)
{
  PetscErrorCode ierr;
  PetscInt i;

  PetscFunctionBegin;
  *X = -1.0;
  for (i=1; i<MyTable->NumberOfEntries; i++)
  {
    if(*X > MyTable->X[i-1] && *X < MyTable->X[i] )
      *Y = MyTable->Y[i-1]+(*X-MyTable->X[i-1])*(MyTable->Y[i]-MyTable->Y[i-1])/(MyTable->X[i]-MyTable->X[i-1]);
  }

  if(ABS(*Y+1.0)<EPSILON) ierr = -1;
  PetscFunctionReturn(0);
}
