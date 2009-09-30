/*
 * RockProperties.c
 *
 *  Created on: Sep 9, 2009
 *      Author: yye00
 */

#include "Defiant.h"

#undef __FUNCT__
#define __FUNCT__ "DefiantUpdatePorosity"
extern PetscErrorCode DefiantUpdatePorosity(BlackOilReservoirSimulation* MySim)
{
  PetscErrorCode ierr;
  Vec TempVec;

  PetscFunctionBegin;

  ierr = VecDuplicate(MySim->Po, &TempVec);CHKERRQ(ierr);
  ierr = VecCopy(MySim->Po, TempVec);CHKERRQ(ierr);

  ierr = VecShift(TempVec, -1.0*MySim->RockCompRefPressure);CHKERRQ(ierr);
  ierr = VecScale(TempVec, MySim->RockCompressibility);CHKERRQ(ierr);
  ierr = VecShift(TempVec, 1.0);CHKERRQ(ierr);
  ierr = VecScale(TempVec, MySim->RockCompRefPorosity);CHKERRQ(ierr);

  /* copy back to porosity */
  ierr = VecCopy(TempVec, MySim->Phi);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
