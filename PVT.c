/*
 * PVT.c
 *
 *  Created on: Sep 6, 2009
 *      Author: yye00
 */

#include "Defiant.h"

#undef __FUNCT__
#define __FUNCT__ "DefiantUpdateDensity"
extern PetscErrorCode DefiantUpdateDensity(BlackOilReservoirSimulation* MySim)
{
  PetscFunctionBegin;

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DefiantUpdateVolumeFactors"
extern PetscErrorCode DefiantUpdateVolumeFactors(BlackOilReservoirSimulation* MySim)
{
  PetscErrorCode ierr;

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DefiantUpdateViscosities"
extern PetscErrorCode DefiantUpdateViscosities(BlackOilReservoirSimulation* MySim)
{
  PetscErrorCode ierr;

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DefiantUpdatePcow"
extern PetscErrorCode DefiantUpdatePcow(BlackOilReservoirSimulation* MySim)
{
  PetscFunctionBegin;
  PetscErrorCode ierr;
  PetscInt i, j, k, mx, my, mz, xm, ym, zm, xs, ys, zs;
  PetscScalar ***LocalFlowMask;
  /* Capillary pressure */
  PetscScalar ***LocalPcow;
  /* Water saturation */
  PetscScalar ***LocalSw;

  PetscFunctionBegin;

  /* Get dimensions and extents of the local vectors */
  ierr = DAGetInfo(MySim->SimDA, 0, &mx, &my, &mz, 0, 0, 0, 0, 0, 0, 0);CHKERRQ(ierr);
  ierr = DAGetCorners(MySim->SimDA, &xs, &ys, &zs, &xm, &ym, &zm);CHKERRQ(ierr);
  /* Grab the data for the flow field */
  ierr = DAVecGetArray(MySim->SimDA, MySim->FlowMask, &LocalFlowMask);CHKERRQ(ierr);
  /* Grab the local data for the capillary pressure */
  ierr = DAVecGetArray(MySim->SimDA, MySim->Pcow, &LocalPcow);CHKERRQ(ierr);
  /* Grab the local data for the Water saturation */
  ierr = DAVecGetArray(MySim->SimDA, MySim->Sw, &LocalSw);CHKERRQ(ierr);

  for (k = zs; k < zs + zm; k++) {
    for (j = ys; j < ys + ym; j++) {
      for (i = xs; i < xs + xm; i++) {
        if (i == 0 || j == 0 || k == 0 || i == mx - 1 || j == my - 1 || k == mz
            - 1) {
        } else if (ABS(LocalFlowMask[k][j][i]-FLUID_FLOW) < EPSILON) {
          DefiantTableXFromY(&MySim->PcowSw, &LocalPcow[k][j][i], &LocalSw[k][j][i]);
        }
      }
    }
  }

  /* Restore the new arrays to their rightful place */
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Pcow, &LocalPcow);CHKERRQ(ierr);

  /* Begin Assembly for vectors */
  ierr = VecAssemblyBegin(MySim->Pcow);CHKERRQ(ierr);
  /* And end Assembly */
  ierr = VecAssemblyEnd(MySim->Pcow);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DefiantUpdatePcog"
extern PetscErrorCode DefiantUpdatePcog(BlackOilReservoirSimulation* MySim)
{
  PetscFunctionBegin;
  PetscErrorCode ierr;
  PetscInt i, j, k, mx, my, mz, xm, ym, zm, xs, ys, zs;
  PetscScalar ***LocalFlowMask;
  /* Capillary pressure */
  PetscScalar ***LocalPcog;
  /* Water saturation */
  PetscScalar ***LocalSg;

  PetscFunctionBegin;

  /* Get dimensions and extents of the local vectors */
  ierr = DAGetInfo(MySim->SimDA, 0, &mx, &my, &mz, 0, 0, 0, 0, 0, 0, 0);CHKERRQ(ierr);
  ierr = DAGetCorners(MySim->SimDA, &xs, &ys, &zs, &xm, &ym, &zm);CHKERRQ(ierr);
  /* Grab the data for the flow field */
  ierr = DAVecGetArray(MySim->SimDA, MySim->FlowMask, &LocalFlowMask);CHKERRQ(ierr);
  /* Grab the local data for the capillary pressure */
  ierr = DAVecGetArray(MySim->SimDA, MySim->Pcog, &LocalPcog);CHKERRQ(ierr);
  /* Grab the local data for the Water saturation */
  ierr = DAVecGetArray(MySim->SimDA, MySim->Sg, &LocalSg);CHKERRQ(ierr);

  for (k = zs; k < zs + zm; k++) {
    for (j = ys; j < ys + ym; j++) {
      for (i = xs; i < xs + xm; i++) {
        if (i == 0 || j == 0 || k == 0 || i == mx - 1 || j == my - 1 || k == mz
            - 1) {
        } else if (ABS(LocalFlowMask[k][j][i]-FLUID_FLOW) < EPSILON) {
          DefiantTableXFromY(&MySim->PcogSg, &LocalPcog[k][j][i], &LocalSg[k][j][i]);
        }
      }
    }
  }

  /* Restore the new arrays to their rightful place */
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Pcog, &LocalPcog);CHKERRQ(ierr);

  /* Begin Assembly for vectors */
  ierr = VecAssemblyBegin(MySim->Pcog);CHKERRQ(ierr);
  /* And end Assembly */
  ierr = VecAssemblyEnd(MySim->Pcog);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DefiantUpdatePwFromPcow"
extern PetscErrorCode DefiantUpdatePwFromPcow(BlackOilReservoirSimulation* MySim)
{
  PetscErrorCode ierr;
  PetscScalar v;
  v = -1.0;

  PetscFunctionBegin;
  ierr = VecWAXPY(MySim->Pw, v, MySim->Pcow, MySim->Po);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DefiantUpdatePgFromPcog"
extern PetscErrorCode DefiantUpdatePgFromPcog(BlackOilReservoirSimulation* MySim)
{
  PetscErrorCode ierr;
  PetscScalar v;
  v = +1.0;

  PetscFunctionBegin;
  ierr = VecWAXPY(MySim->Pg, v, MySim->Pcog, MySim->Po);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
