/*
 * NewtonRaphson3Ph.c
 *
 *  Created on: Sep 20, 2009
 *      Author: yye00
 */

#include "Defiant.h"

#undef __FUNCT__
#define __FUNCT__ "DefiantNewtonRaphson3PhGravity"
extern PetscErrorCode DefiantNewtonRaphson3PhGravity(BlackOilReservoirSimulation* MySim)
{
  /* Here we compute the gravity terms Grav0, Grav1 and Grav 2 */
  /* for oil, water and Gas respectively */
  PetscErrorCode ierr;
  PetscInt i, j, k, mx, my, mz, xm, ym, zm, xs, ys, zs;
  PetscScalar ***LocalFlowMask;
  /* Depth */
  PetscScalar ***Localx3;
  /* Area * Permeability divided by height */
  PetscScalar ***LocalAKHx1m, ***LocalAKHx2m, ***LocalAKHx3m;
  PetscScalar ***LocalAKHx1p, ***LocalAKHx2p, ***LocalAKHx3p;
  /* Density By Viscosity at the faces */
  PetscScalar ***LocalRhoByMuox1m, ***LocalRhoByMuox1p, ***LocalRhoByMuwx1m, ***LocalRhoByMuwx1p, ***LocalRhoByMugx1m, ***LocalRhoByMugx1p;
  PetscScalar ***LocalRhoByMuox2m, ***LocalRhoByMuox2p, ***LocalRhoByMuwx2m, ***LocalRhoByMuwx2p, ***LocalRhoByMugx2m, ***LocalRhoByMugx2p;
  PetscScalar ***LocalRhoByMuox3m, ***LocalRhoByMuox3p, ***LocalRhoByMuwx3m, ***LocalRhoByMuwx3p, ***LocalRhoByMugx3m, ***LocalRhoByMugx3p;
  /* Relative Permeabilities at the faces */
  PetscScalar ***LocalRelPermox1m, ***LocalRelPermox1p, ***LocalRelPermox2m, ***LocalRelPermox2p, ***LocalRelPermox3m, ***LocalRelPermox3p;
  PetscScalar ***LocalRelPermwx1m, ***LocalRelPermwx1p, ***LocalRelPermwx2m, ***LocalRelPermwx2p, ***LocalRelPermwx3m, ***LocalRelPermwx3p;
  PetscScalar ***LocalRelPermgx1m, ***LocalRelPermgx1p, ***LocalRelPermgx2m, ***LocalRelPermgx2p, ***LocalRelPermgx3m, ***LocalRelPermgx3p;
  /* Gas solubility local values at faces */
  PetscScalar ***LocalRsox1p, ***LocalRsox1m, ***LocalRsox2p, ***LocalRsox2m, ***LocalRsox3p, ***LocalRsox3m;
  /* Gravity collected terms */
  PetscScalar ***LocalGrav0, ***LocalGrav1, ***LocalGrav2;

  PetscFunctionBegin;

  /* Get dimensions and extents of the local vectors */
  ierr = DAGetInfo(MySim->SimDA, 0, &mx, &my, &mz, 0, 0, 0, 0, 0, 0, 0);CHKERRQ(ierr);
  ierr = DAGetCorners(MySim->SimDA, &xs, &ys, &zs, &xm, &ym, &zm);CHKERRQ(ierr);
  /* Grab the data for the flow field */
  ierr = DAVecGetArray(MySim->SimDA, MySim->FlowMask, &LocalFlowMask);CHKERRQ(ierr);
  /* Grab the local data for depth */
  ierr = DAVecGetArray(MySim->SimDA, MySim->x3, &Localx3);CHKERRQ(ierr);
  /* Grab the local data for area * permeability divided by height */
  ierr = DAVecGetArray(MySim->SimDA, MySim->AKHx1m, &LocalAKHx1m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->AKHx2m, &LocalAKHx2m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->AKHx3m, &LocalAKHx3m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->AKHx1p, &LocalAKHx1p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->AKHx2p, &LocalAKHx2p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->AKHx3p, &LocalAKHx3p);CHKERRQ(ierr);
  /* Density divided by viscosity at the faces */
  ierr = DAVecGetArray(MySim->SimDA, MySim->RhoByMuox1m, &LocalRhoByMuox1m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->RhoByMuox1p, &LocalRhoByMuox1p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->RhoByMuwx1m, &LocalRhoByMuwx1m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->RhoByMuwx1p, &LocalRhoByMuwx1p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->RhoByMuwx1m, &LocalRhoByMugx1m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->RhoByMuwx1p, &LocalRhoByMugx1p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->RhoByMuox2m, &LocalRhoByMuox2m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->RhoByMuox2p, &LocalRhoByMuox2p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->RhoByMuwx2m, &LocalRhoByMuwx2m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->RhoByMuwx2p, &LocalRhoByMuwx2p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->RhoByMuwx2m, &LocalRhoByMugx2m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->RhoByMuwx2p, &LocalRhoByMugx2p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->RhoByMuox3m, &LocalRhoByMuox3m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->RhoByMuox3p, &LocalRhoByMuox3p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->RhoByMuwx3m, &LocalRhoByMuwx3m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->RhoByMuwx3p, &LocalRhoByMuwx3p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->RhoByMuwx3m, &LocalRhoByMugx3m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->RhoByMuwx3p, &LocalRhoByMugx3p);CHKERRQ(ierr);
  /* Relative permeabilities at the faces */
  ierr = DAVecGetArray(MySim->SimDA, MySim->RelPermox1m, &LocalRelPermox1m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->RelPermox1p, &LocalRelPermox1p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->RelPermox2m, &LocalRelPermox2m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->RelPermox2p, &LocalRelPermox2p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->RelPermox3m, &LocalRelPermox3m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->RelPermox3p, &LocalRelPermox3p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->RelPermwx1m, &LocalRelPermwx1m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->RelPermwx1p, &LocalRelPermwx1p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->RelPermwx2m, &LocalRelPermwx2m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->RelPermwx2p, &LocalRelPermwx2p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->RelPermwx3m, &LocalRelPermwx3m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->RelPermwx3p, &LocalRelPermwx3p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->RelPermgx1m, &LocalRelPermgx1m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->RelPermgx1p, &LocalRelPermgx1p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->RelPermgx2m, &LocalRelPermgx2m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->RelPermgx2p, &LocalRelPermgx2p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->RelPermgx3m, &LocalRelPermgx3m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->RelPermgx3p, &LocalRelPermgx3p);CHKERRQ(ierr);
  /* Gas solubility at the faces */
  ierr = DAVecGetArray(MySim->SimDA, MySim->Rsox1m, &LocalRsox1p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Rsox1p, &LocalRsox1m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Rsox2m, &LocalRsox2p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Rsox2p, &LocalRsox2m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Rsox3m, &LocalRsox3p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Rsox3p, &LocalRsox3m);CHKERRQ(ierr);
  /* Gravity collected terms */
  ierr = DAVecGetArray(MySim->SimDA, MySim->Grav0, &LocalGrav0);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Grav1, &LocalGrav1);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Grav2, &LocalGrav2);CHKERRQ(ierr);

  for (k = zs; k < zs + zm; k++) {
    for (j = ys; j < ys + ym; j++) {
      for (i = xs; i < xs + xm; i++) {
        if (i == 0 || j == 0 || k == 0 || i == mx - 1 || j == my - 1 || k == mz
            - 1) {
        } else if (ABS(LocalFlowMask[k][j][i]-FLUID_FLOW) < EPSILON) {
          LocalGrav0[k][j][i] = MySim->GravAcc * (
                LocalAKHx1p[k][j][i] * LocalRelPermox1p[k][j][i] * LocalRhoByMuox1p[k][j][i] * (Localx3[k][j][i + 1] - Localx3[k][j][i])
              + LocalAKHx2p[k][j][i] * LocalRelPermox2p[k][j][i] * LocalRhoByMuox2p[k][j][i] * (Localx3[k][j + 1][i] - Localx3[k][j][i])
              + LocalAKHx3p[k][j][i] * LocalRelPermox3p[k][j][i] * LocalRhoByMuox3p[k][j][i] * (Localx3[k + 1][j][i] - Localx3[k][j][i])
              - LocalAKHx1m[k][j][i] * LocalRelPermox1m[k][j][i] * LocalRhoByMuox1m[k][j][i] * (Localx3[k][j][i] - Localx3[k][j][i - 1])
              - LocalAKHx2m[k][j][i] * LocalRelPermox2m[k][j][i] * LocalRhoByMuox2m[k][j][i] * (Localx3[k][j][i] - Localx3[k][j - 1][i])
              - LocalAKHx3m[k][j][i] * LocalRelPermox3m[k][j][i] * LocalRhoByMuox3m[k][j][i] * (Localx3[k][j][i] - Localx3[k - 1][j][i])   );

          LocalGrav1[k][j][i] = MySim->GravAcc * (
                LocalAKHx1p[k][j][i] * LocalRelPermwx1p[k][j][i] * LocalRhoByMuwx1p[k][j][i] * (Localx3[k][j][i + 1] - Localx3[k][j][i])
              + LocalAKHx2p[k][j][i] * LocalRelPermwx2p[k][j][i] * LocalRhoByMuwx2p[k][j][i] * (Localx3[k][j + 1][i] - Localx3[k][j][i])
              + LocalAKHx3p[k][j][i] * LocalRelPermwx3p[k][j][i] * LocalRhoByMuwx3p[k][j][i] * (Localx3[k + 1][j][i] - Localx3[k][j][i])
              - LocalAKHx1m[k][j][i] * LocalRelPermwx1m[k][j][i] * LocalRhoByMuwx1m[k][j][i] * (Localx3[k][j][i] - Localx3[k][j][i - 1])
              - LocalAKHx2m[k][j][i] * LocalRelPermwx2m[k][j][i] * LocalRhoByMuwx2m[k][j][i] * (Localx3[k][j][i] - Localx3[k][j - 1][i])
              - LocalAKHx3m[k][j][i] * LocalRelPermwx3m[k][j][i] * LocalRhoByMuwx3m[k][j][i] * (Localx3[k][j][i] - Localx3[k - 1][j][i])    );

          LocalGrav2[k][j][i] = MySim->GravAcc * (
              + LocalAKHx1p[k][j][i] * LocalRelPermgx1p[k][j][i] * LocalRhoByMugx1p[k][j][i] * (Localx3[k][j][i + 1] - Localx3[k][j][i])
              + LocalAKHx2p[k][j][i] * LocalRelPermgx2p[k][j][i] * LocalRhoByMugx2p[k][j][i] * (Localx3[k][j + 1][i] - Localx3[k][j][i])
              + LocalAKHx3p[k][j][i] * LocalRelPermgx3p[k][j][i] * LocalRhoByMugx3p[k][j][i] * (Localx3[k + 1][j][i] - Localx3[k][j][i])
              - LocalAKHx1m[k][j][i] * LocalRelPermgx1m[k][j][i] * LocalRhoByMugx1m[k][j][i] * (Localx3[k][j][i] - Localx3[k][j][i - 1])
              - LocalAKHx2m[k][j][i] * LocalRelPermgx2m[k][j][i] * LocalRhoByMugx2m[k][j][i] * (Localx3[k][j][i] - Localx3[k][j - 1][i])
              - LocalAKHx3m[k][j][i] * LocalRelPermgx3m[k][j][i] * LocalRhoByMugx3m[k][j][i] * (Localx3[k][j][i] - Localx3[k - 1][j][i])

              + LocalRsox1p[k][j][i] * LocalAKHx1p[k][j][i] * LocalRelPermox1p[k][j][i] * LocalRhoByMuox1p[k][j][i] * (Localx3[k][j][i + 1] - Localx3[k][j][i])
              + LocalRsox2p[k][j][i] * LocalAKHx2p[k][j][i] * LocalRelPermox2p[k][j][i] * LocalRhoByMuox2p[k][j][i] * (Localx3[k][j + 1][i] - Localx3[k][j][i])
              + LocalRsox3p[k][j][i] * LocalAKHx3p[k][j][i] * LocalRelPermox3p[k][j][i] * LocalRhoByMuox3p[k][j][i] * (Localx3[k + 1][j][i] - Localx3[k][j][i])
              - LocalRsox1m[k][j][i] * LocalAKHx1m[k][j][i] * LocalRelPermox1m[k][j][i] * LocalRhoByMuox1m[k][j][i] * (Localx3[k][j][i] - Localx3[k][j][i - 1])
              - LocalRsox2m[k][j][i] * LocalAKHx2m[k][j][i] * LocalRelPermox2m[k][j][i] * LocalRhoByMuox2m[k][j][i] * (Localx3[k][j][i] - Localx3[k][j - 1][i])
              - LocalRsox3m[k][j][i] * LocalAKHx3m[k][j][i] * LocalRelPermox3m[k][j][i] * LocalRhoByMuox3m[k][j][i] * (Localx3[k][j][i] - Localx3[k - 1][j][i])  );

        }
      }
    }
  }

  /* Restore the new arrays to their rightful place */
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Grav0, &LocalGrav0);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Grav1, &LocalGrav1);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Grav2, &LocalGrav2);CHKERRQ(ierr);

  /* Begin Assembly for vectors */
  ierr = VecAssemblyBegin(MySim->Grav0);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->Grav1);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->Grav2);CHKERRQ(ierr);
  /* And end Assembly */
  ierr = VecAssemblyEnd(MySim->Grav0);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->Grav1);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->Grav2);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DefiantNewtonRaphson3PhCapillaryPressure"
extern PetscErrorCode DefiantNewtonRaphson3PhCapillaryPressure(BlackOilReservoirSimulation* MySim)
{
  /* Here we compute the capillary pressure terms for oil water  and oil gas */
  PetscErrorCode ierr;
  PetscInt i, j, k, mx, my, mz, xm, ym, zm, xs, ys, zs;
  PetscScalar ***LocalFlowMask;
  /* Transmissibility For Water */
  PetscScalar ***LocalTwx1m, ***LocalTwx2m, ***LocalTwx3m;
  PetscScalar ***LocalTwx1p, ***LocalTwx2p, ***LocalTwx3p;
  /* Transmissibility For Gas */
  PetscScalar ***LocalTgx1m, ***LocalTgx2m, ***LocalTgx3m;
  PetscScalar ***LocalTgx1p, ***LocalTgx2p, ***LocalTgx3p;
  /* Capillary pressure */
  PetscScalar ***LocalPcow, ***LocalPcog;
  /* Capillary pressure collected term */
  PetscScalar ***LocalNRCapPressure, ***LocalNRGasCapPressure;

  PetscFunctionBegin;

  /* Get dimensions and extents of the local vectors */
  ierr = DAGetInfo(MySim->SimDA, 0, &mx, &my, &mz, 0, 0, 0, 0, 0, 0, 0);CHKERRQ(ierr);
  ierr = DAGetCorners(MySim->SimDA, &xs, &ys, &zs, &xm, &ym, &zm);CHKERRQ(ierr);
  /* Grab the data for the flow field */
  ierr = DAVecGetArray(MySim->SimDA, MySim->FlowMask, &LocalFlowMask);CHKERRQ(ierr);
  /* Grab the local data for the water transmissibilities */
  ierr = DAVecGetArray(MySim->SimDA, MySim->Twx1m, &LocalTwx1m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Twx1p, &LocalTwx1p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Twx2m, &LocalTwx2m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Twx2p, &LocalTwx2p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Twx3m, &LocalTwx3m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Twx3p, &LocalTwx3p);CHKERRQ(ierr);
  /* Grab the local data for the gas transmissibilities */
  ierr = DAVecGetArray(MySim->SimDA, MySim->Tgx1m, &LocalTgx1m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Tgx1p, &LocalTgx1p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Tgx2m, &LocalTgx2m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Tgx2p, &LocalTgx2p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Tgx3m, &LocalTgx3m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Tgx3p, &LocalTgx3p);CHKERRQ(ierr);
  /* Grab the local data for the capillary pressures */
  ierr = DAVecGetArray(MySim->SimDA, MySim->Pcow, &LocalPcow);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Pcog, &LocalPcog);CHKERRQ(ierr);
  /* Grab the local data for the collected capillary pressure */
  ierr = DAVecGetArray(MySim->SimDA, MySim->NRCapPressure, &LocalNRCapPressure);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->NRGasCapPressure, &LocalNRGasCapPressure);CHKERRQ(ierr);

  for (k = zs; k < zs + zm; k++) {
    for (j = ys; j < ys + ym; j++) {
      for (i = xs; i < xs + xm; i++) {
        if (i == 0 || j == 0 || k == 0 || i == mx - 1 || j == my - 1 || k == mz
            - 1) {
        } else if (ABS(LocalFlowMask[k][j][i]-FLUID_FLOW) < EPSILON) {
          LocalNRCapPressure[k][j][i] =
              LocalTwx1p[k][j][i] * (LocalPcow[k][j][i + 1] - LocalPcow[k][j][i])
            + LocalTwx2p[k][j][i] * (LocalPcow[k][j + 1][i] - LocalPcow[k][j][i])
            + LocalTwx3p[k][j][i] * (LocalPcow[k + 1][j][i] - LocalPcow[k][j][i])
            - LocalTwx1m[k][j][i] * (LocalPcow[k][j][i] - LocalPcow[k][j][i - 1])
            - LocalTwx2m[k][j][i] * (LocalPcow[k][j][i] - LocalPcow[k][j - 1][i])
            - LocalTwx3m[k][j][i] * (LocalPcow[k][j][i] - LocalPcow[k - 1][j][i])  ;

          LocalNRGasCapPressure[k][j][i] =
              LocalTgx1p[k][j][i] * (LocalPcog[k][j][i + 1] - LocalPcog[k][j][i])
            - LocalTgx2p[k][j][i] * (LocalPcog[k][j + 1][i] - LocalPcog[k][j][i])
            - LocalTgx3p[k][j][i] * (LocalPcog[k + 1][j][i] - LocalPcog[k][j][i])
            + LocalTgx1m[k][j][i] * (LocalPcog[k][j][i] - LocalPcog[k][j][i - 1])
            + LocalTgx2m[k][j][i] * (LocalPcog[k][j][i] - LocalPcog[k][j - 1][i])
            + LocalTgx3m[k][j][i] * (LocalPcog[k][j][i] - LocalPcog[k - 1][j][i])  ;
        }
      }
    }
  }

  /* Restore the new arrays to their rightful place */
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->NRCapPressure, &LocalNRCapPressure);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->NRGasCapPressure, &LocalNRGasCapPressure);CHKERRQ(ierr);

  /* Begin Assembly for vectors */
  ierr = VecAssemblyBegin(MySim->NRCapPressure);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->NRGasCapPressure);CHKERRQ(ierr);
  /* And end Assembly */
  ierr = VecAssemblyEnd(MySim->NRCapPressure);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->NRGasCapPressure);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DefiantNewtonRaphson3PhHandleWells"
extern PetscErrorCode DefiantNewtonRaphson3PhHandleWells(BlackOilReservoirSimulation* MySim)
{
  PetscErrorCode ierr;
  PetscInt mx, my, mz, xm, ym, zm, xs, ys, zs;
  PetscScalar re;
  /* Well handling variables */
  PetscInt PerfIDMine, PerfIDOther, WellID;
  PetscScalar qo, qw, qg;
  PetscInt MyI, MyJ, MyK;
  PetscInt OtherI, OtherJ, OtherK;
  /* Local values for variables*/
  PetscScalar ***LocalFlowMask;
  PetscScalar ***Localx3;
  PetscScalar ***Localh1,   ***Localh2,   ***Localh3;
  PetscScalar ***LocalMuo,  ***LocalMuw,  ***LocalMug;
  PetscScalar ***LocalRhoo, ***LocalRhow, ***LocalRhog;
  PetscScalar ***LocalK11,  ***LocalK22,  ***LocalK33;
  PetscScalar ***LocalKro,  ***LocalKrw,  ***LocalKrg;
  PetscScalar ***LocalPo,   ***LocalPw,   ***LocalPg;
  PetscScalar ***LocalPcow, ***LocalPcog;
  PetscScalar ***LocalBo,   ***LocalBw,   ***LocalBg;
  PetscScalar ***LocalQo,   ***LocalQw,   ***LocalQg;

  PetscFunctionBegin;
  /* Get dimensions and extents of the local vectors */
  ierr = DAGetInfo(MySim->SimDA, 0, &mx, &my, &mz, 0, 0, 0, 0, 0, 0, 0);CHKERRQ(ierr);
  ierr = DAGetCorners(MySim->SimDA, &xs, &ys, &zs, &xm, &ym, &zm);CHKERRQ(ierr);
  /* Grab the data for the flow field */
  ierr = DAVecGetArray(MySim->SimDA, MySim->FlowMask, &LocalFlowMask);CHKERRQ(ierr);
  /* Grab the local depth */
  ierr = DAVecGetArray(MySim->SimDA, MySim->x3, &Localx3);CHKERRQ(ierr);
  /* Grab the local geometry */
  ierr = DAVecGetArray(MySim->SimDA, MySim->h1, &Localh1);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->h2, &Localh2);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->h3, &Localh3);CHKERRQ(ierr);
  /* Grab the Densities */
  ierr = DAVecGetArray(MySim->SimDA, MySim->Rhoo, &LocalRhoo);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Rhow, &LocalRhow);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Rhog, &LocalRhog);CHKERRQ(ierr);
  /* Grab the local viscosities*/
  ierr = DAVecGetArray(MySim->SimDA, MySim->Muo, &LocalMuo);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Muw, &LocalMuw);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Mug, &LocalMug);CHKERRQ(ierr);
  /* Grab the local permeabilities */
  ierr = DAVecGetArray(MySim->SimDA, MySim->K11, &LocalK11);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->K22, &LocalK22);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->K33, &LocalK33);CHKERRQ(ierr);
  /* Grab the local permeabilities */
  ierr = DAVecGetArray(MySim->SimDA, MySim->Kro, &LocalKro);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Krw, &LocalKrw);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Krg, &LocalKrg);CHKERRQ(ierr);
  /* Grab the local data for pressures */
  ierr = DAVecGetArray(MySim->SimDA, MySim->Po, &LocalPo);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Pw, &LocalPw);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Pg, &LocalPg);CHKERRQ(ierr);
  /* Grab the local capillary pressure */
  ierr = DAVecGetArray(MySim->SimDA, MySim->Pcow, &LocalPcow);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Pcog, &LocalPcog);CHKERRQ(ierr);
  /* Grab the local data for volume factors at the cell centers */
  ierr = DAVecGetArray(MySim->SimDA, MySim->Bo, &LocalBo);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Bw, &LocalBw);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Bg, &LocalBg);CHKERRQ(ierr);
  /* Grab the local flow rate */
  ierr = DAVecGetArray(MySim->SimDA, MySim->Qo, &LocalQo);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Qw, &LocalQw);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Qg, &LocalQg);CHKERRQ(ierr);

  for (WellID = 0; WellID < MySim->NumberOfWells; WellID++) {
    for (PerfIDMine = 0; PerfIDMine
        < (MySim->Wells[WellID]).NumberOfPerforations; PerfIDMine++) {
      if (MySim->Wells[WellID].Perforations[PerfIDMine].Constraint == FLOW_RATE_CONSTRAINT){
        MySim->Wells[WellID].Perforations[PerfIDMine].WellIndex = -1.0;
      } else {
        MyI = (MySim->Wells[WellID]).Perforations[PerfIDMine].I;
        MyJ = (MySim->Wells[WellID]).Perforations[PerfIDMine].J;
        MyK = (MySim->Wells[WellID]).Perforations[PerfIDMine].K;
        if ((MySim->Wells[WellID]).Perforations[PerfIDMine].Orientation
            == PERF_ORIENTATION_X1X1) {
          re = 0.14 / 0.5 * sqrt((sqrt(LocalK33[MyK][MyJ][MyI]
              / LocalK22[MyK][MyJ][MyI]) * Localh2[MyK][MyJ][MyI]
              * Localh2[MyK][MyJ][MyI] + sqrt(LocalK22[MyK][MyJ][MyI]
              / LocalK33[MyK][MyJ][MyI]) * Localh3[MyK][MyJ][MyI]
              * Localh3[MyK][MyJ][MyI])) / (pow(LocalK33[MyK][MyJ][MyI]
              / LocalK22[MyK][MyJ][MyI], 0.25) + pow(LocalK22[MyK][MyJ][MyI]
              / LocalK33[MyK][MyJ][MyI], 0.25));
          MySim->Wells[WellID].Perforations[PerfIDMine].WellIndex = 2.0 * PI
              * Localh1[MyK][MyJ][MyI] * sqrt(LocalK11[MyK][MyJ][MyI]
              * LocalK33[MyK][MyJ][MyI]) / (log(re
              / MySim->Wells[WellID].Perforations[PerfIDMine].Rw)
              + MySim->Wells[WellID].Perforations[PerfIDMine].S);
        } else if ((MySim->Wells[WellID]).Perforations[PerfIDMine].Orientation
            == PERF_ORIENTATION_X2X2) {
          re = 0.14 / 0.5 * sqrt((sqrt(LocalK33[MyK][MyJ][MyI]
              / LocalK11[MyK][MyJ][MyI]) * Localh1[MyK][MyJ][MyI]
              * Localh1[MyK][MyJ][MyI] + sqrt(LocalK11[MyK][MyJ][MyI]
              / LocalK33[MyK][MyJ][MyI]) * Localh3[MyK][MyJ][MyI]
              * Localh3[MyK][MyJ][MyI])) / (pow(LocalK33[MyK][MyJ][MyI]
              / LocalK11[MyK][MyJ][MyI], 0.25) + pow(LocalK11[MyK][MyJ][MyI]
              / LocalK33[MyK][MyJ][MyI], 0.25));
          MySim->Wells[WellID].Perforations[PerfIDMine].WellIndex = 2.0 * PI
              * Localh2[MyK][MyJ][MyI] * sqrt(LocalK11[MyK][MyJ][MyI]
              * LocalK33[MyK][MyJ][MyI]) / (log(re
              / MySim->Wells[WellID].Perforations[PerfIDMine].Rw)
              + MySim->Wells[WellID].Perforations[PerfIDMine].S);
        } else if ((MySim->Wells[WellID]).Perforations[PerfIDMine].Orientation
            == PERF_ORIENTATION_X3X3) {
          re = 0.14 / 0.5 * sqrt((sqrt(LocalK22[MyK][MyJ][MyI]
              / LocalK11[MyK][MyJ][MyI]) * Localh1[MyK][MyJ][MyI]
              * Localh1[MyK][MyJ][MyI] + sqrt(LocalK11[MyK][MyJ][MyI]
              / LocalK22[MyK][MyJ][MyI]) * Localh2[MyK][MyJ][MyI]
              * Localh2[MyK][MyJ][MyI])) / (pow(LocalK22[MyK][MyJ][MyI]
              / LocalK11[MyK][MyJ][MyI], 0.25) + pow(LocalK11[MyK][MyJ][MyI]
              / LocalK22[MyK][MyJ][MyI], 0.25));
          MySim->Wells[WellID].Perforations[PerfIDMine].WellIndex = 2.0 * PI
              * Localh3[MyK][MyJ][MyI] * sqrt(LocalK11[MyK][MyJ][MyI]
              * LocalK22[MyK][MyJ][MyI]) / (log(re
              / MySim->Wells[WellID].Perforations[PerfIDMine].Rw)
              + MySim->Wells[WellID].Perforations[PerfIDMine].S);
        }
      }
    }
  }

  for (WellID = 0; WellID < MySim->NumberOfWells; WellID++) {
    for (PerfIDMine = 0; PerfIDMine
        < (MySim->Wells[WellID]).NumberOfPerforations; PerfIDMine++) {
      if (MySim->Wells[WellID].Perforations[PerfIDMine].IsActive == PETSC_TRUE){
        /* grab my location */
        MyI = (MySim->Wells[WellID]).Perforations[PerfIDMine].I;
        MyJ = (MySim->Wells[WellID]).Perforations[PerfIDMine].J;
        MyK = (MySim->Wells[WellID]).Perforations[PerfIDMine].K;

        if (MySim->Wells[WellID].Perforations[PerfIDMine].Constraint == FLOW_RATE_CONSTRAINT){
          if (MySim->Wells[WellID].Perforations[PerfIDMine].WellType == OIL_INJECTOR )
            LocalQo[MyK][MyJ][MyI] = (MySim->Wells[WellID]).Perforations[PerfIDMine].Qo;
          else if (MySim->Wells[WellID].Perforations[PerfIDMine].WellType == WATER_INJECTOR )
            LocalQw[MyK][MyJ][MyI] = (MySim->Wells[WellID]).Perforations[PerfIDMine].Qw;
          else if (MySim->Wells[WellID].Perforations[PerfIDMine].WellType == GAS_INJECTOR )
            LocalQg[MyK][MyJ][MyI] = (MySim->Wells[WellID]).Perforations[PerfIDMine].Qg;

          if (MySim->Wells[WellID].Perforations[PerfIDMine].WellType == OIL_PRODUCER )
            LocalQo[MyK][MyJ][MyI] = (MySim->Wells[WellID]).Perforations[PerfIDMine].Qo;
          else if (MySim->Wells[WellID].Perforations[PerfIDMine].WellType == WATER_PRODUCER )
            LocalQw[MyK][MyJ][MyI] = (MySim->Wells[WellID]).Perforations[PerfIDMine].Qw;
          else if (MySim->Wells[WellID].Perforations[PerfIDMine].WellType == GAS_PRODUCER )
            LocalQg[MyK][MyJ][MyI] = (MySim->Wells[WellID]).Perforations[PerfIDMine].Qg;

        } else if (MySim->Wells[WellID].Perforations[PerfIDMine].Constraint == BHP_CONSTRAINT) {
          if ((MySim->Wells[WellID]).NumberOfPerforations == 1 ) {
            qo  = (MySim->Wells[WellID]).Perforations[PerfIDMine].WellIndex
                * LocalKro[MyK][MyJ][MyI] / LocalMuo[MyK][MyJ][MyI]
                * ((MySim->Wells[WellID]).Perforations[PerfIDMine].BHPo
                - LocalPo[MyK][MyJ][MyI]
                - MySim->GravAcc * LocalRhoo[MyK][MyJ][MyI]
                * ((MySim->Wells[WellID]).Perforations[PerfIDMine].zbh
                - Localx3[MyK][MyJ][MyI])) * Localh1[MyK][MyJ][MyI]
                * Localh2[MyK][MyJ][MyI] * Localh3[MyK][MyJ][MyI]
                / LocalBo[MyK][MyJ][MyI];
            qw  = (MySim->Wells[WellID]).Perforations[PerfIDMine].WellIndex
                * LocalKrw[MyK][MyJ][MyI] / LocalMuw[MyK][MyJ][MyI]
                * ((MySim->Wells[WellID]).Perforations[PerfIDMine].BHPw
                - LocalPw[MyK][MyJ][MyI]
                - MySim->GravAcc * LocalRhow[MyK][MyJ][MyI]
                * ((MySim->Wells[WellID]).Perforations[PerfIDMine].zbh
                - Localx3[MyK][MyJ][MyI])) * Localh1[MyK][MyJ][MyI]
                * Localh2[MyK][MyJ][MyI] * Localh3[MyK][MyJ][MyI]
                / LocalBw[MyK][MyJ][MyI];
            qg  = (MySim->Wells[WellID]).Perforations[PerfIDMine].WellIndex
                * LocalKrg[MyK][MyJ][MyI] / LocalMug[MyK][MyJ][MyI]
                * ((MySim->Wells[WellID]).Perforations[PerfIDMine].BHPg
                - LocalPg[MyK][MyJ][MyI]
                - MySim->GravAcc * LocalRhog[MyK][MyJ][MyI]
                * ((MySim->Wells[WellID]).Perforations[PerfIDMine].zbh
                - Localx3[MyK][MyJ][MyI])) * Localh1[MyK][MyJ][MyI]
                * Localh2[MyK][MyJ][MyI] * Localh3[MyK][MyJ][MyI]
                / LocalBg[MyK][MyJ][MyI];

            /* Store the flow rates and relevant data somewhere useful */
            LocalQo[MyK][MyJ][MyI] = qo;
            LocalQw[MyK][MyJ][MyI] = qw;
            LocalQg[MyK][MyJ][MyI] = qg;
          } else if ((MySim->Wells[WellID]).NumberOfPerforations > 1) {
            /* Set qo and qw to totals, we will be subtracting perfs from other cells from those values */
            qo = MySim->Wells[WellID].TotalQo;
            qw = MySim->Wells[WellID].TotalQw;
            qg = MySim->Wells[WellID].TotalQg;
            for (PerfIDOther = 0; PerfIDOther
                < (MySim->Wells[WellID]).NumberOfPerforations; PerfIDOther++) {
              OtherI = (MySim->Wells[WellID]).Perforations[PerfIDOther].I;
              OtherJ = (MySim->Wells[WellID]).Perforations[PerfIDOther].J;
              OtherK = (MySim->Wells[WellID]).Perforations[PerfIDOther].K;
              if (PerfIDOther != PerfIDMine && MySim->Wells[WellID].Perforations[PerfIDOther].IsActive == PETSC_TRUE) {
                qo  = qo - (MySim->Wells[WellID]).Perforations[PerfIDOther].WellIndex
                    * LocalKro[OtherK][OtherJ][OtherI] / LocalMuo[OtherK][OtherJ][OtherI]
                    * ((MySim->Wells[WellID]).Perforations[PerfIDOther].BHPo
                    - LocalPo[OtherK][OtherJ][OtherI]
                    - MySim->GravAcc * LocalRhoo[OtherK][OtherJ][OtherI]
                    * ((MySim->Wells[WellID]).Perforations[PerfIDOther].zbh
                    - Localx3[OtherK][OtherJ][OtherI])) * Localh1[OtherK][OtherJ][OtherI]
                    * Localh2[OtherK][OtherJ][OtherI] * Localh3[OtherK][OtherJ][OtherI]
                    / LocalBo[OtherK][OtherJ][OtherI];
                qw  = qw - (MySim->Wells[WellID]).Perforations[PerfIDOther].WellIndex
                    * LocalKrw[OtherK][OtherJ][OtherI] / LocalMuw[OtherK][OtherJ][OtherI]
                    * ((MySim->Wells[WellID]).Perforations[PerfIDOther].BHPw
                    - LocalPw[OtherK][OtherJ][OtherI]
                    - MySim->GravAcc * LocalRhow[OtherK][OtherJ][OtherI]
                    * ((MySim->Wells[WellID]).Perforations[PerfIDOther].zbh
                    - Localx3[OtherK][OtherJ][OtherI])) * Localh1[OtherK][OtherJ][OtherI]
                    * Localh2[OtherK][OtherJ][OtherI] * Localh3[OtherK][OtherJ][OtherI]
                    / LocalBw[OtherK][OtherJ][OtherI];
                qg  = qg - (MySim->Wells[WellID]).Perforations[PerfIDOther].WellIndex
                    * LocalKrg[OtherK][OtherJ][OtherI] / LocalMug[OtherK][OtherJ][OtherI]
                    * ((MySim->Wells[WellID]).Perforations[PerfIDOther].BHPg
                    - LocalPg[OtherK][OtherJ][OtherI]
                    - MySim->GravAcc * LocalRhog[OtherK][OtherJ][OtherI]
                    * ((MySim->Wells[WellID]).Perforations[PerfIDOther].zbh
                    - Localx3[OtherK][OtherJ][OtherI])) * Localh1[OtherK][OtherJ][OtherI]
                    * Localh2[OtherK][OtherJ][OtherI] * Localh3[OtherK][OtherJ][OtherI]
                    / LocalBg[OtherK][OtherJ][OtherI];

              }
            }
            /* now that I have the final qo and qw add them to the RHS */
            /* Add values to the RHS for the same column */
            LocalQo[MyK][MyJ][MyI] = qo;
            LocalQw[MyK][MyJ][MyI] = qw;
            LocalQg[MyK][MyJ][MyI] = qg;
          }
        }
      }
    }
  }

  /* Restore the Qo, Qw and Qg */
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Qo,&LocalQo);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Qw,&LocalQw);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Qg,&LocalQg);CHKERRQ(ierr);
  /* Begin assembly the vectors */
  ierr = VecAssemblyBegin(MySim->Qo);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->Qw);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->Qg);CHKERRQ(ierr);
  /* End assembly the vectors */
  ierr = VecAssemblyEnd(MySim->Qo);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->Qw);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->Qg);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/* assemble matrix and RHS */
#undef __FUNCT__
#define __FUNCT__ "DefiantNewtonRaphson3PhAssembleMatrix"
extern PetscErrorCode DefiantNewtonRaphson3PhAssembleMatrix(BlackOilReservoirSimulation* MySim)
{
  /* Assemble the big-ass matrix for fully implicit */
  /* ALWAYS ASSEMBLE BEFOPRE YOU HANDLE WELLS */
  PetscErrorCode ierr;
  PetscInt i, j, k, mx, my, mz, xm, ym, zm, xs, ys, zs;
  PetscInt       M,N,P,m,n,p;
  /* Transmissibility For Oil */
  PetscScalar ***LocalTox1m, ***LocalTox2m, ***LocalTox3m;
  PetscScalar ***LocalTox1p, ***LocalTox2p, ***LocalTox3p;
  PetscScalar ***LocalFlowMask;
  /* Transmissibility For Water */
  PetscScalar ***LocalTwx1m, ***LocalTwx2m, ***LocalTwx3m;
  PetscScalar ***LocalTwx1p, ***LocalTwx2p, ***LocalTwx3p;
  /* Transmissibility For Gas */
  PetscScalar ***LocalTgx1m, ***LocalTgx2m, ***LocalTgx3m;
  PetscScalar ***LocalTgx1p, ***LocalTgx2p, ***LocalTgx3p;
  /* Capillary pressure collected terms */
  PetscScalar ***LocalNRCapPressure, ***LocalNRGasCapPressure;
  /* Gravity collected terms */
  PetscScalar ***LocalGrav0, ***LocalGrav1, ***LocalGrav2;
  /* Flow rates for oil, water, gas */
  PetscScalar ***LocalQo, ***LocalQw, ***LocalQg;
  /* Volume factors for oil, water and gas at cell centers */
  PetscScalar ***LocalBo, ***LocalBw, ***LocalBg;
  /* Gas oil solubility */
  PetscScalar ***LocalRso;
  /* dimensions  for volume */
  PetscScalar ***Localh1, ***Localh2, ***Localh3;
  /* The RHS */
  PetscScalar ***LocalRHS0, ***LocalRHS1, ***LocalRHS2;

  PetscFunctionBegin;
  /* Grab the local data for the oil transmissibilities */
  ierr = DAVecGetArray(MySim->SimDA, MySim->Tox1m, &LocalTgx1m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Tox1p, &LocalTgx1p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Tox2m, &LocalTgx2m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Tox2p, &LocalTgx2p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Tox3m, &LocalTgx3m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Tox3p, &LocalTgx3p);CHKERRQ(ierr);
  /* Grab the local data for the water transmissibilities */
  ierr = DAVecGetArray(MySim->SimDA, MySim->Twx1m, &LocalTwx1m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Twx1p, &LocalTwx1p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Twx2m, &LocalTwx2m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Twx2p, &LocalTwx2p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Twx3m, &LocalTwx3m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Twx3p, &LocalTwx3p);CHKERRQ(ierr);
  /* Grab the local data for the gas transmissibilities */
  ierr = DAVecGetArray(MySim->SimDA, MySim->Tgx1m, &LocalTgx1m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Tgx1p, &LocalTgx1p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Tgx2m, &LocalTgx2m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Tgx2p, &LocalTgx2p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Tgx3m, &LocalTgx3m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Tgx3p, &LocalTgx3p);CHKERRQ(ierr);
  /* Grab the local data for the collected capillary pressure */
  ierr = DAVecGetArray(MySim->SimDA, MySim->NRCapPressure, &LocalNRCapPressure);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->NRGasCapPressure, &LocalNRGasCapPressure);CHKERRQ(ierr);
  /* Gravity collected terms */
  ierr = DAVecGetArray(MySim->SimDA, MySim->Grav0, &LocalGrav0);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Grav1, &LocalGrav1);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Grav2, &LocalGrav2);CHKERRQ(ierr);
  /* Grab the local data for flow rates */
  ierr = DAVecGetArray(MySim->SimDA, MySim->Qo, &LocalBo);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Qw, &LocalBw);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Qg, &LocalBg);CHKERRQ(ierr);
  /* Grab the local data for volume factors at the cell centers */
  ierr = DAVecGetArray(MySim->SimDA, MySim->Bo, &LocalBo);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Bw, &LocalBw);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Bg, &LocalBg);CHKERRQ(ierr);
  /* Grab the local data for the oil gas solubility */
  ierr = DAVecGetArray(MySim->SimDA, MySim->Rso, &LocalRso);CHKERRQ(ierr);
  /* Grab the distances for the volume */
  ierr = DAVecGetArray(MySim->SimDA, MySim->h1, &Localh1);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->h2, &Localh2);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->h3, &Localh3);CHKERRQ(ierr);
  /* Grab the collected right hand sides */
  ierr = DAVecGetArray(MySim->SimDA, MySim->RHS0, &LocalRHS0);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->RHS1, &LocalRHS1);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->RHS2, &LocalRHS2);CHKERRQ(ierr);




  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DefiantNewtonRaphson3PhAssembleRHS"
extern PetscErrorCode DefiantNewtonRaphson3PhAssembleRHS(BlackOilReservoirSimulation* MySim)
{
  /* Sum up the RHS: gravity terms, capillary pressure terms, flow rates and other terms */

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DefiantNewtonRaphson3PhSolve"
extern PetscErrorCode DefiantNewtonRaphson3PhSolve(BlackOilReservoirSimulation* MySim)
{
  /* Perform the actual solve */

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DefiantNewtonRaphsonComputeLegacyTerms"
extern PetscErrorCode DefiantNewtonRaphsonComputeLegacyTerms(BlackOilReservoirSimulation* MySim)
{
  /* Here we compute PhiSoByBoOld, PhiSwByBwOld, PhiSgByBgPlusRsoSoByBoOld
   * at the end of the iteration so that we can retain the values for next
   * iteration
   */

  PetscFunctionReturn(0);
}

/* SNES updator */
#undef __FUNCT__
#define __FUNCT__ "DefiantNewtonRaphson3PhFormFunction"
extern PetscErrorCode DefiantNewtonRaphson3PhFormFunction(SNES snes,Vec X,Vec F,void *ptr)
{
  /* Set the reservoir simulation to the application context */
  BlackOilReservoirSimulation* MySim = (BlackOilReservoirSimulation*) ptr;
  /* Update the SNES RHS */
  PetscErrorCode ierr;
  PetscInt i, j, k, mx, my, mz, xm, ym, zm, xs, ys, zs;
  /* Flowmask and porosity */
  PetscScalar ***LocalFlowMask, ***LocalPhi;
  /* Grab the pressures and saturations */
  PetscScalar ***LocalPo, ***LocalPw, ***LocalPg;
  PetscScalar ***LocalSo, ***LocalSw, ***LocalSg;
  /* Transmissibility For Oil */
  PetscScalar ***LocalTox1m, ***LocalTox2m, ***LocalTox3m;
  PetscScalar ***LocalTox1p, ***LocalTox2p, ***LocalTox3p;
  /* Transmissibility For Water */
  PetscScalar ***LocalTwx1m, ***LocalTwx2m, ***LocalTwx3m;
  PetscScalar ***LocalTwx1p, ***LocalTwx2p, ***LocalTwx3p;
  /* Transmissibility For Gas */
  PetscScalar ***LocalTgx1m, ***LocalTgx2m, ***LocalTgx3m;
  PetscScalar ***LocalTgx1p, ***LocalTgx2p, ***LocalTgx3p;
  /* Gravity collected terms */
  PetscScalar ***LocalGrav0, ***LocalGrav1, ***LocalGrav2;
  /* Flow rates for oil, water, gas */
  PetscScalar ***LocalQo, ***LocalQw, ***LocalQg;
  /* Volume factors for oil, water and gas at cell centers */
  PetscScalar ***LocalBo, ***LocalBw, ***LocalBg;
  /* Gas oil solubility */
  PetscScalar ***LocalRso;
  PetscScalar ***LocalRsox1p, ***LocalRsox2p, ***LocalRsox3p;
  PetscScalar ***LocalRsox1m, ***LocalRsox2m, ***LocalRsox3m;
  /* dimensions  for volume */
  PetscScalar ***Localh1, ***Localh2, ***Localh3;
  /* Temporary vectors */
  PetscScalar ***LocalPhiSoByBoOld, ***LocalPhiSwByBwOld, ***LocalPhiSgByBgPlusRsoSoByBoOld;
  /* The local vectors */
  PetscScalar ****LocalX, ****LocalF;

  PetscFunctionBegin;
  /* Set MySim to the current simulation */
  MySim = CurrentSimulation;
  /* Get dimensions and extents of the local vectors */
  ierr = DAGetInfo(MySim->SimDA, 0, &mx, &my, &mz, 0, 0, 0, 0, 0, 0, 0);CHKERRQ(ierr);
  ierr = DAGetCorners(MySim->SimDA, &xs, &ys, &zs, &xm, &ym, &zm);CHKERRQ(ierr);
  /* Grab the data for the flow field */
  ierr = DAVecGetArray(MySim->SimDA, MySim->FlowMask, &LocalFlowMask);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Phi, &LocalPhi);CHKERRQ(ierr);
  /* Grab the local data for pressures and saturations */
  ierr = DAVecGetArray(MySim->SimDA, MySim->Po, &LocalPo);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->So, &LocalSo);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Sw, &LocalSw);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Sg, &LocalSg);CHKERRQ(ierr);
  /* Grab the local information for X and F */
  ierr = DAVecGetArrayDOF(MySim->NewtonRaphsonDA, X, &LocalX);CHKERRQ(ierr);
  ierr = DAVecGetArrayDOF(MySim->NewtonRaphsonDA, F, &LocalF);CHKERRQ(ierr);

  /* First off, grab the local values from X */
  for (k = zs; k < zs + zm; k++) {
    for (j = ys; j < ys + ym; j++) {
      for (i = xs; i < xs + xm; i++) {
        if (i == 0 || j == 0 || k == 0 || i == mx - 1 || j == my - 1 || k == mz - 1) { }
        else if (ABS(LocalFlowMask[k][j][i]-FLUID_FLOW) < EPSILON) {
          LocalPo[k][j][i] = LocalX[k][j][i][0];
          LocalSo[k][j][i] = LocalX[k][j][i][1];
          LocalSw[k][j][i] = LocalX[k][j][i][2];
          LocalSg[k][j][i] = 1.0 - (LocalSo[k][j][i] + LocalSw[k][j][i]);
        }
      }
    }
  }

  /* We have to push back the Po, So, Sw and Sg */
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Po, &LocalPo);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->So, &LocalSo);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Sw, &LocalSw);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Sg, &LocalSg);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->Po);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->So);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->Sw);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->Sg);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->Po);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->So);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->Sw);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->Sg);CHKERRQ(ierr);


  /* we have new saturations so we calculate new capillary pressures */
  ierr =  DefiantUpdatePcow(MySim);
  ierr =  DefiantUpdatePcog(MySim);
  /* With new capillary pressures we update pressures */
  ierr =  DefiantUpdatePwFromPcow(MySim);
  ierr =  DefiantUpdatePgFromPcog(MySim);
  /* Now we get the Pw and Pg again */
  ierr = DAVecGetArray(MySim->SimDA, MySim->Pw, &LocalPw);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Pg, &LocalPg);CHKERRQ(ierr);
  /* With new pressures we update all properties */
  ierr = DefiantUpdatePorosity(MySim);
  ierr = DefiantUpdateDensity(MySim);
  ierr = DefiantUpdateVolumeFactors(MySim);
  ierr = DefiantUpdateViscosities(MySim);
  ierr = DefiantUpdateRelativePermeability(MySim);
  /* Now we interpolate at the cell faces */
  ierr = DefiantComputeRhoAndMuAtFaces(MySim);
  ierr = DefiantComputeRelativePermsAtFaces(MySim);
  ierr = DefiantComputeRelativePermsAtFaces(MySim);
  /* And now for transmissibilities */
  ierr = DefiantComputeTransmissibilities(MySim);
  /* With transmissibility we can compute gravity terms */
  ierr = DefiantNewtonRaphson3PhGravity(MySim);
  /* Finally the Wells */
  ierr = DefiantNewtonRaphson3PhHandleWells(MySim);

  /* grab the porosity */
  ierr = DAVecGetArray(MySim->SimDA, MySim->Phi, &LocalPhi);CHKERRQ(ierr);
  /* Grab the local data for pressures */
  ierr = DAVecGetArray(MySim->SimDA, MySim->Pw, &LocalPw);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Pg, &LocalPg);CHKERRQ(ierr);
  /* Grab the local data for the oil transmissibilities */
  ierr = DAVecGetArray(MySim->SimDA, MySim->Tox1m, &LocalTox1m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Tox1p, &LocalTox1p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Tox2m, &LocalTox2m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Tox2p, &LocalTox2p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Tox3m, &LocalTox3m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Tox3p, &LocalTox3p);CHKERRQ(ierr);
  /* Grab the local data for the water transmissibilities */
  ierr = DAVecGetArray(MySim->SimDA, MySim->Twx1m, &LocalTwx1m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Twx1p, &LocalTwx1p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Twx2m, &LocalTwx2m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Twx2p, &LocalTwx2p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Twx3m, &LocalTwx3m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Twx3p, &LocalTwx3p);CHKERRQ(ierr);
  /* Grab the local data for the gas transmissibilities */
  ierr = DAVecGetArray(MySim->SimDA, MySim->Tgx1m, &LocalTgx1m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Tgx1p, &LocalTgx1p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Tgx2m, &LocalTgx2m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Tgx2p, &LocalTgx2p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Tgx3m, &LocalTgx3m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Tgx3p, &LocalTgx3p);CHKERRQ(ierr);
  /* Gravity collected terms */
  ierr = DAVecGetArray(MySim->SimDA, MySim->Grav0, &LocalGrav0);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Grav1, &LocalGrav1);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Grav2, &LocalGrav2);CHKERRQ(ierr);
  /* grab the flow rates collected terms */
  ierr = DAVecGetArray(MySim->SimDA, MySim->Qo, &LocalQo);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Qw, &LocalQw);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Qg, &LocalQg);CHKERRQ(ierr);
  /* Grab the local data for volume factors at the cell centers */
  ierr = DAVecGetArray(MySim->SimDA, MySim->Bo, &LocalBo);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Bw, &LocalBw);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Bg, &LocalBg);CHKERRQ(ierr);
  /* Grab the local data for the oil gas solubility */
  ierr = DAVecGetArray(MySim->SimDA, MySim->Rso, &LocalRso);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Rsox1p, &LocalRsox1p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Rsox2p, &LocalRsox2p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Rsox3p, &LocalRsox3p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Rsox1m, &LocalRsox1m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Rsox2m, &LocalRsox2m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Rsox3m, &LocalRsox3m);CHKERRQ(ierr);
  /* Grab the distances for the volume */
  ierr = DAVecGetArray(MySim->SimDA, MySim->h1, &Localh1);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->h2, &Localh2);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->h3, &Localh3);CHKERRQ(ierr);
  /* Grab the local temporary solutions */
  ierr = DAVecGetArray(MySim->SimDA, MySim->PhiSoByBoOld, &LocalPhiSoByBoOld);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->PhiSwByBwOld, &LocalPhiSwByBwOld);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->PhiSgByBgPlusRsoSoByBoOld, &LocalPhiSgByBgPlusRsoSoByBoOld);CHKERRQ(ierr);
  /* Grab the local information for X and F */
  ierr = DAVecGetArrayDOF(MySim->NewtonRaphsonDA, X, &LocalX);CHKERRQ(ierr);
  ierr = DAVecGetArrayDOF(MySim->NewtonRaphsonDA, F, &LocalF);CHKERRQ(ierr);

  /* we have all the info we need for the residual computation */
  for (k = zs; k < zs + zm; k++) {
    for (j = ys; j < ys + ym; j++) {
      for (i = xs; i < xs + xm; i++) {
        if (i == 0 || j == 0 || k == 0 || i == mx - 1 || j == my - 1 || k == mz - 1) { }
        else if (ABS(LocalFlowMask[k][j][i]-FLUID_FLOW) < EPSILON) {
          LocalF[k][j][i][0] = 1.0/MySim->DeltaTS *Localh1[k][j][i]*Localh2[k][j][i]*Localh3[k][j][i]
              * (LocalPhi[k][j][i]*LocalSo[k][j][i]/LocalBo[k][j][i] - LocalPhiSoByBoOld[k][j][i])
              - LocalTox1p[k][j][i]*(LocalPo[k][j][i+1]-LocalPo[k][j][i]) + LocalTox1m[k][j][i]*(LocalPo[k][j][i]-LocalPo[k][j][i-1])
              - LocalTox2p[k][j][i]*(LocalPo[k][j+1][i]-LocalPo[k][j][i]) + LocalTox2m[k][j][i]*(LocalPo[k][j][i]-LocalPo[k][j-1][i])
              - LocalTox3p[k][j][i]*(LocalPo[k+1][j][i]-LocalPo[k][j][i]) + LocalTox3m[k][j][i]*(LocalPo[k][j][i]-LocalPo[k-1][j][i])
              + LocalGrav0[k][j][i] - LocalQo[k][j][i];

          LocalF[k][j][i][1] = 1.0/MySim->DeltaTS *Localh1[k][j][i]*Localh2[k][j][i]*Localh3[k][j][i]
              * (LocalPhi[k][j][i]*LocalSw[k][j][i]/LocalBw[k][j][i] - LocalPhiSwByBwOld[k][j][i])
              - LocalTwx1p[k][j][i]*(LocalPw[k][j][i+1]-LocalPw[k][j][i]) + LocalTwx1m[k][j][i]*(LocalPw[k][j][i]-LocalPw[k][j][i-1])
              - LocalTwx2p[k][j][i]*(LocalPw[k][j+1][i]-LocalPw[k][j][i]) + LocalTwx2m[k][j][i]*(LocalPw[k][j][i]-LocalPw[k][j-1][i])
              - LocalTwx3p[k][j][i]*(LocalPw[k+1][j][i]-LocalPw[k][j][i]) + LocalTwx3m[k][j][i]*(LocalPw[k][j][i]-LocalPw[k-1][j][i])
              + LocalGrav1[k][j][i] - LocalQw[k][j][i];

          LocalF[k][j][i][2] =  1.0/MySim->DeltaTS *Localh1[k][j][i]*Localh2[k][j][i]*Localh3[k][j][i]
              * (LocalPhi[k][j][i]*(LocalSg[k][j][i]/LocalBg[k][j][i] + LocalRso[k][j][i]*LocalSo[k][j][i]/LocalBo[k][j][i]) - LocalPhiSoByBoOld[k][j][i])
              - LocalTgx1p[k][j][i]*(LocalPg[k][j][i+1]-LocalPg[k][j][i]) + LocalTgx1m[k][j][i]*(LocalPg[k][j][i]-LocalPg[k][j][i-1])
              - LocalTgx2p[k][j][i]*(LocalPg[k][j+1][i]-LocalPg[k][j][i]) + LocalTgx2m[k][j][i]*(LocalPg[k][j][i]-LocalPg[k][j-1][i])
              - LocalTgx3p[k][j][i]*(LocalPg[k+1][j][i]-LocalPg[k][j][i]) + LocalTgx3m[k][j][i]*(LocalPg[k][j][i]-LocalPg[k-1][j][i])
              - LocalRsox1p[k][j][i] * LocalTox1p[k][j][i]*(LocalPo[k][j][i+1]-LocalPo[k][j][i]) + LocalRsox1m[k][j][i] * LocalTox1m[k][j][i]*(LocalPo[k][j][i]-LocalPo[k][j][i-1])
              - LocalRsox2p[k][j][i] * LocalTox2p[k][j][i]*(LocalPo[k][j+1][i]-LocalPo[k][j][i]) + LocalRsox2m[k][j][i] * LocalTox2m[k][j][i]*(LocalPo[k][j][i]-LocalPo[k][j-1][i])
              - LocalRsox3p[k][j][i] * LocalTox3p[k][j][i]*(LocalPo[k+1][j][i]-LocalPo[k][j][i]) + LocalRsox3m[k][j][i] * LocalTox3m[k][j][i]*(LocalPo[k][j][i]-LocalPo[k-1][j][i])
              + LocalGrav0[k][j][i] - LocalQo[k][j][i];
        }
      }
    }
  }

  ierr = DAVecRestoreArray(MySim->NewtonRaphsonDA, X, &LocalX);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->NewtonRaphsonDA, F, &LocalF);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(X);
  ierr = VecAssemblyBegin(F);
  ierr = VecAssemblyEnd(X);
  ierr = VecAssemblyEnd(F);

  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "DefiantNewtonRaphsonGetColoring"
extern PetscErrorCode DefiantNewtonRaphsonGetColoring(BlackOilReservoirSimulation* MySim)
{
  PetscErrorCode ierr;
  PetscInt i, j, k, mx, my, mz, xm, ym, zm, xs, ys, zs;
  PetscScalar v[21];
  MatStencil row, col[21];
  /* Well handling variables */
  PetscInt PerfIDMine, PerfIDOther, WellID;
  PetscInt MyI, MyJ, MyK;
  PetscInt OtherI, OtherJ, OtherK;
  /* The temporary matrix we use to get the coloring */
  Mat TempMat;

  PetscFunctionBegin;
  ierr = DAGetInfo(MySim->SimDA, 0, &mx, &my, &mz, 0, 0, 0, 0, 0, 0, 0);CHKERRQ(ierr);
  ierr = DAGetCorners(MySim->SimDA, &xs, &ys, &zs, &xm, &ym, &zm);CHKERRQ(ierr);
  /* create a temporary matrix to get the coloring */
  ierr = DAGetMatrix(MySim->NewtonRaphsonDA ,MATAIJ,&TempMat);CHKERRQ(ierr);

  /* this is what we will do: add ones where we have a column in the row */
  /* First pass is for all first order finite differences */
  for (k = zs; k < zs + zm; k++) {
    for (j = ys; j < ys + ym; j++) {
      for (i = xs; i < xs + xm; i++) {
        /* For oil pressure, it is related to both Sw's and So' */
        row.i = i;row.j = j;row.k = k;row.c = 0;
        v[0] = 1.0; col[0].i = i; col[0].j = j; col[0].k = k - 1; col[0].c = 0;
        v[1] = 1.0; col[1].i = i; col[1].j = j - 1; col[1].k = k; col[0].c = 0;
        v[2] = 1.0; col[2].i = i - 1; col[2].j = j; col[2].k = k; col[0].c = 0;
        /* moved v[3] down to sum over all other v's */
        v[4] = 1.0; col[4].i = i + 1; col[4].j = j; col[4].k = k; col[0].c = 0;
        v[5] = 1.0; col[5].i = i; col[5].j = j + 1; col[5].k = k; col[0].c = 0;
        v[6] = 1.0; col[6].i = i; col[6].j = j; col[6].k = k + 1; col[0].c = 0;
        v[3] = 1.0;
        col[3].i = row.i; col[3].j = row.j; col[3].k = row.k; col[3].c = 0;
        /* Pressure relating to oil saturation */
        v[7] = 1.0; col[7].i = i; col[7].j = j; col[7].k = k - 1; col[7].c = 1;
        v[8] = 1.0; col[8].i = i; col[8].j = j - 1; col[8].k = k; col[8].c = 1;
        v[9] = 1.0; col[9].i = i - 1; col[9].j = j; col[9].k = k; col[9].c = 1;
        /* moved v[3] down to sum over all other v's */
        v[11] = 1.0; col[11].i = i + 1; col[11].j = j; col[11].k = k; col[11].c = 1;
        v[12] = 1.0; col[12].i = i; col[12].j = j + 1; col[12].k = k; col[12].c = 1;
        v[13] = 1.0; col[13].i = i; col[13].j = j; col[13].k = k + 1; col[13].c = 1;
        v[10] = 1.0;
        col[10].i = row.i; col[10].j = row.j; col[10].k = row.k; col[10].c = 1;
        /* Pressure relating to water saturation */
        v[14] = 1.0; col[14].i = i; col[14].j = j; col[14].k = k - 1; col[14].c = 2;
        v[15] = 1.0; col[15].i = i; col[15].j = j - 1; col[15].k = k; col[15].c = 2;
        v[16] = 1.0; col[16].i = i - 1; col[16].j = j; col[16].k = k; col[16].c = 2;
        /* moved v[3] down to sum over all other v's */
        v[18] = 1.0; col[18].i = i + 1; col[18].j = j; col[18].k = k; col[18].c = 2;
        v[19] = 1.0; col[19].i = i; col[19].j = j + 1; col[19].k = k; col[19].c = 2;
        v[20] = 1.0; col[20].i = i; col[20].j = j; col[20].k = k + 1; col[20].c = 2;
        v[17] = 1.0;
        col[17].i = row.i; col[17].j = row.j; col[17].k = row.k; col[17].c = 2;
        ierr = MatSetValuesStencil(TempMat, 1, &row, 21, col, v, INSERT_VALUES);CHKERRQ(ierr);


        /* For oil Saturation, it is related to both Sw's and So' */
        row.i = i;row.j = j;row.k = k;row.c = 1;
        v[0] = 1.0; col[0].i = i; col[0].j = j; col[0].k = k - 1; col[0].c = 0;
        v[1] = 1.0; col[1].i = i; col[1].j = j - 1; col[1].k = k; col[0].c = 0;
        v[2] = 1.0; col[2].i = i - 1; col[2].j = j; col[2].k = k; col[0].c = 0;
        /* moved v[3] down to sum over all other v's */
        v[4] = 1.0; col[4].i = i + 1; col[4].j = j; col[4].k = k; col[0].c = 0;
        v[5] = 1.0; col[5].i = i; col[5].j = j + 1; col[5].k = k; col[0].c = 0;
        v[6] = 1.0; col[6].i = i; col[6].j = j; col[6].k = k + 1; col[0].c = 0;
        v[3] = 1.0;
        col[3].i = row.i; col[3].j = row.j; col[3].k = row.k; col[3].c = 0;
        /* Pressure relating to oil saturation */
        v[7] = 1.0; col[7].i = i; col[7].j = j; col[7].k = k - 1; col[7].c = 1;
        v[8] = 1.0; col[8].i = i; col[8].j = j - 1; col[8].k = k; col[8].c = 1;
        v[9] = 1.0; col[9].i = i - 1; col[9].j = j; col[9].k = k; col[9].c = 1;
        /* moved v[3] down to sum over all other v's */
        v[11] = 1.0; col[11].i = i + 1; col[11].j = j; col[11].k = k; col[11].c = 1;
        v[12] = 1.0; col[12].i = i; col[12].j = j + 1; col[12].k = k; col[12].c = 1;
        v[13] = 1.0; col[13].i = i; col[13].j = j; col[13].k = k + 1; col[13].c = 1;
        v[10] = 1.0;
        col[10].i = row.i; col[10].j = row.j; col[10].k = row.k; col[10].c = 1;
        /* Pressure relating to water saturation */
        v[14] = 1.0; col[14].i = i; col[14].j = j; col[14].k = k - 1; col[14].c = 2;
        v[15] = 1.0; col[15].i = i; col[15].j = j - 1; col[15].k = k; col[15].c = 2;
        v[16] = 1.0; col[16].i = i - 1; col[16].j = j; col[16].k = k; col[16].c = 2;
        /* moved v[3] down to sum over all other v's */
        v[18] = 1.0; col[18].i = i + 1; col[18].j = j; col[18].k = k; col[18].c = 2;
        v[19] = 1.0; col[19].i = i; col[19].j = j + 1; col[19].k = k; col[19].c = 2;
        v[20] = 1.0; col[20].i = i; col[20].j = j; col[20].k = k + 1; col[20].c = 2;
        v[17] = 1.0;
        col[17].i = row.i; col[17].j = row.j; col[17].k = row.k; col[17].c = 2;
        ierr = MatSetValuesStencil(TempMat, 1, &row, 21, col, v, INSERT_VALUES);CHKERRQ(ierr);

        /* For oil Saturation, it is related to both Sw's and So' */
        row.i = i;row.j = j;row.k = k;row.c = 2;
        v[0] = 1.0; col[0].i = i; col[0].j = j; col[0].k = k - 1; col[0].c = 0;
        v[1] = 1.0; col[1].i = i; col[1].j = j - 1; col[1].k = k; col[0].c = 0;
        v[2] = 1.0; col[2].i = i - 1; col[2].j = j; col[2].k = k; col[0].c = 0;
        /* moved v[3] down to sum over all other v's */
        v[4] = 1.0; col[4].i = i + 1; col[4].j = j; col[4].k = k; col[0].c = 0;
        v[5] = 1.0; col[5].i = i; col[5].j = j + 1; col[5].k = k; col[0].c = 0;
        v[6] = 1.0; col[6].i = i; col[6].j = j; col[6].k = k + 1; col[0].c = 0;
        v[3] = 1.0;
        col[3].i = row.i; col[3].j = row.j; col[3].k = row.k; col[3].c = 0;
        /* Pressure relating to oil saturation */
        v[7] = 1.0; col[7].i = i; col[7].j = j; col[7].k = k - 1; col[7].c = 1;
        v[8] = 1.0; col[8].i = i; col[8].j = j - 1; col[8].k = k; col[8].c = 1;
        v[9] = 1.0; col[9].i = i - 1; col[9].j = j; col[9].k = k; col[9].c = 1;
        /* moved v[3] down to sum over all other v's */
        v[11] = 1.0; col[11].i = i + 1; col[11].j = j; col[11].k = k; col[11].c = 1;
        v[12] = 1.0; col[12].i = i; col[12].j = j + 1; col[12].k = k; col[12].c = 1;
        v[13] = 1.0; col[13].i = i; col[13].j = j; col[13].k = k + 1; col[13].c = 1;
        v[10] = 1.0;
        col[10].i = row.i; col[10].j = row.j; col[10].k = row.k; col[10].c = 1;
        /* Pressure relating to water saturation */
        v[14] = 1.0; col[14].i = i; col[14].j = j; col[14].k = k - 1; col[14].c = 2;
        v[15] = 1.0; col[15].i = i; col[15].j = j - 1; col[15].k = k; col[15].c = 2;
        v[16] = 1.0; col[16].i = i - 1; col[16].j = j; col[16].k = k; col[16].c = 2;
        /* moved v[3] down to sum over all other v's */
        v[18] = 1.0; col[18].i = i + 1; col[18].j = j; col[18].k = k; col[18].c = 2;
        v[19] = 1.0; col[19].i = i; col[19].j = j + 1; col[19].k = k; col[19].c = 2;
        v[20] = 1.0; col[20].i = i; col[20].j = j; col[20].k = k + 1; col[20].c = 2;
        v[17] = 1.0;
        col[17].i = row.i; col[17].j = row.j; col[17].k = row.k; col[17].c = 2;
        ierr = MatSetValuesStencil(TempMat, 1, &row, 21, col, v, INSERT_VALUES);CHKERRQ(ierr);
      }
    }
  }

  /* Now handle the wells */
  for (WellID = 0; WellID < MySim->NumberOfWells; WellID++) {
    for (PerfIDMine = 0; PerfIDMine
        < (MySim->Wells[WellID]).NumberOfPerforations; PerfIDMine++) {
      if (MySim->Wells[WellID].Perforations[PerfIDMine].IsActive == PETSC_TRUE){
        /* grab my location */
        MyI = (MySim->Wells[WellID]).Perforations[PerfIDMine].I;
        MyJ = (MySim->Wells[WellID]).Perforations[PerfIDMine].J;
        MyK = (MySim->Wells[WellID]).Perforations[PerfIDMine].K;

        if (MySim->Wells[WellID].Perforations[PerfIDMine].Constraint == BHP_CONSTRAINT) {
          if ((MySim->Wells[WellID]).NumberOfPerforations > 1) {
            for (PerfIDOther = 0; PerfIDOther < (MySim->Wells[WellID]).NumberOfPerforations; PerfIDOther++) {
              OtherI = (MySim->Wells[WellID]).Perforations[PerfIDOther].I;
              OtherJ = (MySim->Wells[WellID]).Perforations[PerfIDOther].J;
              OtherK = (MySim->Wells[WellID]).Perforations[PerfIDOther].K;
              if (PerfIDOther != PerfIDMine && MySim->Wells[WellID].Perforations[PerfIDOther].IsActive == PETSC_TRUE) {
                /* we correlate pressure at my location with all pressures and saturations of
                 * perforations belonging to the same well */
                row.i = MyI;row.j = MyJ;row.k = MyK;row.c = 0;
                v[0] = 1.0;
                col[0].i = OtherI; col[0].j = OtherJ; col[0].k = OtherK; col[0].c = 0;
                /* Pressure relating to oil saturation */
                v[1] = 1.0;
                col[1].i = OtherI; col[1].j = OtherJ; col[1].k = OtherK; col[1].c = 1;
                v[2] = 1.0;
                col[2].i = OtherI; col[2].j = OtherJ; col[2].k = OtherK; col[2].c = 2;
                ierr = MatSetValuesStencil(TempMat, 1, &row, 3, col, v, INSERT_VALUES);CHKERRQ(ierr);

                /* we correlate oil saturation at my location with all pressures and saturations of
                 * perforations belonging to the same well */
                row.i = MyI;row.j = MyJ;row.k = MyK;row.c = 1;
                v[0] = 1.0;
                col[0].i = OtherI; col[0].j = OtherJ; col[0].k = OtherK; col[0].c = 0;
                /* Pressure relating to oil saturation */
                v[1] = 1.0;
                col[1].i = OtherI; col[1].j = OtherJ; col[1].k = OtherK; col[1].c = 1;
                v[2] = 1.0;
                col[2].i = OtherI; col[2].j = OtherJ; col[2].k = OtherK; col[2].c = 2;
                ierr = MatSetValuesStencil(TempMat, 1, &row, 3, col, v, INSERT_VALUES);CHKERRQ(ierr);

                /* we correlate water saturation at my location with all pressures and saturations of
                 * perforations belonging to the same well */
                row.i = MyI;row.j = MyJ;row.k = MyK;row.c = 2;
                v[0] = 1.0;
                col[0].i = OtherI; col[0].j = OtherJ; col[0].k = OtherK; col[0].c = 0;
                /* Pressure relating to oil saturation */
                v[1] = 1.0;
                col[1].i = OtherI; col[1].j = OtherJ; col[1].k = OtherK; col[1].c = 1;
                v[2] = 1.0;
                col[2].i = OtherI; col[2].j = OtherJ; col[2].k = OtherK; col[2].c = 2;
                ierr = MatSetValuesStencil(TempMat, 1, &row, 3, col, v, INSERT_VALUES);CHKERRQ(ierr);

              }
            }
          }
        }
      }
    }
  }

  ierr = MatAssemblyBegin(TempMat, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(TempMat, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  /* Check the actual matrix */


  PetscFunctionReturn(0);
}

/* Iterating routines */
#undef __FUNCT__
#define __FUNCT__ "DefiantNewtonRaphson3PhIterate"
extern PetscErrorCode DefiantNewtonRaphson3PhIterate(BlackOilReservoirSimulation* MySim)
{

  PetscFunctionReturn(0);
}
