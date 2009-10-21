/*
 * IMPES2Ph.c
 *
 *  Created on: Sep 6, 2009
 *      Author: yye00
 */

#include "Defiant.h"


#undef __FUNCT__
#define __FUNCT__ "Defiant2PhGravity"
extern PetscErrorCode Defiant2PhGravity(BlackOilReservoirSimulation* MySim) {
  PetscErrorCode ierr;
  PetscInt i, j, k, mx, my, mz, xm, ym, zm, xs, ys, zs;
  PetscScalar ***LocalFlowMask;
  /* Depth */
  PetscScalar ***Localx3;
  /* Area * Permeability divided by height */
  PetscScalar ***LocalAKHx1m, ***LocalAKHx2m, ***LocalAKHx3m;
  PetscScalar ***LocalAKHx1p, ***LocalAKHx2p, ***LocalAKHx3p;
  /* Density By Viscosity at the faces */
  PetscScalar ***LocalRhoByMuox1m, ***LocalRhoByMuox1p, ***LocalRhoByMuwx1m,
      ***LocalRhoByMuwx1p;
  PetscScalar ***LocalRhoByMuox2m, ***LocalRhoByMuox2p, ***LocalRhoByMuwx2m,
      ***LocalRhoByMuwx2p;
  PetscScalar ***LocalRhoByMuox3m, ***LocalRhoByMuox3p, ***LocalRhoByMuwx3m,
      ***LocalRhoByMuwx3p;
  /* Relative Permeabilities at the faces */
  PetscScalar ***LocalRelPermox1m, ***LocalRelPermox1p, ***LocalRelPermox2m,
      ***LocalRelPermox2p, ***LocalRelPermox3m, ***LocalRelPermox3p;
  PetscScalar ***LocalRelPermwx1m, ***LocalRelPermwx1p, ***LocalRelPermwx2m,
      ***LocalRelPermwx2p, ***LocalRelPermwx3m, ***LocalRelPermwx3p;
  /* Gravity collected term */
  PetscScalar ***LocalGravity;

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
  ierr = DAVecGetArray(MySim->SimDA, MySim->RhoByMuox2m, &LocalRhoByMuox2m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->RhoByMuox2p, &LocalRhoByMuox2p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->RhoByMuwx2m, &LocalRhoByMuwx2m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->RhoByMuwx2p, &LocalRhoByMuwx2p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->RhoByMuox3m, &LocalRhoByMuox3m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->RhoByMuox3p, &LocalRhoByMuox3p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->RhoByMuwx3m, &LocalRhoByMuwx3m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->RhoByMuwx3p, &LocalRhoByMuwx3p);CHKERRQ(ierr);
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
  /* Gravity collected term */
  ierr = DAVecGetArray(MySim->SimDA, MySim->Gravity, &LocalGravity);CHKERRQ(ierr);

  for (k = zs; k < zs + zm; k++) {
    for (j = ys; j < ys + ym; j++) {
      for (i = xs; i < xs + xm; i++) {
        if (i == 0 || j == 0 || k == 0 || i == mx - 1 || j == my - 1 || k == mz
            - 1) {
        } else if (ABS(LocalFlowMask[k][j][i]-FLUID_FLOW) < EPSILON) {
          LocalGravity[k][j][i] = MySim->GravAcc * (

               LocalAKHx1m[k][j][i] * (
                 LocalRelPermox1m[k][j][i] * LocalRhoByMuox1m[k][j][i]
               + LocalRelPermwx1m[k][j][i] * LocalRhoByMuwx1m[k][j][i])
                 * (Localx3[k][j][i] - Localx3[k][j][i - 1])

             + LocalAKHx2m[k][j][i] * (
                 LocalRelPermox2m[k][j][i] * LocalRhoByMuox2m[k][j][i]
               + LocalRelPermwx2m[k][j][i] * LocalRhoByMuwx2m[k][j][i])
                 * (Localx3[k][j][i] - Localx3[k][j - 1][i])

             + LocalAKHx3m[k][j][i] * (
                 LocalRelPermox3m[k][j][i] * LocalRhoByMuox3m[k][j][i]
               + LocalRelPermwx3m[k][j][i] * LocalRhoByMuwx3m[k][j][i])
                 * (Localx3[k][j][i] - Localx3[k - 1][j][i])

             - LocalAKHx1p[k][j][i] * (
                 LocalRelPermox1p[k][j][i] * LocalRhoByMuox1p[k][j][i]
               + LocalRelPermwx1p[k][j][i] * LocalRhoByMuwx1p[k][j][i])
                 * (Localx3[k][j][i + 1] - Localx3[k][j][i])

             - LocalAKHx2p[k][j][i] * (
                 LocalRelPermox2p[k][j][i] * LocalRhoByMuox2p[k][j][i]
               + LocalRelPermwx2p[k][j][i] * LocalRhoByMuwx2p[k][j][i])
                 * (Localx3[k][j + 1][i] - Localx3[k][j][i])

             - LocalAKHx3p[k][j][i] * (
                 LocalRelPermox3p[k][j][i] * LocalRhoByMuox3p[k][j][i]
               + LocalRelPermwx3p[k][j][i] * LocalRhoByMuwx3p[k][j][i])
                  * (Localx3[k + 1][j][i] - Localx3[k][j][i])  );
        }
      }
    }
  }

  /* Restore the new arrays to their rightful place */
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Gravity, &LocalGravity);CHKERRQ(ierr);

  /* Begin Assembly for vectors */
  ierr = VecAssemblyBegin(MySim->Gravity);CHKERRQ(ierr);
  /* And end Assembly */
  ierr = VecAssemblyEnd(MySim->Gravity);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "Defiant2PhCapillaryPressure"
extern PetscErrorCode Defiant2PhCapillaryPressure(
    BlackOilReservoirSimulation* MySim) {
  PetscErrorCode ierr;
  PetscInt i, j, k, mx, my, mz, xm, ym, zm, xs, ys, zs;
  PetscScalar ***LocalFlowMask;
  /* Transmissibility For Water*/
  PetscScalar ***LocalTwx1m, ***LocalTwx2m, ***LocalTwx3m;
  PetscScalar ***LocalTwx1p, ***LocalTwx2p, ***LocalTwx3p;
  /* Capillary pressure */
  PetscScalar ***LocalPcow;
  /* Capillary pressure collected term */
  PetscScalar ***LocalCapillaryPressure;

  PetscFunctionBegin;

  /* Get dimensions and extents of the local vectors */
  ierr = DAGetInfo(MySim->SimDA, 0, &mx, &my, &mz, 0, 0, 0, 0, 0, 0, 0);CHKERRQ(ierr);
  ierr = DAGetCorners(MySim->SimDA, &xs, &ys, &zs, &xm, &ym, &zm);CHKERRQ(ierr);
  /* Grab the data for the flow field */
  ierr = DAVecGetArray(MySim->SimDA, MySim->FlowMask, &LocalFlowMask);CHKERRQ(ierr);
  /* Grab the local data for the transmissibilities */
  ierr = DAVecGetArray(MySim->SimDA, MySim->Twx1m, &LocalTwx1m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Twx1p, &LocalTwx1p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Twx2m, &LocalTwx2m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Twx2p, &LocalTwx2p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Twx3m, &LocalTwx3m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Twx3p, &LocalTwx3p);CHKERRQ(ierr);
  /* Grab the local data for the capillary pressure */
  ierr = DAVecGetArray(MySim->SimDA, MySim->Pcow, &LocalPcow);CHKERRQ(ierr);
  /* Grab the local data for the collected capillary pressure */
  ierr = DAVecGetArray(MySim->SimDA, MySim->CapillaryPressure, &LocalCapillaryPressure);CHKERRQ(ierr);

  for (k = zs; k < zs + zm; k++) {
    for (j = ys; j < ys + ym; j++) {
      for (i = xs; i < xs + xm; i++) {
        if (i == 0 || j == 0 || k == 0 || i == mx - 1 || j == my - 1 || k == mz
            - 1) {
        } else if (ABS(LocalFlowMask[k][j][i]-FLUID_FLOW) < EPSILON) {
          LocalCapillaryPressure[k][j][i] = LocalTwx1m[k][j][i]
              * (LocalPcow[k][j][i] - LocalPcow[k][j][i - 1])
              + LocalTwx2m[k][j][i] * (LocalPcow[k][j][i]
                  - LocalPcow[k][j - 1][i]) + LocalTwx3m[k][j][i]
              * (LocalPcow[k][j][i] - LocalPcow[k - 1][j][i])
              - LocalTwx1p[k][j][i] * (LocalPcow[k][j][i + 1]
                  - LocalPcow[k][j][i]) - LocalTwx2p[k][j][i] * (LocalPcow[k][j
              + 1][i] - LocalPcow[k][j][i]) - LocalTwx3p[k][j][i]
              * (LocalPcow[k + 1][j][i] - LocalPcow[k][j][i]);
        }
      }
    }
  }

  /* Restore the new arrays to their rightful place */
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->CapillaryPressure,
      &LocalCapillaryPressure);CHKERRQ(ierr);

  /* Begin Assembly for vectors */
  ierr = VecAssemblyBegin(MySim->CapillaryPressure);CHKERRQ(ierr);
  /* And end Assembly */
  ierr = VecAssemblyEnd(MySim->CapillaryPressure);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DefiantIMPES2PhAssembleMatrix"
extern PetscErrorCode DefiantIMPES2PhAssembleMatrix(BlackOilReservoirSimulation* MySim) {
  PetscErrorCode ierr;
  PetscInt i, j, k, mx, my, mz, xm, ym, zm, xs, ys, zs;
  PetscScalar v[7];
  MatStencil row, col[7];
  /* Transmissibility For Oil*/
  PetscScalar ***LocalTox1p, ***LocalTox2p, ***LocalTox3p;
  PetscScalar ***LocalTox1m, ***LocalTox2m, ***LocalTox3m;
  /* Transmissibility For Water*/
  PetscScalar ***LocalTwx1p, ***LocalTwx2p, ***LocalTwx3p;
  PetscScalar ***LocalTwx1m, ***LocalTwx2m, ***LocalTwx3m;

  PetscFunctionBegin;
  /* Get dimensions and extents of the local vectors */
  ierr = DAGetInfo(MySim->SimDA, 0, &mx, &my, &mz, 0, 0, 0, 0, 0, 0, 0);CHKERRQ(ierr);
  ierr = DAGetCorners(MySim->SimDA, &xs, &ys, &zs, &xm, &ym, &zm);CHKERRQ(ierr);

  /* Grab the local data for the transmissibilities */
  ierr = DAVecGetArray(MySim->SimDA, MySim->Tox1m, &LocalTox1m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Tox1p, &LocalTox1p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Twx1m, &LocalTwx1m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Twx1p, &LocalTwx1p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Tox2m, &LocalTox2m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Tox2p, &LocalTox2p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Twx2m, &LocalTwx2m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Twx2p, &LocalTwx2p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Tox3m, &LocalTox3m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Tox3p, &LocalTox3p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Twx3m, &LocalTwx3m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Twx3p, &LocalTwx3p);CHKERRQ(ierr);
  for (k = zs; k < zs + zm; k++) {
    for (j = ys; j < ys + ym; j++) {
      for (i = xs; i < xs + xm; i++) {
        row.i = i;
        row.j = j;
        row.k = k;
        if (i == 0 || j == 0 || k == 0 || i == mx - 1 || j == my - 1 || k == mz
            - 1) {
          v[0] = 0.0;
          ierr = MatSetValuesStencil(MySim->A, 1, &row, 1, &row, v,
              INSERT_VALUES);
          CHKERRQ(ierr);
        } else {
          v[0] = -(LocalTox3m[k][j][i] + LocalTwx3m[k][j][i]);
          col[0].i = i;
          col[0].j = j;
          col[0].k = k - 1;
          v[1] = -(LocalTox2m[k][j][i] + LocalTwx2m[k][j][i]);
          col[1].i = i;
          col[1].j = j - 1;
          col[1].k = k;
          v[2] = -(LocalTox1m[k][j][i] + LocalTwx1m[k][j][i]);
          col[2].i = i - 1;
          col[2].j = j;
          col[2].k = k;
          /* moved v[3] down to sum over all other v's */
          v[4] = -(LocalTox1p[k][j][i] + LocalTwx1p[k][j][i]);
          col[4].i = i + 1;
          col[4].j = j;
          col[4].k = k;
          v[5] = -(LocalTox2p[k][j][i] + LocalTwx2p[k][j][i]);
          col[5].i = i;
          col[5].j = j + 1;
          col[5].k = k;
          v[6] = -(LocalTox3p[k][j][i] + LocalTwx3p[k][j][i]);
          col[6].i = i;
          col[6].j = j;
          col[6].k = k + 1;
          /* v[3] is the sum of all other column values */
          v[3] = -1.0*(v[0] + v[1] + v[2] + v[4] + v[5] + v[6]);
          col[3].i = row.i;
          col[3].j = row.j;
          col[3].k = row.k;
          ierr = MatSetValuesStencil(MySim->A, 1, &row, 7, col, v,
              INSERT_VALUES);
          CHKERRQ(ierr);
        }
      }
    }
  }
  ierr = MatAssemblyBegin(MySim->A, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(MySim->A, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DefiantIMPES2PhAssembleRHS"
extern PetscErrorCode DefiantIMPES2PhAssembleRHS(BlackOilReservoirSimulation* MySim)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;

  /* Add gravity terms to RHS */
  ierr = VecAXPY(MySim->RHS, 1.0, MySim->Gravity);CHKERRQ(ierr);
  /* Add capillary pressure to RHS */
  ierr = VecAXPY(MySim->RHS, 1.0, MySim->CapillaryPressure);CHKERRQ(ierr);
  /* Assemble the RHS vector */
  ierr = VecAssemblyBegin(MySim->RHS);CHKERRQ(ierr);
  /* And end Assembly */
  ierr = VecAssemblyEnd(MySim->RHS);CHKERRQ(ierr);
  PetscFunctionReturn(0);

}


#undef __FUNCT__
#define __FUNCT__ "DefiantIMPES2PhHandleWells"
extern PetscErrorCode DefiantIMPES2PhHandleWells(BlackOilReservoirSimulation* MySim)
{
  PetscErrorCode ierr;
  PetscInt mx, my, mz, xm, ym, zm, xs, ys, zs;
  PetscScalar v;
  MatStencil row, col;
  /* Well handling variables */
  PetscInt PerfIDMine, PerfIDOther, WellID;
  PetscInt MyI, MyJ, MyK;
  PetscInt OtherI, OtherJ, OtherK;
  PetscScalar apo, apw;
  PetscScalar qo, qw;
  /* Local values for variables*/
  PetscScalar ***LocalRHS;
  PetscScalar ***LocalQo, ***LocalQw;

  PetscFunctionBegin;
  /* Get dimensions and extents of the local vectors */
  ierr = DAGetInfo(MySim->SimDA, 0, &mx, &my, &mz, 0, 0, 0, 0, 0, 0, 0);CHKERRQ(ierr);
  ierr = DAGetCorners(MySim->SimDA, &xs, &ys, &zs, &xm, &ym, &zm);CHKERRQ(ierr);
  /* Grab the local RHS */
  ierr = DAVecGetArray(MySim->SimDA, MySim->RHS, &LocalRHS);CHKERRQ(ierr);
  /* Grab the local flow rate */
  ierr = DAVecGetArray(MySim->SimDA, MySim->Qw, &LocalQw);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Qo, &LocalQo);CHKERRQ(ierr);

  for (WellID = 0; WellID < MySim->NumberOfWells; WellID++) {
    for (PerfIDMine = 0; PerfIDMine
        < (MySim->Wells[WellID]).NumberOfPerforations; PerfIDMine++) {
      MyI = (MySim->Wells[WellID]).Perforations[PerfIDMine].I;
      MyJ = (MySim->Wells[WellID]).Perforations[PerfIDMine].J;
      MyK = (MySim->Wells[WellID]).Perforations[PerfIDMine].K;

      if (MyI >= xs && MyI < xs+xm && MyJ >= ys && MyJ < ys+ym && MyK >= zs && MyK < zs+zm)
      {
        if (MySim->Wells[WellID].Perforations[PerfIDMine].IsActive == PETSC_TRUE ){
          /* initialize temp variables for flow rates */
          qo = 0.0; qw = 0.0;
          if (MySim->Wells[WellID].Perforations[PerfIDMine].Constraint == FLOW_RATE_CONSTRAINT){

            if (MySim->Wells[WellID].Perforations[PerfIDMine].WellType == WATER_INJECTOR )
            {
              qw = (MySim->Wells[WellID]).Perforations[PerfIDMine].Qw;
              LocalRHS[MyK][MyJ][MyI] = LocalRHS[MyK][MyJ][MyI] + qw;
              LocalQw[MyK][MyJ][MyI] = qw;
            }
            else if (MySim->Wells[WellID].Perforations[PerfIDMine].WellType == OIL_INJECTOR )
            {
              qo = (MySim->Wells[WellID]).Perforations[PerfIDMine].Qo;
              LocalRHS[MyK][MyJ][MyI] = LocalRHS[MyK][MyJ][MyI] + qo;
              LocalQo[MyK][MyJ][MyI] = qo;
            }
            else if (MySim->Wells[WellID].Perforations[PerfIDMine].WellType == WATER_PRODUCER )
            {
              qw = (MySim->Wells[WellID]).Perforations[PerfIDMine].Qw;
              LocalRHS[MyK][MyJ][MyI] = LocalRHS[MyK][MyJ][MyI] + qw;
              LocalQw[MyK][MyJ][MyI] = qw;
            }
            else if (MySim->Wells[WellID].Perforations[PerfIDMine].WellType == OIL_PRODUCER )
            {
              qo = (MySim->Wells[WellID]).Perforations[PerfIDMine].Qo;
              LocalRHS[MyK][MyJ][MyI] = LocalRHS[MyK][MyJ][MyI] + qo;
              LocalQo[MyK][MyJ][MyI] = qo;
            }
          } else if (MySim->Wells[WellID].Perforations[PerfIDMine].Constraint == BHP_CONSTRAINT) {
            /* set the row for the insert */
            row.i = MyI; row.j = MyJ; row.k = MyK;
            /* initialize temp variables for flow rates */
            qo = 0.0; qw = 0.0;
            /* initialize temp variable for Ap */
            apo = 0.0, apw = 0.0;
            if ((MySim->Wells[WellID]).NumberOfPerforations == 1) {
              /* for all other perforations belonging to the same well do this */
              qo  = MySim->Wells[WellID].Perforations[PerfIDMine].WellIndex
                  * MySim->Wells[WellID].Perforations[PerfIDMine].Kro / MySim->Wells[WellID].Perforations[PerfIDMine].Muo
                  * ((MySim->Wells[WellID]).Perforations[PerfIDMine].BHPo
                  - MySim->GravAcc * MySim->Wells[WellID].Perforations[PerfIDMine].Rhoo
                  * (MySim->Wells[WellID].Perforations[PerfIDMine].zbh
                  - MySim->Wells[WellID].Perforations[PerfIDMine].x3)) * MySim->Wells[WellID].Perforations[PerfIDMine].h1
                  * MySim->Wells[WellID].Perforations[PerfIDMine].h2 * MySim->Wells[WellID].Perforations[PerfIDMine].h3
                  / MySim->Wells[WellID].Perforations[PerfIDMine].Bo;
              qw  = MySim->Wells[WellID].Perforations[PerfIDMine].WellIndex
                  * MySim->Wells[WellID].Perforations[PerfIDMine].Krw / MySim->Wells[WellID].Perforations[PerfIDMine].Muw
                  * (MySim->Wells[WellID].Perforations[PerfIDMine].BHPw
                  + MySim->Wells[WellID].Perforations[PerfIDMine].Pcow - MySim->GravAcc
                  * MySim->Wells[WellID].Perforations[PerfIDMine].Rhow
                  * (MySim->Wells[WellID].Perforations[PerfIDMine].zbh
                  - MySim->Wells[WellID].Perforations[PerfIDMine].x3)) * MySim->Wells[WellID].Perforations[PerfIDMine].h1
                  * MySim->Wells[WellID].Perforations[PerfIDMine].h2 * MySim->Wells[WellID].Perforations[PerfIDMine].h3
                  / MySim->Wells[WellID].Perforations[PerfIDMine].Bw;
              apo = MySim->Wells[WellID].Perforations[PerfIDMine].WellIndex
                  * MySim->Wells[WellID].Perforations[PerfIDMine].Kro / MySim->Wells[WellID].Perforations[PerfIDMine].Muo
                  * MySim->Wells[WellID].Perforations[PerfIDMine].h1 * MySim->Wells[WellID].Perforations[PerfIDMine].h2 * MySim->Wells[WellID].Perforations[PerfIDMine].h3
                  / MySim->Wells[WellID].Perforations[PerfIDMine].Bo;
              apw = MySim->Wells[WellID].Perforations[PerfIDMine].WellIndex
                  * MySim->Wells[WellID].Perforations[PerfIDMine].Krw / MySim->Wells[WellID].Perforations[PerfIDMine].Muw
                  * MySim->Wells[WellID].Perforations[PerfIDMine].h1 * MySim->Wells[WellID].Perforations[PerfIDMine].h2 * MySim->Wells[WellID].Perforations[PerfIDMine].h3
                  / MySim->Wells[WellID].Perforations[PerfIDMine].Bw;

              /* Since one perforation, the total well is the value here */
              MySim->Wells[WellID].TotalQo = qo;
              MySim->Wells[WellID].TotalQw = qw;
              /* Store the flow rates and relevant data somewhere useful */
              MySim->Wells[WellID].Perforations[PerfIDMine].Apo = apo;
              MySim->Wells[WellID].Perforations[PerfIDMine].Apw = apw;
              LocalQo[MyK][MyJ][MyI] = qo;
              LocalQw[MyK][MyJ][MyI] = qw;
              /* Add values to the Ap of the matrix */
              row.i = MyI; row.j = MyJ; row.k = MyK;
              v = apo+apw;
              ierr = MatSetValuesStencil(MySim->A, 1, &row, 1, &row, &v, ADD_VALUES);CHKERRQ(ierr);
              /* Add values to the RHS for the same column */
              LocalRHS[MyK][MyJ][MyI] = LocalRHS[MyK][MyJ][MyI] + qo + qw;
            } else if ((MySim->Wells[WellID]).NumberOfPerforations > 1) {
              /* Set qo and qw to totals, we will be subtracting perfs from other cells from those values */
              qo = MySim->Wells[WellID].TotalQo;
              qw = MySim->Wells[WellID].TotalQw;
              for (PerfIDOther = 0; PerfIDOther
                  < (MySim->Wells[WellID]).NumberOfPerforations; PerfIDOther++) {
                if (PerfIDOther != PerfIDMine && MySim->Wells[WellID].Perforations[PerfIDOther].IsActive == PETSC_TRUE) {
                  OtherI = (MySim->Wells[WellID]).Perforations[PerfIDOther].I;
                  OtherJ = (MySim->Wells[WellID]).Perforations[PerfIDOther].J;
                  OtherK = (MySim->Wells[WellID]).Perforations[PerfIDOther].K;
                  /* for all other perforations belonging to the same well do this */
                  qo  = qo - (MySim->Wells[WellID]).Perforations[PerfIDOther].WellIndex
                      * MySim->Wells[WellID].Perforations[PerfIDOther].Kro / MySim->Wells[WellID].Perforations[PerfIDOther].Muo
                      * ((MySim->Wells[WellID]).Perforations[PerfIDOther].BHPo
                      - MySim->GravAcc * MySim->Wells[WellID].Perforations[PerfIDOther].Rhoo
                      * ((MySim->Wells[WellID]).Perforations[PerfIDOther].zbh
                      - MySim->Wells[WellID].Perforations[PerfIDOther].x3)) * MySim->Wells[WellID].Perforations[PerfIDOther].h1
                      * MySim->Wells[WellID].Perforations[PerfIDOther].h2 * MySim->Wells[WellID].Perforations[PerfIDOther].h3
                      / MySim->Wells[WellID].Perforations[PerfIDOther].Bo;
                  qw  = qw - (MySim->Wells[WellID]).Perforations[PerfIDOther].WellIndex
                      * MySim->Wells[WellID].Perforations[PerfIDOther].Krw / MySim->Wells[WellID].Perforations[PerfIDOther].Muw
                      * ((MySim->Wells[WellID]).Perforations[PerfIDOther].BHPw
                      + MySim->Wells[WellID].Perforations[PerfIDOther].Pcow - MySim->GravAcc
                      * MySim->Wells[WellID].Perforations[PerfIDOther].Rhow
                      * ((MySim->Wells[WellID]).Perforations[PerfIDOther].zbh
                      - MySim->Wells[WellID].Perforations[PerfIDOther].x3)) * MySim->Wells[WellID].Perforations[PerfIDOther].h1
                      * MySim->Wells[WellID].Perforations[PerfIDOther].h2 * MySim->Wells[WellID].Perforations[PerfIDOther].h3
                      / MySim->Wells[WellID].Perforations[PerfIDOther].Bw;
                  apo = -(MySim->Wells[WellID]).Perforations[PerfIDOther].WellIndex
                      * MySim->Wells[WellID].Perforations[PerfIDOther].Kro / MySim->Wells[WellID].Perforations[PerfIDOther].Muo
                      * MySim->Wells[WellID].Perforations[PerfIDOther].h1 * MySim->Wells[WellID].Perforations[PerfIDOther].h2 * MySim->Wells[WellID].Perforations[PerfIDOther].h3
                      / MySim->Wells[WellID].Perforations[PerfIDOther].Bo;
                  apw = - (MySim->Wells[WellID]).Perforations[PerfIDOther].WellIndex
                      * MySim->Wells[WellID].Perforations[PerfIDOther].Krw / MySim->Wells[WellID].Perforations[PerfIDOther].Muw
                      * MySim->Wells[WellID].Perforations[PerfIDOther].h1 * MySim->Wells[WellID].Perforations[PerfIDOther].h2 * MySim->Wells[WellID].Perforations[PerfIDOther].h3
                      / MySim->Wells[WellID].Perforations[PerfIDOther].Bw;

                  /* Add values to the Ap of the matrix */
                  /* we are at the same row but at a different column */
                  row.i = MyI; row.j = MyJ; row.k = MyK;
                  col.i = OtherI; col.j =OtherJ; col.k = OtherK;
                  v = apo+apw;
                  ierr = MatSetValuesStencil(MySim->A, 1, &row, 1, &col, &v, ADD_VALUES);CHKERRQ(ierr);
                }
              }
              /* now that I have the final qo and qw add them to the RHS */
              /* Add values to the RHS for the same column */
              LocalRHS[MyK][MyJ][MyI] = LocalRHS[MyK][MyJ][MyI] + qo + qw;
              LocalQo[MyK][MyJ][MyI] = qo;
              LocalQw[MyK][MyJ][MyI] = qw;
            }
          }
        }
      }
    }
  }

  /* Restore the RHS to it's place */
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->RHS,&LocalRHS);CHKERRQ(ierr);
  /* Restore the Qo and Qw */
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Qo,&LocalQo);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Qw,&LocalQw);CHKERRQ(ierr);
  /* Begin assembly the vectors */
  ierr = VecAssemblyBegin(MySim->RHS);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->Qo);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->Qw);CHKERRQ(ierr);
  /* End assembly the vectors */
  ierr = VecAssemblyEnd(MySim->RHS);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->Qo);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->Qw);CHKERRQ(ierr);

  /* Finalize by assembling the matrix */
  ierr = MatAssemblyBegin(MySim->A, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(MySim->A, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "DefiantIMPES2PhHandleWellsForSaturation"
extern PetscErrorCode DefiantIMPES2PhHandleWellsForSaturation(BlackOilReservoirSimulation* MySim)
{
  PetscErrorCode ierr;
  PetscInt mx, my, mz, xm, ym, zm, xs, ys, zs;
  /* Well handling variables */
  PetscInt PerfIDMine, PerfIDOther, WellID;
  PetscScalar apo, apw;
  PetscScalar qo, qw;
  PetscInt MyI, MyJ, MyK;
  PetscInt OtherI, OtherJ, OtherK;
  /* Local values for variables*/
  PetscScalar ***LocalRHS;
  PetscScalar ***LocalQo, ***LocalQw;

  PetscFunctionBegin;
  /* Get dimensions and extents of the local vectors */
  ierr = DAGetInfo(MySim->SimDA, 0, &mx, &my, &mz, 0, 0, 0, 0, 0, 0, 0);CHKERRQ(ierr);
  ierr = DAGetCorners(MySim->SimDA, &xs, &ys, &zs, &xm, &ym, &zm);CHKERRQ(ierr);
  /* Grab the local RHS */
  ierr = DAVecGetArray(MySim->SimDA, MySim->RHS, &LocalRHS);CHKERRQ(ierr);
  /* Grab the local flow rate */
  ierr = DAVecGetArray(MySim->SimDA, MySim->Qw, &LocalQw);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Qo, &LocalQo);CHKERRQ(ierr);

  for (WellID = 0; WellID < MySim->NumberOfWells; WellID++) {
    for (PerfIDMine = 0; PerfIDMine
        < (MySim->Wells[WellID]).NumberOfPerforations; PerfIDMine++) {
      if (MySim->Wells[WellID].Perforations[PerfIDMine].IsActive == PETSC_TRUE){
        if (MySim->Wells[WellID].Perforations[PerfIDMine].Constraint == FLOW_RATE_CONSTRAINT){
          /* grab my index */
          MyI = (MySim->Wells[WellID]).Perforations[PerfIDMine].I;
          MyJ = (MySim->Wells[WellID]).Perforations[PerfIDMine].J;
          MyK = (MySim->Wells[WellID]).Perforations[PerfIDMine].K;

          if (MySim->Wells[WellID].Perforations[PerfIDMine].WellType == WATER_INJECTOR )
          {
            qw = (MySim->Wells[WellID]).Perforations[PerfIDMine].Qw;
            LocalQw[MyK][MyJ][MyI] = qw;
          }
          else if (MySim->Wells[WellID].Perforations[PerfIDMine].WellType == OIL_INJECTOR )
          {
            qo = (MySim->Wells[WellID]).Perforations[PerfIDMine].Qo;
            LocalQw[MyK][MyJ][MyI] = qo;
          }
          if (MySim->Wells[WellID].Perforations[PerfIDMine].WellType == WATER_PRODUCER )
          {
            qw = (MySim->Wells[WellID]).Perforations[PerfIDMine].Qw;
            LocalQw[MyK][MyJ][MyI] = qw;
          }
          else if (MySim->Wells[WellID].Perforations[PerfIDMine].WellType == OIL_PRODUCER )
          {
            qo = (MySim->Wells[WellID]).Perforations[PerfIDMine].Qo;
            LocalQw[MyK][MyJ][MyI] = qo;
          }
        } else if (MySim->Wells[WellID].Perforations[PerfIDMine].Constraint == BHP_CONSTRAINT) {
          /* initialize temp variables for flow rates */
          qo = 0.0; qw = 0.0;
          /* initialize temp variable for Ap */
          apo = 0.0, apw = 0.0;
          if ((MySim->Wells[WellID]).NumberOfPerforations == 1) {
            /* for all other perforations belonging to the same well do this */
            MyI = (MySim->Wells[WellID]).Perforations[PerfIDMine].I;
            MyJ = (MySim->Wells[WellID]).Perforations[PerfIDMine].J;
            MyK = (MySim->Wells[WellID]).Perforations[PerfIDMine].K;
            qo  = MySim->Wells[WellID].Perforations[PerfIDMine].WellIndex
                * MySim->Wells[WellID].Perforations[PerfIDMine].Kro / MySim->Wells[WellID].Perforations[PerfIDMine].Muo
                * (MySim->Wells[WellID].Perforations[PerfIDMine].BHPo
                - MySim->Wells[WellID].Perforations[PerfIDMine].Po
                - MySim->GravAcc * MySim->Wells[WellID].Perforations[PerfIDMine].Rhoo
                * (MySim->Wells[WellID].Perforations[PerfIDMine].zbh
                - MySim->Wells[WellID].Perforations[PerfIDMine].x3)) * MySim->Wells[WellID].Perforations[PerfIDMine].h1
                * MySim->Wells[WellID].Perforations[PerfIDMine].h2 * MySim->Wells[WellID].Perforations[PerfIDMine].h3
                / MySim->Wells[WellID].Perforations[PerfIDMine].Bo;
            qw  = MySim->Wells[WellID].Perforations[PerfIDMine].WellIndex
                * MySim->Wells[WellID].Perforations[PerfIDMine].Krw / MySim->Wells[WellID].Perforations[PerfIDMine].Muw
                * (MySim->Wells[WellID].Perforations[PerfIDMine].BHPw
                - MySim->Wells[WellID].Perforations[PerfIDMine].Pw
                + MySim->Wells[WellID].Perforations[PerfIDMine].Pcow - MySim->GravAcc
                * MySim->Wells[WellID].Perforations[PerfIDMine].Rhow
                * (MySim->Wells[WellID].Perforations[PerfIDMine].zbh
                - MySim->Wells[WellID].Perforations[PerfIDMine].x3)) * MySim->Wells[WellID].Perforations[PerfIDMine].h1
                * MySim->Wells[WellID].Perforations[PerfIDMine].h2 * MySim->Wells[WellID].Perforations[PerfIDMine].h3
                / MySim->Wells[WellID].Perforations[PerfIDMine].Bw;
            /* Store the flow rates and relevant data somewhere useful */
            LocalQo[MyK][MyJ][MyI] = qo;
            LocalQw[MyK][MyJ][MyI] = qw;
          } else if ((MySim->Wells[WellID]).NumberOfPerforations > 1) {
            /* Set qo and qw to totals, we will be subtracting perfs from other cells from those values */
            qo = MySim->Wells[WellID].TotalQo;
            qw = MySim->Wells[WellID].TotalQw;
            for (PerfIDOther = 0; PerfIDOther
                < (MySim->Wells[WellID]).NumberOfPerforations; PerfIDOther++) {
              OtherI = (MySim->Wells[WellID]).Perforations[PerfIDOther].I;
              OtherJ = (MySim->Wells[WellID]).Perforations[PerfIDOther].J;
              OtherK = (MySim->Wells[WellID]).Perforations[PerfIDOther].K;
              if (PerfIDOther != PerfIDMine && MySim->Wells[WellID].Perforations[PerfIDOther].IsActive == PETSC_TRUE) {
                /* for all other perforations belonging to the same well do this */
                qo  = qo - MySim->Wells[WellID].Perforations[PerfIDOther].WellIndex
                    * MySim->Wells[WellID].Perforations[PerfIDOther].Kro / MySim->Wells[WellID].Perforations[PerfIDOther].Muo
                    * (MySim->Wells[WellID].Perforations[PerfIDOther].BHPo
                    - MySim->Wells[WellID].Perforations[PerfIDOther].Po
                    - MySim->GravAcc * MySim->Wells[WellID].Perforations[PerfIDOther].Rhoo
                    * (MySim->Wells[WellID].Perforations[PerfIDOther].zbh
                    - MySim->Wells[WellID].Perforations[PerfIDOther].x3)) * MySim->Wells[WellID].Perforations[PerfIDOther].h1
                    * MySim->Wells[WellID].Perforations[PerfIDOther].h2 * MySim->Wells[WellID].Perforations[PerfIDOther].h3
                    / MySim->Wells[WellID].Perforations[PerfIDOther].Bo;
                qw  = qw - MySim->Wells[WellID].Perforations[PerfIDOther].WellIndex
                    * MySim->Wells[WellID].Perforations[PerfIDOther].Krw / MySim->Wells[WellID].Perforations[PerfIDOther].Muw
                    * (MySim->Wells[WellID].Perforations[PerfIDOther].BHPw
                    - MySim->Wells[WellID].Perforations[PerfIDOther].Pw
                    + MySim->Wells[WellID].Perforations[PerfIDOther].Pcow - MySim->GravAcc
                    * MySim->Wells[WellID].Perforations[PerfIDOther].Rhow
                    * (MySim->Wells[WellID].Perforations[PerfIDOther].zbh
                    - MySim->Wells[WellID].Perforations[PerfIDOther].x3)) * MySim->Wells[WellID].Perforations[PerfIDOther].h1
                    * MySim->Wells[WellID].Perforations[PerfIDOther].h2 * MySim->Wells[WellID].Perforations[PerfIDOther].h3
                    / MySim->Wells[WellID].Perforations[PerfIDOther].Bw;
              }
            }
            LocalQo[MyK][MyJ][MyI] = qo;
            LocalQw[MyK][MyJ][MyI] = qw;
          }
        }
      }
    }
  }

  /* Restore the Qo and Qw */
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Qo,&LocalQo);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Qw,&LocalQw);CHKERRQ(ierr);
  /* Begin assembly the vectors */
  ierr = VecAssemblyBegin(MySim->Qo);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->Qw);CHKERRQ(ierr);
  /* End assembly the vectors */
  ierr = VecAssemblyEnd(MySim->Qo);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->Qw);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "DefiantIMPES2PhUpdateSaturations"
extern PetscErrorCode DefiantIMPES2PhUpdateSaturations(BlackOilReservoirSimulation* MySim)
{
  PetscErrorCode ierr;
  PetscInt i, j, k, mx, my, mz, xm, ym, zm, xs, ys, zs;
  PetscScalar ***LocalFlowMask;
  /* Some geometry */
  PetscScalar ***Localh1, ***Localh2, ***Localh3;
  PetscScalar ***Localx3;
  /* Transmissibility For Water*/
  PetscScalar ***LocalTwx1p, ***LocalTwx2p, ***LocalTwx3p;
  PetscScalar ***LocalTwx1m, ***LocalTwx2m, ***LocalTwx3m;
  /* Densities at the cell faces */
  PetscScalar ***LocalRhowx1p, ***LocalRhowx2p, ***LocalRhowx3p;
  PetscScalar ***LocalRhowx1m, ***LocalRhowx2m, ***LocalRhowx3m;
  /* Volume factors at cell center */
  PetscScalar ***LocalBw;
  /* Pressure and porosity */
  PetscScalar ***LocalPw;
  PetscScalar ***LocalPhi;
  /* Flow rates */
  PetscScalar ***LocalQw;
  /* Saturations */
  PetscScalar ***LocalSo, ***LocalSw;
  /* Temporary vector */
  Vec vecLocalPw;

  /* Temporary IMPES variables */
  PetscInt TempInt;
  PetscReal TempSw;
  Vec  TempDSDT;
  PetscScalar ***LocalTempDSDT;

  PetscFunctionBegin;
  /* First create the temporary vector TempDSDT */
  ierr = VecDuplicate(MySim->Sw, &TempDSDT);CHKERRQ(ierr);
  /* Get dimensions and extents of the local vectors */
  ierr = DAGetInfo(MySim->SimDA, 0, &mx, &my, &mz, 0, 0, 0, 0, 0, 0, 0);CHKERRQ(ierr);
  ierr = DAGetCorners(MySim->SimDA, &xs, &ys, &zs, &xm, &ym, &zm);CHKERRQ(ierr);
  /* Grab the local data for the temporary Sw */
  ierr = DAVecGetArray(MySim->SimDA, TempDSDT, &LocalTempDSDT);CHKERRQ(ierr);
  /* Grab the data for the flow field */
  ierr = DAVecGetArray(MySim->SimDA, MySim->FlowMask, &LocalFlowMask);CHKERRQ(ierr);
  /* Grab local geometry */
  ierr = DAVecGetArray(MySim->SimDA, MySim->h1, &Localh1);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->h2, &Localh2);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->h3, &Localh3);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->x3, &Localx3);CHKERRQ(ierr);
  /* Grab transmissibilities */
  ierr = DAVecGetArray(MySim->SimDA, MySim->Twx1p, &LocalTwx1p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Twx2p, &LocalTwx2p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Twx3p, &LocalTwx3p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Twx1m, &LocalTwx1m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Twx2m, &LocalTwx2m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Twx3m, &LocalTwx3m);CHKERRQ(ierr);
  /* Grab the densities */
  ierr = DAVecGetArray(MySim->SimDA, MySim->Rhowx1p, &LocalRhowx1p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Rhowx2p, &LocalRhowx2p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Rhowx3p, &LocalRhowx3p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Rhowx1m, &LocalRhowx1m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Rhowx2m, &LocalRhowx2m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Rhowx3m, &LocalRhowx3m);CHKERRQ(ierr);
  /* Grab the volume factor */
  ierr = DAVecGetArray(MySim->SimDA, MySim->Bw, &LocalBw);CHKERRQ(ierr);
  /* Grab the local data for water pressure */
  ierr = DAGetLocalVector(MySim->SimDA, &vecLocalPw);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGlobalToLocalBegin(MySim->SimDA,MySim->Pw,INSERT_VALUES,vecLocalPw);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGlobalToLocalEnd(MySim->SimDA,MySim->Pw,INSERT_VALUES,vecLocalPw);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAVecGetArray(MySim->SimDA, vecLocalPw, &LocalPw);CHKERRQ(ierr);
  /* Get the porosity */
  ierr = DAVecGetArray(MySim->SimDA, MySim->Phi, &LocalPhi);CHKERRQ(ierr);
  /* Grab the local flow rate */
  ierr = DAVecGetArray(MySim->SimDA, MySim->Qw, &LocalQw);CHKERRQ(ierr);
  /* Grab the local data for Saturations */
  ierr = DAVecGetArray(MySim->SimDA, MySim->So, &LocalSo);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Sw, &LocalSw);CHKERRQ(ierr);

  for (k = zs; k < zs + zm; k++) {
    for (j = ys; j < ys + ym; j++) {
      for (i = xs; i < xs + xm; i++) {
        if (i == 0 || j == 0 || k == 0 || i == mx - 1 || j == my - 1 || k == mz - 1) {}
        else if (ABS(LocalFlowMask[k][j][i]-FLUID_FLOW) < EPSILON) {
          TempSw = LocalSw[k][j][i]
               + MySim->DeltaTS*LocalBw[k][j][i]/LocalPhi[k][j][i]
               / (Localh1[k][j][i]*Localh2[k][j][i]*Localh3[k][j][i])
               * (LocalTwx1p[k][j][i]*(LocalPw[k][j][i+1]-LocalPw[k][j][i])-LocalTwx1m[k][j][i]*(LocalPw[k][j][i]-LocalPw[k][j][i-1])
               +  LocalTwx2p[k][j][i]*(LocalPw[k][j+1][i]-LocalPw[k][j][i])-LocalTwx2m[k][j][i]*(LocalPw[k][j][i]-LocalPw[k][j-1][i])
               +  LocalTwx3p[k][j][i]*(LocalPw[k+1][j][i]-LocalPw[k][j][i])-LocalTwx3m[k][j][i]*(LocalPw[k][j][i]-LocalPw[k-1][j][i])
               -  LocalTwx1p[k][j][i]*MySim->GravAcc*LocalRhowx1p[k][j][i]*(Localx3[k][j][i+1]-Localx3[k][j][i])
               +  LocalTwx1m[k][j][i]*MySim->GravAcc*LocalRhowx1m[k][j][i]*(Localx3[k][j][i]-Localx3[k][j][i-1])
               -  LocalTwx2p[k][j][i]*MySim->GravAcc*LocalRhowx2p[k][j][i]*(Localx3[k][j+1][i]-Localx3[k][j][i])
               +  LocalTwx2m[k][j][i]*MySim->GravAcc*LocalRhowx2m[k][j][i]*(Localx3[k][j][i]-Localx3[k][j-1][i])
               -  LocalTwx3p[k][j][i]*MySim->GravAcc*LocalRhowx3p[k][j][i]*(Localx3[k+1][j][i]-Localx3[k][j][i])
               +  LocalTwx3m[k][j][i]*MySim->GravAcc*LocalRhowx3m[k][j][i]*(Localx3[k][j][i]-Localx3[k-1][j][i])
               +  LocalQw[k][j][i] );
          LocalTempDSDT[k][j][i] = ABS(TempSw - LocalSw[k][j][i])/MySim->DeltaTS;
        }
      }
    }
  }

  /* Restore the temp dsdt vector from local values */
  ierr = DAVecRestoreArray(MySim->SimDA, TempDSDT, &LocalTempDSDT);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(TempDSDT);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(TempDSDT);CHKERRQ(ierr);
  ierr = VecMax(TempDSDT, &TempInt, &MySim->DSDTmax);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "\nMaximum DSDTmax is:%f \n", MySim->DSDTmax);CHKERRQ(ierr);
  /* Compute the desired Delta T and set it to the simulation DeltaT */
  /* Check if we are doing adaptive time-step IMPES */
  if(MySim->AdaptiveTimeStep==PETSC_TRUE)
    MySim->DeltaTS = MySim->DSmax/MySim->DSDTmax;

  if (ABS(MySim->DSmax) < EPSILON) MySim->DeltaTS = MySim->DeltaTP;

  /* Now perform the proper update */
  for (k = zs; k < zs + zm; k++) {
    for (j = ys; j < ys + ym; j++) {
      for (i = xs; i < xs + xm; i++) {
        if (i == 0 || j == 0 || k == 0 || i == mx - 1 || j == my - 1 || k == mz - 1) {}
        else if (ABS(LocalFlowMask[k][j][i]-FLUID_FLOW) < EPSILON) {
          LocalSw[k][j][i] = LocalSw[k][j][i]
               + MySim->DeltaTS*LocalBw[k][j][i]/LocalPhi[k][j][i]
               / (Localh1[k][j][i]*Localh2[k][j][i]*Localh3[k][j][i])
               * (LocalTwx1p[k][j][i]*(LocalPw[k][j][i+1]-LocalPw[k][j][i])-LocalTwx1m[k][j][i]*(LocalPw[k][j][i]-LocalPw[k][j][i-1])
               +  LocalTwx2p[k][j][i]*(LocalPw[k][j+1][i]-LocalPw[k][j][i])-LocalTwx2m[k][j][i]*(LocalPw[k][j][i]-LocalPw[k][j-1][i])
               +  LocalTwx3p[k][j][i]*(LocalPw[k+1][j][i]-LocalPw[k][j][i])-LocalTwx3m[k][j][i]*(LocalPw[k][j][i]-LocalPw[k-1][j][i])
               -  LocalTwx1p[k][j][i]*MySim->GravAcc*LocalRhowx1p[k][j][i]*(Localx3[k][j][i+1]-Localx3[k][j][i])
               +  LocalTwx1m[k][j][i]*MySim->GravAcc*LocalRhowx1m[k][j][i]*(Localx3[k][j][i]-Localx3[k][j][i-1])
               -  LocalTwx2p[k][j][i]*MySim->GravAcc*LocalRhowx2p[k][j][i]*(Localx3[k][j+1][i]-Localx3[k][j][i])
               +  LocalTwx2m[k][j][i]*MySim->GravAcc*LocalRhowx2m[k][j][i]*(Localx3[k][j][i]-Localx3[k][j-1][i])
               -  LocalTwx3p[k][j][i]*MySim->GravAcc*LocalRhowx3p[k][j][i]*(Localx3[k+1][j][i]-Localx3[k][j][i])
               +  LocalTwx3m[k][j][i]*MySim->GravAcc*LocalRhowx3m[k][j][i]*(Localx3[k][j][i]-Localx3[k-1][j][i])
               +  LocalQw[k][j][i] );
          if (LocalSw[k][j][i] < MySim->Swc) LocalSw[k][j][i] = MySim->Swc;
          else if (LocalSw[k][j][i] > 1.0 - MySim->Sor) LocalSw[k][j][i] = 1.0 - MySim->Sor;
          LocalSo[k][j][i] = 1.0 - LocalSw[k][j][i];
        }
      }
    }
  }

  /* Restore the new arrays to their rightful place */
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->So, &LocalSo);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Sw, &LocalSw);CHKERRQ(ierr);

  /* Begin Assembly for vectors */
  ierr = VecAssemblyBegin(MySim->So);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->Sw);CHKERRQ(ierr);

  /* And end Assembly */
  ierr = VecAssemblyEnd(MySim->So);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->Sw);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DefiantClearMatrixRHS"
PetscErrorCode DefiantClearMatrixRHS(BlackOilReservoirSimulation* MySim)
{
  PetscErrorCode ierr;
  ierr = VecZeroEntries(MySim->RHS);CHKERRQ(ierr);CHKMEMQ;
  ierr = MatZeroEntries(MySim->A);CHKERRQ(ierr);CHKMEMQ;
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "DefiantIMPES2PhSolve"
extern PetscErrorCode DefiantIMPES2PhSolve(BlackOilReservoirSimulation* MySim)
{
  PetscErrorCode ierr;
  PetscReal norm;
  Vec TempVec;

  PetscFunctionBegin;

  /* Set the current simulation to MySIM */
  CurrentSimulation = MySim;
  /* duplicate Temp with X */
  ierr = VecDuplicate(MySim->Po, &TempVec);CHKERRQ(ierr);CHKMEMQ;
  /* Set the KSP */
  ierr = KSPSetOperators(MySim->ksp, MySim->A, MySim->A, DIFFERENT_NONZERO_PATTERN);
  ierr = KSPGetPC(MySim->ksp,&MySim->pc);CHKERRQ(ierr);
  ierr = PCSetType(MySim->pc,PCJACOBI);CHKERRQ(ierr);
  ierr = KSPSetTolerances(MySim->ksp,1.e-7,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
  ierr = KSPSetFromOptions(MySim->ksp);CHKERRQ(ierr);
  ierr = KSPSolve(MySim->ksp,MySim->RHS,MySim->Po);CHKERRQ(ierr);
  /* solver info */
  //ierr = KSPView(MySim->ksp,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  /* residual info */
  ierr = MatMult(MySim->A, MySim->Po, TempVec);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecAXPY(TempVec, -1.0, MySim->RHS);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecNorm(TempVec, NORM_2, &norm);CHKERRQ(ierr);CHKMEMQ;
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Residual L-2 norm is: %G\n", norm);CHKERRQ(ierr);CHKMEMQ;

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DefiantIMPES2PhIterate"
extern PetscErrorCode DefiantIMPES2PhIterate(BlackOilReservoirSimulation* MySim)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;

  /* Now we start */
  MySim->CurrentTimeP = MySim->StartTime;
  MySim->CurrentTimeS = MySim->StartTime;

  /* Before we start iterating sync the location of the wells */
  ierr = DefiantBlackOilSyncPerfOwners(MySim);CHKERRQ(ierr);CHKMEMQ;

  /* compute some geometric and static factors*/
  ierr = DefiantComputeKAByHAtFaces(MySim);CHKERRQ(ierr);CHKMEMQ;
  /* update the relative permeability and capillary pressure, saturation based dependent variables */
  ierr = DefiantUpdateRelativePermeability(MySim);CHKERRQ(ierr);CHKMEMQ;
  ierr = DefiantUpdatePcow(MySim);CHKERRQ(ierr);CHKMEMQ;

  while(MySim->CurrentTimeP < MySim->EndTime)
  {
    /* we already have the values at cell faces, we interpolate at faces */
    ierr = DefiantComputeRhoAndMuAtFaces(MySim);CHKERRQ(ierr);CHKMEMQ;
    ierr = DefiantComputeRelativePermsAtFaces(MySim);CHKERRQ(ierr);CHKMEMQ;

#if DEFIANT_DEBUG
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Output AKHX1m is: ");
    ierr = VecView(MySim->AKHx1m ,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Output AKHX1p is: ");
    ierr = VecView(MySim->AKHx1p ,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Output RelPermox1m is: ");
    ierr = VecView(MySim->RelPermox1m ,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Output RelPermox1p is: ");
    ierr = VecView(MySim->RelPermox1p ,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Output RelPermwx1m is: ");
    ierr = VecView(MySim->RelPermwx1m ,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Output RelPermwx1p is: ");
    ierr = VecView(MySim->RelPermwx1p ,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
#endif

    ierr = DefiantComputeVolumeFactorsAtFaces(MySim);CHKERRQ(ierr);CHKMEMQ;

    /* we have everything at the faces so we compute transmissbility */
    ierr = DefiantComputeTransmissibilities(MySim);CHKERRQ(ierr);CHKMEMQ;

#if DEFIANT_DEBUG
    ierr = PetscPrintf(PETSC_COMM_WORLD,"RelPermox1m is: ");
    ierr = VecView(MySim->RelPermox1m,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"RelPermox1p is: ");
    ierr = VecView(MySim->RelPermox1p,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"RelPermwx1m is: ");
    ierr = VecView(MySim->RelPermwx1m,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"RelPermwx1p is: ");
    ierr = VecView(MySim->RelPermwx1p,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;

    ierr = PetscPrintf(PETSC_COMM_WORLD,"RelPermox2m is: ");
    ierr = VecView(MySim->RelPermox2m,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"RelPermox2p is: ");
    ierr = VecView(MySim->RelPermox2p,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"RelPermwx2m is: ");
    ierr = VecView(MySim->RelPermwx2m,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"RelPermwx2p is: ");
    ierr = VecView(MySim->RelPermwx2p,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;

    ierr = PetscPrintf(PETSC_COMM_WORLD,"RelPermox3m is: ");
    ierr = VecView(MySim->RelPermox3m,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"RelPermox3p is: ");
    ierr = VecView(MySim->RelPermox3p,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"RelPermwx3m is: ");
    ierr = VecView(MySim->RelPermwx3m,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"RelPermwx3p is: ");
    ierr = VecView(MySim->RelPermwx3p,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;

    ierr = PetscPrintf(PETSC_COMM_WORLD,"Tox1m is: ");
    ierr = VecView(MySim->Tox1m,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Tox1p is: ");
    ierr = VecView(MySim->Tox1p,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Twx1m is: ");
    ierr = VecView(MySim->Twx1m,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Twx1p is: ");
    ierr = VecView(MySim->Twx1p,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;

    ierr = PetscPrintf(PETSC_COMM_WORLD,"Tox2m is: ");
    ierr = VecView(MySim->Tox2m,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Tox2p is: ");
    ierr = VecView(MySim->Tox2p,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Twx2m is: ");
    ierr = VecView(MySim->Twx2m,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Twx2p is: ");
    ierr = VecView(MySim->Twx2p,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;

    ierr = PetscPrintf(PETSC_COMM_WORLD,"Tox3m is: ");
    ierr = VecView(MySim->Tox3m,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Tox3p is: ");
    ierr = VecView(MySim->Tox3p,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Twx3m is: ");
    ierr = VecView(MySim->Twx3m,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Twx3p is: ");
    ierr = VecView(MySim->Twx3p,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
#endif

    /* IMPES 2Ph functions */
    ierr = Defiant2PhGravity(MySim);CHKERRQ(ierr);CHKMEMQ;
    ierr = Defiant2PhCapillaryPressure(MySim);CHKERRQ(ierr);CHKMEMQ;
    /* we clear the values in the matrix and RHS before we assemble */
    ierr = DefiantClearMatrixRHS(MySim);CHKERRQ(ierr);CHKMEMQ;
    ierr = DefiantIMPES2PhAssembleMatrix(MySim);CHKERRQ(ierr);CHKMEMQ;
    /* assemble RHS */
    ierr = DefiantIMPES2PhAssembleRHS(MySim);CHKERRQ(ierr);CHKMEMQ;
    /* Sync the well information */
    ierr = DefiantBlackOilComputePerfIndicesSyncPerfs(MySim);CHKERRQ(ierr);CHKMEMQ;
    /* handling the wells is done AFTER we have the matrix and the RHS */
    ierr = DefiantIMPES2PhHandleWells(MySim);CHKERRQ(ierr);CHKMEMQ;

#if DEFIANT_DEBUG
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Qo is: ");
    ierr = VecView(MySim->Qo,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Qw is: ");
    ierr = VecView(MySim->Qw,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
    /* Diagnostic output */
    ierr = MatView(MySim->A,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
    ierr = VecView(MySim->RHS,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
#endif

    /* we now are ready to solve */
    ierr = DefiantIMPES2PhSolve(MySim);CHKERRQ(ierr);CHKMEMQ;
    ierr = DefiantUpdatePwFromPcow(MySim);CHKERRQ(ierr);CHKMEMQ;

#if DEFIANT_DEBUG
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Output pressure is: ");
    ierr = VecView(MySim->Po,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
#endif

    /* With the new pressure, we update sats then capillary pressures then pressures
     * Then update relative perms, porosity, density, volume factors and viscosities */
    /* update pressure dependent properties */
    ierr = DefiantUpdatePorosity(MySim);CHKERRQ(ierr);CHKMEMQ;
    ierr = DefiantUpdateDensity(MySim);CHKERRQ(ierr);CHKMEMQ;
    ierr = DefiantUpdateVolumeFactors(MySim);CHKERRQ(ierr);CHKMEMQ;
    ierr = DefiantUpdateViscosities(MySim);CHKERRQ(ierr);CHKMEMQ;

    /* interpolate pressure dependent properties */
    ierr = DefiantComputeRhoAndMuAtFaces(MySim);CHKERRQ(ierr);CHKMEMQ;
    ierr = DefiantComputeVolumeFactorsAtFaces(MySim);CHKERRQ(ierr);CHKMEMQ;

    /* Now w eare ready to update the saturations using improved IMPES */
    MySim->CurrentTimeS = MySim->CurrentTimeP;

    while(MySim->CurrentTimeS <= MySim->CurrentTimeP + MySim->DeltaTP)
    {
      ierr = DefiantIMPES2PhUpdateSaturations(MySim);CHKERRQ(ierr);CHKMEMQ;
      /* update saturation based properties */
      ierr = DefiantUpdateRelativePermeability(MySim);CHKERRQ(ierr);CHKMEMQ;
      ierr = DefiantUpdatePcow(MySim);CHKERRQ(ierr);CHKMEMQ;
      /* update pw now that we have a new pcow */
      ierr = DefiantUpdatePwFromPcow(MySim);CHKERRQ(ierr);CHKMEMQ;
      /* with new pcow we have new pw, so new properties */
      ierr = DefiantUpdatePorosity(MySim);CHKERRQ(ierr);CHKMEMQ;
      ierr = DefiantUpdateDensity(MySim);CHKERRQ(ierr);CHKMEMQ;
      ierr = DefiantUpdateVolumeFactors(MySim);CHKERRQ(ierr);CHKMEMQ;
      ierr = DefiantUpdateViscosities(MySim);CHKERRQ(ierr);CHKMEMQ;
      /* update those new properties at the faces with interpolation */
      ierr = DefiantComputeRhoAndMuAtFaces(MySim);CHKERRQ(ierr);CHKMEMQ;
      ierr = DefiantComputeRelativePermsAtFaces(MySim);CHKERRQ(ierr);CHKMEMQ;
      ierr = DefiantComputeVolumeFactorsAtFaces(MySim);CHKERRQ(ierr);CHKMEMQ;

      /* Transmissibility calculation routines */
      ierr = DefiantComputeTransmissibilities(MySim);CHKERRQ(ierr);CHKMEMQ;

      /* Sync the well information */
      ierr = DefiantBlackOilComputePerfIndicesSyncPerfs(MySim);CHKERRQ(ierr);CHKMEMQ;
      /* Finally update the well information */
      ierr = DefiantIMPES2PhHandleWellsForSaturation(MySim);CHKERRQ(ierr);CHKMEMQ;

#if DEFIANT_DEBUG
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Output saturations is: ");
      ierr = VecView(MySim->Sw,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Pcow is: ");
      ierr = VecView(MySim->Pcow,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Pw is: ");
      ierr = VecView(MySim->Pw,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Krw is: ");
      ierr = VecView(MySim->Krw,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
      ierr = PetscPrintf(PETSC_COMM_WORLD,"RelPermwx1p is: ");
      ierr = VecView(MySim->RelPermwx1p,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
      ierr = PetscPrintf(PETSC_COMM_WORLD,"RelPermwx1m is: ");
      ierr = VecView(MySim->RelPermwx1m,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
      ierr = PetscPrintf(PETSC_COMM_WORLD,"RelPermox1p is: ");
      ierr = VecView(MySim->RelPermox1p,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
      ierr = PetscPrintf(PETSC_COMM_WORLD,"RelPermox1m is: ");
      ierr = VecView(MySim->RelPermox1m,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
#endif

      MySim->CurrentTimeS = MySim->CurrentTimeS + MySim->DeltaTS;
      ierr = PetscPrintf(PETSC_COMM_WORLD,"\n Current DeltaTS is: %f ", MySim->DeltaTS);
      ierr = PetscPrintf(PETSC_COMM_WORLD,"\n Current CurrentTimeS is: %f ", MySim->CurrentTimeS);
      ierr = PetscPrintf(PETSC_COMM_WORLD,"\n Current DeltaTP is: %f ", MySim->DeltaTP);
      ierr = PetscPrintf(PETSC_COMM_WORLD,"\n Current CurrentTimeP is: %f ", MySim->CurrentTimeP);

#if DEFIANT_DEBUG
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Sw is: ");
      ierr = VecView(MySim->Sw,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
#endif

    }

    MySim->CurrentTimeP = MySim->CurrentTimeP + MySim->DeltaTP;

    ierr = PetscPrintf(PETSC_COMM_WORLD,"\n Current time is: %f ", MySim->CurrentTimeP);
    if (MySim->CurrentTimeP == 1500){
      ierr = PetscPrintf(PETSC_COMM_WORLD,"FINAL FINAL Output saturations is: ");
      ierr = VecView(MySim->Sw,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
      break;
    }
  }

  PetscFunctionReturn(0);
}


/**********************************************
 *
 * MultiGrid Functions go here
 *
 *********************************************/
#undef __FUNCT__
#define __FUNCT__ "DefiantIMPES2PhDMMGComputeRHS"
extern PetscErrorCode DefiantIMPES2PhDMMGComputeRHS(DMMG dmmg, Vec b);

#undef __FUNCT__
#define __FUNCT__ "DefiantIMPES2PhDMMGComputeRHS"
extern PetscErrorCode DefiantIMPES2PhDMMGComputeMatrix(DMMG dmmg, Mat jac, Mat B);

#undef __FUNCT__
#define __FUNCT__ "DefiantIMPES2PhDMMGSolve"
extern PetscErrorCode DefiantIMPES2PhDMMGSolve(BlackOilReservoirSimulation* MySim);

#undef __FUNCT__
#define __FUNCT__ "DefiantIMPES2PhDMMGIterate"
extern PetscErrorCode DefiantIMPES2PhDMMGIterate(BlackOilReservoirSimulation* MySim)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;

  /* Now we start */
  MySim->CurrentTimeP = MySim->StartTime;
  MySim->CurrentTimeS = MySim->StartTime;

  /* compute some geometric and static factors*/
  ierr = DefiantComputeKAByHAtFaces(MySim);CHKERRQ(ierr);CHKMEMQ;
  /* update the relative permeability and capillary pressure, saturation based dependent variables */
  ierr = DefiantUpdateRelativePermeability(MySim);CHKERRQ(ierr);CHKMEMQ;
  ierr = DefiantUpdatePcow(MySim);CHKERRQ(ierr);CHKMEMQ;

  while(MySim->CurrentTimeP < MySim->EndTime)
  {
    /* we already have the values at cell faces, we interpolate at faces */
    ierr = DefiantComputeRhoAndMuAtFaces(MySim);CHKERRQ(ierr);CHKMEMQ;
    ierr = DefiantComputeRelativePermsAtFaces(MySim);CHKERRQ(ierr);CHKMEMQ;
    ierr = DefiantComputeVolumeFactorsAtFaces(MySim);CHKERRQ(ierr);CHKMEMQ;

    /* we have everything at the faces so we compute transmissbility */
    ierr = DefiantComputeTransmissibilities(MySim);CHKERRQ(ierr);CHKMEMQ;

    /* IMPES 2Ph functions */
    ierr = Defiant2PhGravity(MySim);CHKERRQ(ierr);CHKMEMQ;
    ierr = Defiant2PhCapillaryPressure(MySim);CHKERRQ(ierr);CHKMEMQ;
    /* we clear the values in the matrix and RHS before we assemble */
    ierr = DefiantClearMatrixRHS(MySim);CHKERRQ(ierr);CHKMEMQ;

    /* we now are ready to solve with DMMG*/
    ierr = DefiantIMPES2PhDMMGSolve(MySim);CHKERRQ(ierr);CHKMEMQ;
    ierr = DefiantUpdatePwFromPcow(MySim);CHKERRQ(ierr);CHKMEMQ;

    ierr = PetscPrintf(PETSC_COMM_WORLD,"Output pressure is: ");
    ierr = VecView(MySim->Po,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;

    /* With the new pressure, we update sats then capillary pressures then pressures
     * Then update relative perms, porosity, density, volume factors and viscosities */
    /* update pressure dependent properties */
    ierr = DefiantUpdatePorosity(MySim);CHKERRQ(ierr);CHKMEMQ;
    ierr = DefiantUpdateDensity(MySim);CHKERRQ(ierr);CHKMEMQ;
    ierr = DefiantUpdateVolumeFactors(MySim);CHKERRQ(ierr);CHKMEMQ;
    ierr = DefiantUpdateViscosities(MySim);CHKERRQ(ierr);CHKMEMQ;

    /* interpolate pressure dependent properties */
    ierr = DefiantComputeRhoAndMuAtFaces(MySim);CHKERRQ(ierr);CHKMEMQ;
    ierr = DefiantComputeVolumeFactorsAtFaces(MySim);CHKERRQ(ierr);CHKMEMQ;

    /* Now w eare ready to update the saturations using improved IMPES */
    while(MySim->CurrentTimeS < MySim->CurrentTimeP + MySim->DeltaTP)
    {
      ierr = DefiantIMPES2PhUpdateSaturations(MySim);CHKERRQ(ierr);CHKMEMQ;
      /* update saturation based properties */
      ierr = DefiantUpdateRelativePermeability(MySim);CHKERRQ(ierr);CHKMEMQ;
      ierr = DefiantUpdatePcow(MySim);CHKERRQ(ierr);CHKMEMQ;
      /* update pw now that we have a new pcow */
      ierr = DefiantUpdatePwFromPcow(MySim);CHKERRQ(ierr);CHKMEMQ;
      /* with new pcow we have new pw, so new properties */
      ierr = DefiantUpdatePorosity(MySim);CHKERRQ(ierr);CHKMEMQ;
      ierr = DefiantUpdateDensity(MySim);CHKERRQ(ierr);CHKMEMQ;
      ierr = DefiantUpdateVolumeFactors(MySim);CHKERRQ(ierr);CHKMEMQ;
      ierr = DefiantUpdateViscosities(MySim);CHKERRQ(ierr);CHKMEMQ;
      /* update those new properties at the faces with interpolation */
      ierr = DefiantComputeRhoAndMuAtFaces(MySim);CHKERRQ(ierr);CHKMEMQ;
      ierr = DefiantComputeRelativePermsAtFaces(MySim);CHKERRQ(ierr);CHKMEMQ;
      ierr = DefiantComputeVolumeFactorsAtFaces(MySim);CHKERRQ(ierr);CHKMEMQ;

      /* Transmissibility calculation routines */
      ierr = DefiantComputeTransmissibilities(MySim);CHKERRQ(ierr);CHKMEMQ;

      /* Finally update the well information */
      ierr = DefiantIMPES2PhHandleWellsForSaturation(MySim);CHKERRQ(ierr);CHKMEMQ;

      ierr = PetscPrintf(PETSC_COMM_WORLD,"Output saturations is: ");
      ierr = VecView(MySim->Sw,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;

      ierr = PetscPrintf(PETSC_COMM_WORLD,"Pcow is: ");
      ierr = VecView(MySim->Pcow,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Pw is: ");
      ierr = VecView(MySim->Pw,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;

      ierr = PetscPrintf(PETSC_COMM_WORLD,"Krw is: ");
      ierr = VecView(MySim->Krw,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;


      ierr = PetscPrintf(PETSC_COMM_WORLD,"RelPermwx1p is: ");
      ierr = VecView(MySim->RelPermwx1p,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
      ierr = PetscPrintf(PETSC_COMM_WORLD,"RelPermwx1m is: ");
      ierr = VecView(MySim->RelPermwx1m,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Transmissibility Twx1p is: ");
      ierr = VecView(MySim->Twx1p,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Transmissibility Twx1m is: ");
      ierr = VecView(MySim->Twx1m,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;

      MySim->CurrentTimeS = MySim->CurrentTimeS + MySim->DeltaTS;
      ierr = DefiantViewSaturationsASCII(MySim);CHKERRQ(ierr);CHKMEMQ;
      break;
    }

    MySim->CurrentTimeP = MySim->CurrentTimeP + MySim->DeltaTP;
    break;
  }

  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "DefiantIMPES2PhDMMGSolve"
extern PetscErrorCode DefiantIMPES2PhDMMGSolve(BlackOilReservoirSimulation* MySim)
{
  PetscErrorCode ierr;
  PetscReal norm;

  PetscFunctionBegin;

  /* Set the current simulation to MySIM */
  CurrentSimulation = MySim;

  /* Set the KSP */
  ierr = DMMGSetKSP(MySim->SimDMMG, DefiantIMPES2PhDMMGComputeRHS, DefiantIMPES2PhDMMGComputeMatrix);CHKERRQ(ierr);
  /* Do the solve */
  ierr = DMMGSolve(MySim->SimDMMG);CHKERRQ(ierr);CHKMEMQ;
  ierr = MatMult(DMMGGetJ(MySim->SimDMMG), DMMGGetx(MySim->SimDMMG), DMMGGetr(MySim->SimDMMG));CHKERRQ(ierr);CHKMEMQ;
  ierr = VecAXPY(DMMGGetr(MySim->SimDMMG), -1.0, DMMGGetRHS(MySim->SimDMMG));CHKERRQ(ierr);CHKMEMQ;
  ierr = VecNorm(DMMGGetr(MySim->SimDMMG), NORM_2, &norm);CHKERRQ(ierr);CHKMEMQ;
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Residual norm %G\n", norm);CHKERRQ(ierr);CHKMEMQ;

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DefiantIMPES2PhDMMGComputeRHS"
PetscErrorCode DefiantIMPES2PhDMMGComputeRHS(DMMG dmmg, Vec b)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  /* Copy the RHS from the current simulation to the RHS */
  ierr = VecCopy(CurrentSimulation->RHS, b);CHKERRQ(ierr);CHKMEMQ;
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "DefiantIMPES2PhDMMGComputeMatrix"
PetscErrorCode DefiantIMPES2PhDMMGComputeMatrix(DMMG dmmg, Mat jac, Mat B)
{
  PetscErrorCode ierr;
  BlackOilReservoirSimulation* MySim = CurrentSimulation;
  PetscInt i, j, k, mx, my, mz, xm, ym, zm, xs, ys, zs;
  PetscScalar v[7];
  MatStencil row, col[7];
  /* Transmissibility For Oil*/
  PetscScalar ***LocalTox1p, ***LocalTox2p, ***LocalTox3p;
  PetscScalar ***LocalTox1m, ***LocalTox2m, ***LocalTox3m;
  /* Transmissibility For Water*/
  PetscScalar ***LocalTwx1p, ***LocalTwx2p, ***LocalTwx3p;
  PetscScalar ***LocalTwx1m, ***LocalTwx2m, ***LocalTwx3m;

  PetscFunctionBegin;
  /* Get dimensions and extents of the local vectors */
  ierr = DAGetInfo((DA)dmmg->dm, 0, &mx, &my, &mz, 0, 0, 0, 0, 0, 0, 0);CHKERRQ(ierr);
  ierr = DAGetCorners((DA)dmmg->dm, &xs, &ys, &zs, &xm, &ym, &zm);CHKERRQ(ierr);

  /* Grab the local data for the transmissibilities */
  ierr = DAVecGetArray(MySim->SimDA, MySim->Tox1m, &LocalTox1m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Tox1p, &LocalTox1p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Twx1m, &LocalTwx1m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Twx1p, &LocalTwx1p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Tox2m, &LocalTox2m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Tox2p, &LocalTox2p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Twx2m, &LocalTwx2m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Twx2p, &LocalTwx2p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Tox3m, &LocalTox3m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Tox3p, &LocalTox3p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Twx3m, &LocalTwx3m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Twx3p, &LocalTwx3p);CHKERRQ(ierr);

  for (k = zs; k < zs + zm; k++) {
    for (j = ys; j < ys + ym; j++) {
      for (i = xs; i < xs + xm; i++) {
        row.i = i;
        row.j = j;
        row.k = k;
        if (i == 0 || j == 0 || k == 0 || i == mx - 1 || j == my - 1 || k == mz
            - 1) {
          v[0] = 1.0;
          ierr = MatSetValuesStencil(B, 1, &row, 1, &row, v,
              INSERT_VALUES);
          CHKERRQ(ierr);
        } else {
          v[0] = -(LocalTox3m[k][j][i] + LocalTwx3m[k][j][i]);
          col[0].i = i;
          col[0].j = j;
          col[0].k = k - 1;
          v[1] = -(LocalTox2m[k][j][i] + LocalTwx2m[k][j][i]);
          col[1].i = i;
          col[1].j = j - 1;
          col[1].k = k;
          v[2] = -(LocalTox1m[k][j][i] + LocalTwx1m[k][j][i]);
          col[2].i = i - 1;
          col[2].j = j;
          col[2].k = k;
          /* moved v[3] down to sum over all other v's */
          v[4] = -(LocalTox1p[k][j][i] + LocalTwx1p[k][j][i]);
          col[4].i = i + 1;
          col[4].j = j;
          col[4].k = k;
          v[5] = -(LocalTox2p[k][j][i] + LocalTwx2p[k][j][i]);
          col[5].i = i;
          col[5].j = j + 1;
          col[5].k = k;
          v[6] = -(LocalTox3p[k][j][i] + LocalTwx3p[k][j][i]);
          col[6].i = i;
          col[6].j = j;
          col[6].k = k + 1;
          /* v[3] is the sum of all other column values */
          v[3] = (v[0] + v[1] + v[2] + v[4] + v[5] + v[6]);
          col[3].i = row.i;
          col[3].j = row.j;
          col[3].k = row.k;
          ierr = MatSetValuesStencil(B, 1, &row, 7, col, v,
              INSERT_VALUES);
          CHKERRQ(ierr);
        }
      }
    }
  }
  ierr = MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
