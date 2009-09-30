/*
 * RelativePermeability.c
 *
 *  Created on: Sep 6, 2009
 *      Author: yye00
 */

#include "Defiant.h"


#undef __FUNCT__
#define __FUNCT__ "DefiantUpdateRelativePermeability"
extern PetscErrorCode DefiantUpdateRelativePermeability(BlackOilReservoirSimulation* MySim)
{
  PetscErrorCode ierr;

  /* This obviously needs to be parametrised */
  ierr = DefiantComputeRelativePermsFromSatsCorey(MySim);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DefiantComputeRelativePermsFromSatsCorey"
extern PetscErrorCode DefiantComputeRelativePermsFromSatsCorey(
    BlackOilReservoirSimulation* MySim) {
  PetscErrorCode ierr;
  PetscInt i, j, k, mx, my, mz, xm, ym, zm, xs, ys, zs;
  PetscScalar ***LocalFlowMask;
  PetscScalar ***LocalKro, ***LocalKrw;
  PetscScalar ***LocalSo, ***LocalSw;
  PetscScalar Snw;

  PetscFunctionBegin;
  /* Get dimensions and extents of the local vectors */
  ierr = DAGetInfo(MySim->SimDA, 0, &mx, &my, &mz, 0, 0, 0, 0, 0, 0, 0);CHKERRQ(ierr);
  ierr = DAGetCorners(MySim->SimDA, &xs, &ys, &zs, &xm, &ym, &zm);CHKERRQ(ierr);
  /* Grab the data for the flow field */
  ierr = DAVecGetArray(MySim->SimDA, MySim->FlowMask, &LocalFlowMask);CHKERRQ(ierr);
  /* Grab the local data for Saturations */
  ierr = DAVecGetArray(MySim->SimDA, MySim->So, &LocalSo);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Sw, &LocalSw);CHKERRQ(ierr);
  /* Grab the local data for the relative permeabilities */
  ierr = DAVecGetArray(MySim->SimDA, MySim->Kro, &LocalKro);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Krw, &LocalKrw);CHKERRQ(ierr);

  for (k = zs; k < zs + zm; k++) {
    for (j = ys; j < ys + ym; j++) {
      for (i = xs; i < xs + xm; i++) {
        if (i == 0 || j == 0 || k == 0 || i == mx - 1 || j == my - 1 || k == mz
            - 1) {
        } else if (ABS(LocalFlowMask[k][j][i]-FLUID_FLOW) < EPSILON) {
          Snw = (LocalSw[k][j][i] - MySim->Swc) / (1.0 - MySim->Swc);
          LocalKro[k][j][i] = pow(1.0 - Snw, 2.0) * (1.0 - pow(Snw, 2.0));
          LocalKrw[k][j][i] = pow(Snw, 4.0);
        }
      }
    }
  }

  /* Restore the new arrays to their rightful place */
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Kro, &LocalKro);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Krw, &LocalKrw);CHKERRQ(ierr);

  /* Begin Assembly for vectors */
  ierr = VecAssemblyBegin(MySim->Kro);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->Krw);CHKERRQ(ierr);

  /* And end Assembly */
  ierr = VecAssemblyEnd(MySim->Kro);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->Krw);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DefiantComputeRelativePermsFromSatsNaar"
extern PetscErrorCode DefiantComputeRelativePermsFromSatsNaar(
    BlackOilReservoirSimulation* MySim) {
  PetscErrorCode ierr;
  PetscInt i, j, k, mx, my, mz, xm, ym, zm, xs, ys, zs;
  PetscScalar ***LocalFlowMask;
  PetscScalar ***LocalKro, ***LocalKrw;
  PetscScalar ***LocalSo, ***LocalSw;
  PetscScalar Snw;

  PetscFunctionBegin;
  /* Get dimensions and extents of the local vectors */
  ierr = DAGetInfo(MySim->SimDA, 0, &mx, &my, &mz, 0, 0, 0, 0, 0, 0, 0);CHKERRQ(ierr);
  ierr = DAGetCorners(MySim->SimDA, &xs, &ys, &zs, &xm, &ym, &zm);CHKERRQ(ierr);
  /* Grab the data for the flow field */
  ierr = DAVecGetArray(MySim->SimDA, MySim->FlowMask, &LocalFlowMask);CHKERRQ(ierr);
  /* Grab the local data for Saturations */
  ierr = DAVecGetArray(MySim->SimDA, MySim->So, &LocalSo);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Sw, &LocalSw);CHKERRQ(ierr);
  /* Grab the local data for the relative permeabilities */
  ierr = DAVecGetArray(MySim->SimDA, MySim->Kro, &LocalKro);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Krw, &LocalKrw);CHKERRQ(ierr);

  for (k = zs; k < zs + zm; k++) {
    for (j = ys; j < ys + ym; j++) {
      for (i = xs; i < xs + xm; i++) {
        if (i == 0 || j == 0 || k == 0 || i == mx - 1 || j == my - 1 || k == mz
            - 1) {
        } else if (ABS(LocalFlowMask[k][j][i]-FLUID_FLOW) < EPSILON) {
          Snw = (LocalSw[k][j][i] - MySim->Swc) / (1.0 - MySim->Swc);
          LocalKro[k][j][i] = pow(1.0 - 2.0 * Snw, 1.5) * (2 - pow(1 - 2.0
              * Snw, 0.5));
          LocalKrw[k][j][i] = pow(Snw, 4);
        }
      }
    }
  }

  /* Restore the new arrays to their reightful place */
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Kro, &LocalKro);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Krw, &LocalKrw);CHKERRQ(ierr);

  /* Begin Assembly for vectors */
  ierr = VecAssemblyBegin(MySim->Kro);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->Krw);CHKERRQ(ierr);

  /* And end Assembly */
  ierr = VecAssemblyEnd(MySim->Kro);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->Krw);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DefiantComputeRelativePermsFromSatsWygal"
extern PetscErrorCode DefiantComputeRelativePermsFromSatsWygal(
    BlackOilReservoirSimulation* MySim) {
  PetscErrorCode ierr;
  PetscInt i, j, k, mx, my, mz, xm, ym, zm, xs, ys, zs;
  PetscScalar ***LocalFlowMask;
  PetscScalar ***LocalKro, ***LocalKrw, ***LocalKrg;
  PetscScalar ***LocalSo, ***LocalSw, ***LocalSg;

  PetscFunctionBegin;
  /* Get dimensions and extents of the local vectors */
  ierr = DAGetInfo(MySim->SimDA, 0, &mx, &my, &mz, 0, 0, 0, 0, 0, 0, 0);CHKERRQ(ierr);
  ierr = DAGetCorners(MySim->SimDA, &xs, &ys, &zs, &xm, &ym, &zm);CHKERRQ(ierr);
  /* Grab the data for the flow field */
  ierr = DAVecGetArray(MySim->SimDA, MySim->FlowMask, &LocalFlowMask);CHKERRQ(ierr);
  /* Grab the local data for Saturations */
  ierr = DAVecGetArray(MySim->SimDA, MySim->So, &LocalSo);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Sw, &LocalSw);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Sg, &LocalSg);CHKERRQ(ierr);
  /* Grab the local data for the relative permeabilities */
  ierr = DAVecGetArray(MySim->SimDA, MySim->Kro, &LocalKro);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Krw, &LocalKrw);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Krg, &LocalKrg);CHKERRQ(ierr);

  for (k = zs; k < zs + zm; k++) {
    for (j = ys; j < ys + ym; j++) {
      for (i = xs; i < xs + xm; i++) {
        if (i == 0 || j == 0 || k == 0 || i == mx - 1 || j == my - 1 || k == mz
            - 1) {
        } else if (ABS(LocalFlowMask[k][j][i]-FLUID_FLOW) < EPSILON) {
          LocalKro[k][j][i] = pow(LocalSo[k][j][i], 3) * (1 - LocalSg[k][j][i]
              + 2.0 * LocalSw[k][j][i] - 3.0 * MySim->Swc) / pow(1.0
              - MySim->Swc, 4);
          LocalKrw[k][j][i] = pow((LocalSw[k][j][i] - MySim->Swc) / (1.0
              - MySim->Swc), 4);
          LocalKrg[k][j][i] = pow(LocalSg[k][j][i], 3) * (2 - LocalSg[k][j][i]
              + 2.0 * MySim->Swc) / pow(1.0 - MySim->Swc, 4);

          if (LocalSo[k][j][i] <= MySim->Sor)
            LocalKro[k][j][i] = 0.0;
          if (LocalSg[k][j][i] <= MySim->Sgc)
            LocalKrg[k][j][i] = 0.0;
        }
      }
    }
  }

  /* Restore the new arrays to their reightful place */
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Kro, &LocalKro);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Krw, &LocalKrw);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Krg, &LocalKrg);CHKERRQ(ierr);

  /* Begin Assembly for vectors */
  ierr = VecAssemblyBegin(MySim->Kro);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->Krw);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->Krg);CHKERRQ(ierr);

  /* And end Assembly */
  ierr = VecAssemblyEnd(MySim->Kro);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->Krw);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->Krg);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DefiantComputeRelativePermsFromSatsStonesI"
extern PetscErrorCode DefiantComputeRelativePermsFromSatsStonesI(
    BlackOilReservoirSimulation* MySim) {
  PetscErrorCode ierr;
  PetscInt i, j, k, mx, my, mz, xm, ym, zm, xs, ys, zs;
  PetscScalar Sno, Snw, Sng, Betaw, Betag;
  PetscScalar ***LocalFlowMask;
  PetscScalar ***LocalKro, ***LocalKrw, ***LocalKrg;
  PetscScalar ***LocalSo, ***LocalSw, ***LocalSg;
  PetscScalar Krow, Krog;

  PetscFunctionBegin;
  /* Get dimensions and extents of the local vectors */
  ierr = DAGetInfo(MySim->SimDA, 0, &mx, &my, &mz, 0, 0, 0, 0, 0, 0, 0);CHKERRQ(ierr);
  ierr = DAGetCorners(MySim->SimDA, &xs, &ys, &zs, &xm, &ym, &zm);CHKERRQ(ierr);
  /* Grab the data for the flow field */
  ierr = DAVecGetArray(MySim->SimDA, MySim->FlowMask, &LocalFlowMask);CHKERRQ(ierr);
  /* Grab the local data for Saturations */
  ierr = DAVecGetArray(MySim->SimDA, MySim->So, &LocalSo);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Sw, &LocalSw);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Sg, &LocalSg);CHKERRQ(ierr);
  /* Grab the local data for the relative permeabilities */
  ierr = DAVecGetArray(MySim->SimDA, MySim->Kro, &LocalKro);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Krw, &LocalKrw);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Krg, &LocalKrg);CHKERRQ(ierr);

  for (k = zs; k < zs + zm; k++) {
    for (j = ys; j < ys + ym; j++) {
      for (i = xs; i < xs + xm; i++) {
        if (i == 0 || j == 0 || k == 0 || i == mx - 1 || j == my - 1 || k == mz
            - 1) {
        } else if (ABS(LocalFlowMask[k][j][i]-FLUID_FLOW) < EPSILON) {
          if (LocalSo[k][j][i] >= MySim->Sor)
            Sno = (LocalSo[k][j][i] - MySim->Sor) / (1.0 - MySim->Swc
                - MySim->Sor);
          else
            Sno = 0.0;

          if (LocalSw[k][j][i] >= MySim->Swc)
            Snw = (LocalSw[k][j][i] - MySim->Swc) / (1.0 - MySim->Swc
                - MySim->Sor);
          else
            Snw = 0.0;

          Sng = LocalSg[k][j][i] / (1.0 - MySim->Swc - MySim->Sor);
          /* Find the relative oil-water and oil-gas permeabilities */
          ierr = DefiantComputeKrow(&Krow, LocalSw[k][j][i]);
          CHKERRQ(ierr);
          ierr = DefiantComputeKrog(&Krog, LocalSg[k][j][i]);
          CHKERRQ(ierr);
          /* Compute Betas */
          Betaw = Krow / (1.0 - Snw);
          Betag = Krog / (1.0 - Sng);
          /* Compute relative permeability of oil */
          LocalKro[k][j][i] = Sno * Betaw * Betag;
          /* Retrieve the relative permeability for water and gas from their functions */
          ierr = DefiantComputeKrw(&(LocalKrw[k][j][i]), LocalSw[k][j][i]);
          ierr = DefiantComputeKrg(&(LocalKrg[k][j][i]), LocalSg[k][j][i]);
        }
      }
    }
  }

  /* Restore the new arrays to their reightful place */
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Kro, &LocalKro);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Krw, &LocalKrw);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Krg, &LocalKrg);CHKERRQ(ierr);

  /* Begin Assembly for vectors */
  ierr = VecAssemblyBegin(MySim->Kro);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->Krw);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->Krg);CHKERRQ(ierr);

  /* And end Assembly */
  ierr = VecAssemblyEnd(MySim->Kro);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->Krw);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->Krg);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DefiantComputeRelativePermsFromSatsStonesII"
extern PetscErrorCode DefiantComputeRelativePermsFromSatsStonesII(
    BlackOilReservoirSimulation* MySim) {
  PetscErrorCode ierr;
  PetscInt i, j, k, mx, my, mz, xm, ym, zm, xs, ys, zs;
  PetscScalar ***LocalFlowMask;
  PetscScalar ***LocalKro, ***LocalKrw, ***LocalKrg;
  PetscScalar ***LocalSo, ***LocalSw, ***LocalSg;
  PetscScalar Krow, Krog;

  PetscFunctionBegin;
  /* Get dimensions and extents of the local vectors */
  ierr = DAGetInfo(MySim->SimDA, 0, &mx, &my, &mz, 0, 0, 0, 0, 0, 0, 0);CHKERRQ(ierr);
  ierr = DAGetCorners(MySim->SimDA, &xs, &ys, &zs, &xm, &ym, &zm);CHKERRQ(ierr);
  /* Grab the data for the flow field */
  ierr = DAVecGetArray(MySim->SimDA, MySim->FlowMask, &LocalFlowMask);CHKERRQ(ierr);
  /* Grab the local data for Saturations */
  ierr = DAVecGetArray(MySim->SimDA, MySim->So, &LocalSo);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Sw, &LocalSw);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Sg, &LocalSg);CHKERRQ(ierr);
  /* Grab the local data for the relative permeabilities */
  ierr = DAVecGetArray(MySim->SimDA, MySim->Kro, &LocalKro);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Krw, &LocalKrw);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Krg, &LocalKrg);CHKERRQ(ierr);

  for (k = zs; k < zs + zm; k++) {
    for (j = ys; j < ys + ym; j++) {
      for (i = xs; i < xs + xm; i++) {
        if (i == 0 || j == 0 || k == 0 || i == mx - 1 || j == my - 1 || k == mz
            - 1) {
        } else if (ABS(LocalFlowMask[k][j][i]-FLUID_FLOW) < EPSILON) {
          /* Find the relative oil-water and oil-gas permeabilities */
          ierr = DefiantComputeKrow(&Krow, LocalSw[k][j][i]);
          CHKERRQ(ierr);
          ierr = DefiantComputeKrog(&Krog, LocalSg[k][j][i]);
          CHKERRQ(ierr);
          /* Retrieve the relative permeability for water and gas from their functions */
          ierr = DefiantComputeKrw(&(LocalKrw[k][j][i]), LocalSw[k][j][i]);
          ierr = DefiantComputeKrg(&(LocalKrg[k][j][i]), LocalSg[k][j][i]);

          /* Compute relative permeability of oil */
          LocalKro[k][j][i] = (Krow + LocalKrw[k][j][i]) * (Krog
              + LocalKrg[k][j][i]) - (LocalKrw[k][j][i] + LocalKrg[k][j][i]);

        }
      }
    }
  }

  /* Restore the new arrays to their reightful place */
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Kro, &LocalKro);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Krw, &LocalKrw);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Krg, &LocalKrg);CHKERRQ(ierr);

  /* Begin Assembly for vectors */
  ierr = VecAssemblyBegin(MySim->Kro);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->Krw);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->Krg);CHKERRQ(ierr);

  /* And end Assembly */
  ierr = VecAssemblyEnd(MySim->Kro);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->Krw);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->Krg);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}


/*
 *
 * These functions are used to compute relative permeabilities from saturations
 *
 */

#undef __FUNCT__
#define __FUNCT__ "DefiantComputeKrw"
extern PetscErrorCode DefiantComputeKrw(PetscScalar * Krw, PetscScalar Sw) {
  PetscFunctionBegin;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DefiantComputeKrw"
extern PetscErrorCode DefiantComputeKrg(PetscScalar * Krg, PetscScalar Sg) {
  PetscFunctionBegin;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DefiantComputeKrow"
extern PetscErrorCode DefiantComputeKrow(PetscScalar * Krow, PetscScalar Sw) {
  PetscFunctionBegin;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DefiantComputeKrog"
extern PetscErrorCode DefiantComputeKrog(PetscScalar * Krog, PetscScalar Sg) {
  PetscFunctionBegin;
  PetscFunctionReturn(0);
}
