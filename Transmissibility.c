/*
 * Transmissibility.c
 *
 *  Created on: Sep 6, 2009
 *      Author: yye00
 */

#include "Defiant.h"


#undef __FUNCT__
#define __FUNCT__ "DefiantComputeTransmissibilities"
extern PetscErrorCode DefiantComputeTransmissibilities(
    BlackOilReservoirSimulation* MySim) {
  PetscErrorCode ierr;
  PetscInt i, j, k, mx, my, mz, xm, ym, zm, xs, ys, zs;
  PetscScalar ***LocalFlowMask;
  /* Area * Permeability divided by height */
  PetscScalar ***LocalAKHx1m, ***LocalAKHx2m, ***LocalAKHx3m;
  PetscScalar ***LocalAKHx1p, ***LocalAKHx2p, ***LocalAKHx3p;
  /* Local viscosity at the faces */
  PetscScalar ***LocalMuox1m, ***LocalMuox1p, ***LocalMuwx1m, ***LocalMuwx1p,
      ***LocalMugx1m, ***LocalMugx1p;
  PetscScalar ***LocalMuox2m, ***LocalMuox2p, ***LocalMuwx2m, ***LocalMuwx2p,
      ***LocalMugx2m, ***LocalMugx2p;
  PetscScalar ***LocalMuox3m, ***LocalMuox3p, ***LocalMuwx3m, ***LocalMuwx3p,
      ***LocalMugx3m, ***LocalMugx3p;
  /* Relative Permeabilities at the faces */
  PetscScalar ***LocalRelPermox1m, ***LocalRelPermox1p, ***LocalRelPermox2m,
      ***LocalRelPermox2p, ***LocalRelPermox3m, ***LocalRelPermox3p;
  PetscScalar ***LocalRelPermwx1m, ***LocalRelPermwx1p, ***LocalRelPermwx2m,
      ***LocalRelPermwx2p, ***LocalRelPermwx3m, ***LocalRelPermwx3p;
  PetscScalar ***LocalRelPermgx1m, ***LocalRelPermgx1p, ***LocalRelPermgx2m,
      ***LocalRelPermgx2p, ***LocalRelPermgx3m, ***LocalRelPermgx3p;
  /* Volume factors at cell center */
  PetscScalar ***LocalBo, ***LocalBw, ***LocalBg;
  /* Volume factors at cell faces */
  PetscScalar ***LocalBox1p, ***LocalBox2p, ***LocalBox3p;
  PetscScalar ***LocalBox1m, ***LocalBox2m, ***LocalBox3m;
  PetscScalar ***LocalBwx1p, ***LocalBwx2p, ***LocalBwx3p;
  PetscScalar ***LocalBwx1m, ***LocalBwx2m, ***LocalBwx3m;
  PetscScalar ***LocalBgx1p, ***LocalBgx2p, ***LocalBgx3p;
  PetscScalar ***LocalBgx1m, ***LocalBgx2m, ***LocalBgx3m;
  /* Transmissibility For Oil*/
  PetscScalar ***LocalTox1p, ***LocalTox2p, ***LocalTox3p;
  PetscScalar ***LocalTox1m, ***LocalTox2m, ***LocalTox3m;
  /* Transmissibility For Water*/
  PetscScalar ***LocalTwx1p, ***LocalTwx2p, ***LocalTwx3p;
  PetscScalar ***LocalTwx1m, ***LocalTwx2m, ***LocalTwx3m;
  /* Transmissibility For Gas*/
  PetscScalar ***LocalTgx1p, ***LocalTgx2p, ***LocalTgx3p;
  PetscScalar ***LocalTgx1m, ***LocalTgx2m, ***LocalTgx3m;

  PetscFunctionBegin;
  /* Get dimensions and extents of the local vectors */
  ierr = DAGetInfo(MySim->SimDA, 0, &mx, &my, &mz, 0, 0, 0, 0, 0, 0, 0);CHKERRQ(ierr);
  ierr = DAGetCorners(MySim->SimDA, &xs, &ys, &zs, &xm, &ym, &zm);CHKERRQ(ierr);
  /* Grab the data for the flow field */
  ierr = DAVecGetArray(MySim->SimDA, MySim->FlowMask, &LocalFlowMask);CHKERRQ(ierr);
  /* Grab the local data for area * permeability divided by height */
  ierr = DAVecGetArray(MySim->SimDA, MySim->AKHx1m, &LocalAKHx1m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->AKHx2m, &LocalAKHx2m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->AKHx3m, &LocalAKHx3m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->AKHx1p, &LocalAKHx1p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->AKHx2p, &LocalAKHx2p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->AKHx3p, &LocalAKHx3p);CHKERRQ(ierr);
  /* Grab the local data for the face viscosities */
  ierr = DAVecGetArray(MySim->SimDA, MySim->Muox1m, &LocalMuox1m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Muox1p, &LocalMuox1p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Muwx1m, &LocalMuwx1m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Muwx1p, &LocalMuwx1p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Mugx1m, &LocalMugx1m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Mugx1p, &LocalMugx1p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Muox2m, &LocalMuox2m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Muox2p, &LocalMuox2p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Muwx2m, &LocalMuwx2m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Muwx2p, &LocalMuwx2p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Mugx2m, &LocalMugx2m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Mugx2p, &LocalMugx2p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Muox3m, &LocalMuox3m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Muox3p, &LocalMuox3p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Muwx3m, &LocalMuwx3m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Muwx3p, &LocalMuwx3p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Mugx3m, &LocalMugx3m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Mugx3p, &LocalMugx3p);CHKERRQ(ierr);
  /* Grab the local data for the face Relative Permeabilities */
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
  /* Grab the local data for volume factors at the cell centers */
  ierr = DAVecGetArray(MySim->SimDA, MySim->Bo, &LocalBo);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Bw, &LocalBw);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Bg, &LocalBg);CHKERRQ(ierr);
  /* Grab the local data for Volume factors at the faces */
  ierr = DAVecGetArray(MySim->SimDA, MySim->Box1p, &LocalBox1p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Box2p, &LocalBox2p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Box3p, &LocalBox3p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Box1m, &LocalBox1m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Box2m, &LocalBox2m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Box3m, &LocalBox3m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Bwx1p, &LocalBwx1p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Bwx2p, &LocalBwx2p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Bwx3p, &LocalBwx3p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Bwx1m, &LocalBwx1m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Bwx2m, &LocalBwx2m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Bwx3m, &LocalBwx3m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Bgx1p, &LocalBgx1p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Bgx2p, &LocalBgx2p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Bgx3p, &LocalBgx3p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Bgx1m, &LocalBgx1m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Bgx2m, &LocalBgx2m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Bgx3m, &LocalBgx3m);CHKERRQ(ierr);
  /* Grab the local data for the transmissibilities */
  ierr = DAVecGetArray(MySim->SimDA, MySim->Tox1m, &LocalTox1m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Tox1p, &LocalTox1p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Twx1m, &LocalTwx1m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Twx1p, &LocalTwx1p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Tgx1m, &LocalTgx1m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Tgx1p, &LocalTgx1p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Tox2m, &LocalTox2m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Tox2p, &LocalTox2p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Twx2m, &LocalTwx2m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Twx2p, &LocalTwx2p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Tgx2m, &LocalTgx2m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Tgx2p, &LocalTgx2p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Tox3m, &LocalTox3m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Tox3p, &LocalTox3p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Twx3m, &LocalTwx3m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Twx3p, &LocalTwx3p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Tgx3m, &LocalTgx3m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Tgx3p, &LocalTgx3p);CHKERRQ(ierr);

  for (k = zs; k < zs + zm; k++) {
    for (j = ys; j < ys + ym; j++) {
      for (i = xs; i < xs + xm; i++) {
        if (i == 0 || j == 0 || k == 0 || i == mx - 1 || j == my - 1 || k == mz
            - 1) {
        } else if (ABS(LocalFlowMask[k][j][i]-FLUID_FLOW) < EPSILON) {
          /* Transmissibility */
          LocalTox1m[k][j][i] = LocalRelPermox1m[k][j][i]
              * LocalAKHx1m[k][j][i] / LocalMuox1m[k][j][i] / LocalBox1m[k][j][i];
          LocalTox1p[k][j][i] = LocalRelPermox1p[k][j][i]
              * LocalAKHx1p[k][j][i] / LocalMuox1p[k][j][i] / LocalBox1p[k][j][i];
          LocalTox2m[k][j][i] = LocalRelPermox2m[k][j][i]
              * LocalAKHx2m[k][j][i] / LocalMuox2m[k][j][i] / LocalBox2m[k][j][i];
          LocalTox2p[k][j][i] = LocalRelPermox2p[k][j][i]
              * LocalAKHx2p[k][j][i] / LocalMuox2p[k][j][i] / LocalBox2p[k][j][i];
          LocalTox3m[k][j][i] = LocalRelPermox3m[k][j][i]
              * LocalAKHx3m[k][j][i] / LocalMuox3m[k][j][i] / LocalBox3m[k][j][i];
          LocalTox3p[k][j][i] = LocalRelPermox3p[k][j][i]
              * LocalAKHx3p[k][j][i] / LocalMuox3p[k][j][i] / LocalBox3p[k][j][i];

          LocalTwx1m[k][j][i] = LocalRelPermwx1m[k][j][i]
              * LocalAKHx1m[k][j][i] / LocalMuwx1m[k][j][i] / LocalBwx1m[k][j][i];
          LocalTwx1p[k][j][i] = LocalRelPermwx1p[k][j][i]
              * LocalAKHx1p[k][j][i] / LocalMuwx1p[k][j][i] / LocalBwx1p[k][j][i];
          LocalTwx2m[k][j][i] = LocalRelPermwx2m[k][j][i]
              * LocalAKHx2m[k][j][i] / LocalMuwx2m[k][j][i] / LocalBwx2m[k][j][i];
          LocalTwx2p[k][j][i] = LocalRelPermwx2p[k][j][i]
              * LocalAKHx2p[k][j][i] / LocalMuwx2p[k][j][i] / LocalBwx2p[k][j][i];
          LocalTwx3m[k][j][i] = LocalRelPermwx3m[k][j][i]
              * LocalAKHx3m[k][j][i] / LocalMuwx3m[k][j][i] / LocalBwx3m[k][j][i];
          LocalTwx3p[k][j][i] = LocalRelPermwx3p[k][j][i]
              * LocalAKHx3p[k][j][i] / LocalMuwx3p[k][j][i] / LocalBwx3p[k][j][i];

          LocalTgx1m[k][j][i] = LocalRelPermgx1m[k][j][i]
              * LocalAKHx1m[k][j][i] / LocalMugx1m[k][j][i] / LocalBgx1m[k][j][i];
          LocalTgx1p[k][j][i] = LocalRelPermgx1p[k][j][i]
              * LocalAKHx1p[k][j][i] / LocalMugx1p[k][j][i] / LocalBgx1p[k][j][i];
          LocalTgx2m[k][j][i] = LocalRelPermgx2m[k][j][i]
              * LocalAKHx2m[k][j][i] / LocalMugx2m[k][j][i] / LocalBgx2m[k][j][i];
          LocalTgx2p[k][j][i] = LocalRelPermgx2p[k][j][i]
              * LocalAKHx2p[k][j][i] / LocalMugx2p[k][j][i] / LocalBgx2p[k][j][i];
          LocalTgx3m[k][j][i] = LocalRelPermgx3m[k][j][i]
              * LocalAKHx3m[k][j][i] / LocalMugx3m[k][j][i] / LocalBgx3m[k][j][i];
          LocalTgx3p[k][j][i] = LocalRelPermgx3p[k][j][i]
              * LocalAKHx3p[k][j][i] / LocalMugx3p[k][j][i] / LocalBgx3p[k][j][i];
        }
      }
    }
  }

  /* Restore the new arrays to their reightful place */
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Tox1m, &LocalTox1m);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Tox1p, &LocalTox1p);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Twx1m, &LocalTwx1m);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Twx1p, &LocalTwx1p);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Tgx1m, &LocalTgx1m);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Tgx1p, &LocalTgx1p);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Tox2m, &LocalTox2m);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Tox2p, &LocalTox2p);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Twx2m, &LocalTwx2m);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Twx2p, &LocalTwx2p);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Tgx2m, &LocalTgx2m);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Tgx2p, &LocalTgx2p);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Tox3m, &LocalTox3m);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Tox3p, &LocalTox3p);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Twx3m, &LocalTwx3m);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Twx3p, &LocalTwx3p);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Tgx3m, &LocalTgx3m);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Tgx3p, &LocalTgx3p);CHKERRQ(ierr);

  /* Begin Assembly for vectors */
  ierr = VecAssemblyBegin(MySim->Tox1p);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->Tox2p);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->Tox3p);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->Tox1m);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->Tox2m);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->Tox3m);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->Twx1p);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->Twx2p);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->Twx3p);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->Twx1m);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->Twx2m);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->Twx3m);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->Tgx1p);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->Tgx2p);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->Tgx3p);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->Tgx1m);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->Tgx2m);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->Tgx3m);CHKERRQ(ierr);

  /* And end Assembly */
  ierr = VecAssemblyEnd(MySim->Tox1p);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->Tox2p);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->Tox3p);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->Tox1m);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->Tox2m);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->Tox3m);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->Twx1p);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->Twx2p);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->Twx3p);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->Twx1m);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->Twx2m);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->Twx3m);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->Tgx1p);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->Tgx2p);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->Tgx3p);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->Tgx1m);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->Tgx2m);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->Tgx3m);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
