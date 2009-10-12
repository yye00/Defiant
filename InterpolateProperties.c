/*
 * InterpolateProperties.c
 *
 *  Created on: Sep 6, 2009
 *      Author: yye00
 */

#include "Defiant.h"

#undef __FUNCT__
#define __FUNCT__ "DefiantComputeKAByHAtFaces"
PetscErrorCode DefiantComputeKAByHAtFaces(BlackOilReservoirSimulation* MySim) {
  PetscErrorCode ierr;
  PetscInt i, j, k, mx, my, mz, xm, ym, zm, xs, ys, zs;
  PetscScalar ***LocalFlowMask;
  PetscScalar ***LocalA1, ***LocalA2, ***LocalA3; /* Areas          */
  PetscScalar ***LocalK11, ***LocalK22, ***LocalK33; /* Permeabilities */
  PetscScalar ***Localh1, ***Localh2, ***Localh3; /* heights        */
  /* Area * Permeability divided by height */
  PetscScalar ***LocalAKHx1m, ***LocalAKHx2m, ***LocalAKHx3m;
  PetscScalar ***LocalAKHx1p, ***LocalAKHx2p, ***LocalAKHx3p;

  /* Local Temporary vectors */
  Vec  vecLocalFlowMask;
  Vec  vecLocalA1, vecLocalA2, vecLocalA3; /* Areas          */
  Vec  vecLocalK11, vecLocalK22, vecLocalK33; /* Permeabilities */
  Vec  vecLocalh1, vecLocalh2, vecLocalh3; /* heights        */

  PetscFunctionBegin;

  /* Get dimensions and extents of the local vectors */
  ierr = DAGetInfo(MySim->SimDA, 0, &mx, &my, &mz, 0, 0, 0, 0, 0, 0, 0);CHKERRQ(ierr);
  ierr = DAGetCorners(MySim->SimDA, &xs, &ys, &zs, &xm, &ym, &zm);CHKERRQ(ierr);

  /* Get the local vectors */
  ierr = DAGetLocalVector(MySim->SimDA, &vecLocalFlowMask);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGetLocalVector(MySim->SimDA, &vecLocalA1);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGetLocalVector(MySim->SimDA, &vecLocalA2);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGetLocalVector(MySim->SimDA, &vecLocalA3);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGetLocalVector(MySim->SimDA, &vecLocalK11);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGetLocalVector(MySim->SimDA, &vecLocalK22);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGetLocalVector(MySim->SimDA, &vecLocalK33);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGetLocalVector(MySim->SimDA, &vecLocalh1);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGetLocalVector(MySim->SimDA, &vecLocalh2);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGetLocalVector(MySim->SimDA, &vecLocalh3);CHKERRQ(ierr);CHKMEMQ;

  /* scatter global vectors to local ones */
  ierr = DAGlobalToLocalBegin(MySim->SimDA, MySim-> FlowMask  ,INSERT_VALUES, vecLocalFlowMask  );CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGlobalToLocalEnd(MySim->SimDA, MySim-> FlowMask  ,INSERT_VALUES, vecLocalFlowMask  );CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGlobalToLocalBegin(MySim->SimDA, MySim-> A1  ,INSERT_VALUES, vecLocalA1  );CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGlobalToLocalEnd(MySim->SimDA, MySim-> A1  ,INSERT_VALUES, vecLocalA1  );CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGlobalToLocalBegin(MySim->SimDA, MySim-> A2  ,INSERT_VALUES, vecLocalA2  );CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGlobalToLocalEnd(MySim->SimDA, MySim-> A2  ,INSERT_VALUES, vecLocalA2  );CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGlobalToLocalBegin(MySim->SimDA, MySim-> A3  ,INSERT_VALUES, vecLocalA3  );CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGlobalToLocalEnd(MySim->SimDA, MySim-> A3  ,INSERT_VALUES, vecLocalA3  );CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGlobalToLocalBegin(MySim->SimDA, MySim-> K11 ,INSERT_VALUES, vecLocalK11 );CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGlobalToLocalEnd(MySim->SimDA, MySim-> K11 ,INSERT_VALUES, vecLocalK11 );CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGlobalToLocalBegin(MySim->SimDA, MySim-> K22 ,INSERT_VALUES, vecLocalK22 );CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGlobalToLocalEnd(MySim->SimDA, MySim-> K22 ,INSERT_VALUES, vecLocalK22 );CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGlobalToLocalBegin(MySim->SimDA, MySim-> K33 ,INSERT_VALUES, vecLocalK33 );CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGlobalToLocalEnd(MySim->SimDA, MySim-> K33 ,INSERT_VALUES, vecLocalK33 );CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGlobalToLocalBegin(MySim->SimDA, MySim-> h1  ,INSERT_VALUES, vecLocalh1  );CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGlobalToLocalEnd(MySim->SimDA, MySim-> h1  ,INSERT_VALUES, vecLocalh1  );CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGlobalToLocalBegin(MySim->SimDA, MySim-> h2  ,INSERT_VALUES, vecLocalh2  );CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGlobalToLocalEnd(MySim->SimDA, MySim-> h2  ,INSERT_VALUES, vecLocalh2  );CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGlobalToLocalBegin(MySim->SimDA, MySim-> h3  ,INSERT_VALUES, vecLocalh3  );CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGlobalToLocalEnd(MySim->SimDA, MySim-> h3  ,INSERT_VALUES, vecLocalh3  );CHKERRQ(ierr);CHKMEMQ;

  /* These variables need ghost-zones */
  ierr = DAVecGetArray(MySim->SimDA, vecLocalFlowMask, &LocalFlowMask);CHKERRQ(ierr);
  /* Grab the local data for areas */
  ierr = DAVecGetArray(MySim->SimDA, vecLocalA1, &LocalA1);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, vecLocalA2, &LocalA2);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, vecLocalA3, &LocalA3);CHKERRQ(ierr);
  /* Grab the local data for the permeabilities */
  ierr = DAVecGetArray(MySim->SimDA, vecLocalK11, &LocalK11);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, vecLocalK22, &LocalK22);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, vecLocalK33, &LocalK33);CHKERRQ(ierr);
  /* Grab the local data for heights */
  ierr = DAVecGetArray(MySim->SimDA, vecLocalh1, &Localh1);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, vecLocalh2, &Localh2);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, vecLocalh3, &Localh3);CHKERRQ(ierr);

  /* Grab the local data for area * permeability divided by height */
  ierr = DAVecGetArray(MySim->SimDA, MySim->AKHx1m, &LocalAKHx1m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->AKHx2m, &LocalAKHx2m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->AKHx3m, &LocalAKHx3m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->AKHx1p, &LocalAKHx1p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->AKHx2p, &LocalAKHx2p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->AKHx3p, &LocalAKHx3p);CHKERRQ(ierr);

  for (k = zs; k < zs + zm; k++) {
    for (j = ys; j < ys + ym; j++) {
      for (i = xs; i < xs + xm; i++) {
        if (i == 0 || j == 0 || k == 0 || i == mx - 1 || j == my - 1 || k == mz
            - 1) {
          LocalAKHx1m[k][j][i] = 0.0;
          LocalAKHx1p[k][j][i] = 0.0;
          LocalAKHx2m[k][j][i] = 0.0;
          LocalAKHx2p[k][j][i] = 0.0;
          LocalAKHx3m[k][j][i] = 0.0;
          LocalAKHx3p[k][j][i] = 0.0;
        } else if (ABS(LocalFlowMask[k][j][i]-FLUID_FLOW) < EPSILON) {
          LocalAKHx1m[k][j][i] = 2.0 * LocalA1[k][j][i] * LocalK11[k][j][i]
              * LocalA1[k][j][i - 1] * LocalK11[k][j][i - 1]
              / (LocalA1[k][j][i] * LocalK11[k][j][i] * Localh1[k][j][i - 1]
                  + LocalA1[k][j][i - 1] * LocalK11[k][j][i - 1]
                      * Localh1[k][j][i]);
          LocalAKHx1p[k][j][i] = 2.0 * LocalA1[k][j][i] * LocalK11[k][j][i]
              * LocalA1[k][j][i + 1] * LocalK11[k][j][i + 1]
              / (LocalA1[k][j][i] * LocalK11[k][j][i] * Localh1[k][j][i + 1]
                  + LocalA1[k][j][i + 1] * LocalK11[k][j][i + 1]
                      * Localh1[k][j][i]);
          LocalAKHx2m[k][j][i] = 2.0 * LocalA2[k][j][i] * LocalK22[k][j][i]
              * LocalA2[k][j - 1][i] * LocalK22[k][j - 1][i]
              / (LocalA2[k][j][i] * LocalK22[k][j][i] * Localh2[k][j - 1][i]
                  + LocalA2[k][j - 1][i] * LocalK22[k][j - 1][i]
                      * Localh2[k][j][i]);
          LocalAKHx2p[k][j][i] = 2.0 * LocalA2[k][j][i] * LocalK22[k][j][i]
              * LocalA2[k][j + 1][i] * LocalK22[k][j + 1][i]
              / (LocalA2[k][j][i] * LocalK22[k][j][i] * Localh2[k][j + 1][i]
                  + LocalA2[k][j + 1][i] * LocalK22[k][j + 1][i]
                      * Localh2[k][j][i]);
          LocalAKHx3m[k][j][i] = 2.0 * LocalA3[k][j][i] * LocalK33[k][j][i]
              * LocalA3[k - 1][j][i] * LocalK33[k - 1][j][i]
              / (LocalA3[k][j][i] * LocalK33[k][j][i] * Localh3[k - 1][j][i]
                  + LocalA3[k - 1][j][i] * LocalK33[k - 1][j][i]
                      * Localh3[k][j][i]);
          LocalAKHx3p[k][j][i] = 2.0 * LocalA3[k][j][i] * LocalK33[k][j][i]
              * LocalA3[k + 1][j][i] * LocalK33[k + 1][j][i]
              / (LocalA3[k][j][i] * LocalK33[k][j][i] * Localh3[k + 1][j][i]
                  + LocalA3[k + 1][j][i] * LocalK33[k + 1][j][i]
                      * Localh3[k][j][i]);
        } else {
          LocalAKHx1m[k][j][i] = 0.0;
          LocalAKHx1p[k][j][i] = 0.0;
          LocalAKHx2m[k][j][i] = 0.0;
          LocalAKHx2p[k][j][i] = 0.0;
          LocalAKHx3m[k][j][i] = 0.0;
          LocalAKHx3p[k][j][i] = 0.0;
        }
      }
    }
  }

  /* Restore the new arrays to their rightful place */
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->AKHx1m, &LocalAKHx1m);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->AKHx2m, &LocalAKHx2m);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->AKHx3m, &LocalAKHx3m);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->AKHx1p, &LocalAKHx1p);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->AKHx2p, &LocalAKHx2p);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->AKHx3p, &LocalAKHx3p);CHKERRQ(ierr);

  /* Begin Assembly for vectors */
  ierr = VecAssemblyBegin(MySim->AKHx1m);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->AKHx2m);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->AKHx3m);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->AKHx1p);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->AKHx2p);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->AKHx3p);CHKERRQ(ierr);
  /* And end Assembly */
  ierr = VecAssemblyEnd(MySim->AKHx1m);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->AKHx2m);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->AKHx3m);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->AKHx1p);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->AKHx2p);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->AKHx3p);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DefiantComputeRhoAndMuAtFaces"
PetscErrorCode DefiantComputeRhoAndMuAtFaces(BlackOilReservoirSimulation* MySim) {
  PetscErrorCode ierr;
  PetscInt i, j, k, mx, my, mz, xm, ym, zm, xs, ys, zs;
  PetscScalar Beta;
  PetscScalar ***LocalFlowMask;
  PetscScalar ***LocalPhi; /* Porosity    */
  PetscScalar ***LocalRhoo, ***LocalRhow, ***LocalRhog; /* Densities   */
  PetscScalar ***LocalMuo, ***LocalMuw, ***LocalMug; /* Viscosities */
  /* Local Density at the faces */
  PetscScalar ***LocalRhoox1m, ***LocalRhoox1p, ***LocalRhowx1m,
      ***LocalRhowx1p, ***LocalRhogx1m, ***LocalRhogx1p;
  PetscScalar ***LocalRhoox2m, ***LocalRhoox2p, ***LocalRhowx2m,
      ***LocalRhowx2p, ***LocalRhogx2m, ***LocalRhogx2p;
  PetscScalar ***LocalRhoox3m, ***LocalRhoox3p, ***LocalRhowx3m,
      ***LocalRhowx3p, ***LocalRhogx3m, ***LocalRhogx3p;
  /* Local viscosity at the faces */
  PetscScalar ***LocalMuox1m, ***LocalMuox1p, ***LocalMuwx1m, ***LocalMuwx1p,
      ***LocalMugx1m, ***LocalMugx1p;
  PetscScalar ***LocalMuox2m, ***LocalMuox2p, ***LocalMuwx2m, ***LocalMuwx2p,
      ***LocalMugx2m, ***LocalMugx2p;
  PetscScalar ***LocalMuox3m, ***LocalMuox3p, ***LocalMuwx3m, ***LocalMuwx3p,
      ***LocalMugx3m, ***LocalMugx3p;
  /* Density By Viscosity at the faces */
  PetscScalar ***LocalRhoByMuox1m, ***LocalRhoByMuox1p, ***LocalRhoByMuwx1m,
      ***LocalRhoByMuwx1p, ***LocalRhoByMugx1m, ***LocalRhoByMugx1p;
  PetscScalar ***LocalRhoByMuox2m, ***LocalRhoByMuox2p, ***LocalRhoByMuwx2m,
      ***LocalRhoByMuwx2p, ***LocalRhoByMugx2m, ***LocalRhoByMugx2p;
  PetscScalar ***LocalRhoByMuox3m, ***LocalRhoByMuox3p, ***LocalRhoByMuwx3m,
      ***LocalRhoByMuwx3p, ***LocalRhoByMugx3m, ***LocalRhoByMugx3p;

  /* Vectors for ghosting */
  Vec vecLocalFlowMask;
  Vec vecLocalPhi; /* Porosity    */
  Vec vecLocalRhoo, vecLocalRhow, vecLocalRhog; /* Densities   */
  Vec vecLocalMuo, vecLocalMuw, vecLocalMug; /* Viscosities */

  PetscFunctionBegin;

  /* Get dimensions and extents of the local vectors */
  ierr = DAGetInfo(MySim->SimDA, 0, &mx, &my, &mz, 0, 0, 0, 0, 0, 0, 0);CHKERRQ(ierr);
  ierr = DAGetCorners(MySim->SimDA, &xs, &ys, &zs, &xm, &ym, &zm);CHKERRQ(ierr);

  ierr = DAGetLocalVector(MySim->SimDA, &vecLocalFlowMask);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGetLocalVector(MySim->SimDA, &vecLocalPhi);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGetLocalVector(MySim->SimDA, &vecLocalRhoo);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGetLocalVector(MySim->SimDA, &vecLocalRhow);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGetLocalVector(MySim->SimDA, &vecLocalRhog);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGetLocalVector(MySim->SimDA, &vecLocalMuo);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGetLocalVector(MySim->SimDA, &vecLocalMuw);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGetLocalVector(MySim->SimDA, &vecLocalMug);CHKERRQ(ierr);CHKMEMQ;

  ierr = DAGlobalToLocalBegin(MySim->SimDA, MySim-> FlowMask  ,INSERT_VALUES, vecLocalFlowMask  );CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGlobalToLocalEnd(MySim->SimDA, MySim-> FlowMask  ,INSERT_VALUES, vecLocalFlowMask  );CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGlobalToLocalBegin(MySim->SimDA, MySim-> Phi  ,INSERT_VALUES, vecLocalPhi  );CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGlobalToLocalEnd(MySim->SimDA, MySim-> Phi  ,INSERT_VALUES, vecLocalPhi  );CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGlobalToLocalBegin(MySim->SimDA, MySim-> Rhoo  ,INSERT_VALUES, vecLocalRhoo  );CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGlobalToLocalEnd(MySim->SimDA, MySim-> Rhoo  ,INSERT_VALUES, vecLocalRhoo  );CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGlobalToLocalBegin(MySim->SimDA, MySim-> Rhow  ,INSERT_VALUES, vecLocalRhow  );CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGlobalToLocalEnd(MySim->SimDA, MySim-> Rhow  ,INSERT_VALUES, vecLocalRhow  );CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGlobalToLocalBegin(MySim->SimDA, MySim-> Rhog  ,INSERT_VALUES, vecLocalRhog  );CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGlobalToLocalEnd(MySim->SimDA, MySim-> Rhog  ,INSERT_VALUES, vecLocalRhog  );CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGlobalToLocalBegin(MySim->SimDA, MySim-> Muo  ,INSERT_VALUES, vecLocalMuo  );CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGlobalToLocalEnd(MySim->SimDA, MySim-> Muo  ,INSERT_VALUES, vecLocalMuo  );CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGlobalToLocalBegin(MySim->SimDA, MySim-> Muw  ,INSERT_VALUES, vecLocalMuw  );CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGlobalToLocalEnd(MySim->SimDA, MySim-> Muw  ,INSERT_VALUES, vecLocalMuw  );CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGlobalToLocalBegin(MySim->SimDA, MySim-> Mug  ,INSERT_VALUES, vecLocalMug  );CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGlobalToLocalEnd(MySim->SimDA, MySim-> Mug  ,INSERT_VALUES, vecLocalMug  );CHKERRQ(ierr);CHKMEMQ;

  /* Grab the data for the flow field */
  ierr = DAVecGetArray(MySim->SimDA, vecLocalFlowMask, &LocalFlowMask);CHKERRQ(ierr);
  /* Grab the local data for porosity */
  ierr = DAVecGetArray(MySim->SimDA, vecLocalPhi, &LocalPhi);CHKERRQ(ierr);
  /* Grab the Densities */
  ierr = DAVecGetArray(MySim->SimDA, vecLocalRhoo, &LocalRhoo);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, vecLocalRhow, &LocalRhow);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, vecLocalRhog, &LocalRhog);CHKERRQ(ierr);
  /* Grab the viscosities */
  ierr = DAVecGetArray(MySim->SimDA, vecLocalMuo, &LocalMuo);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, vecLocalMuw, &LocalMuw);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, vecLocalMug, &LocalMug);CHKERRQ(ierr);


  /* Grab the local data for the face densities */
  ierr = DAVecGetArray(MySim->SimDA, MySim->Rhoox1m, &LocalRhoox1m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Rhoox1p, &LocalRhoox1p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Rhowx1m, &LocalRhowx1m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Rhowx1p, &LocalRhowx1p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Rhogx1m, &LocalRhogx1m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Rhogx1p, &LocalRhogx1p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Rhoox2m, &LocalRhoox2m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Rhoox2p, &LocalRhoox2p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Rhowx2m, &LocalRhowx2m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Rhowx2p, &LocalRhowx2p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Rhogx2m, &LocalRhogx2m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Rhogx2p, &LocalRhogx2p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Rhoox3m, &LocalRhoox3m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Rhoox3p, &LocalRhoox3p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Rhowx3m, &LocalRhowx3m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Rhowx3p, &LocalRhowx3p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Rhogx3m, &LocalRhogx3m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Rhogx3p, &LocalRhogx3p);CHKERRQ(ierr);
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
  /* Density divided by viscosity at the faces */
  ierr = DAVecGetArray(MySim->SimDA, MySim->RhoByMuox1m, &LocalRhoByMuox1m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->RhoByMuox1p, &LocalRhoByMuox1p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->RhoByMuwx1m, &LocalRhoByMuwx1m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->RhoByMuwx1p, &LocalRhoByMuwx1p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->RhoByMugx1m, &LocalRhoByMugx1m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->RhoByMugx1p, &LocalRhoByMugx1p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->RhoByMuox2m, &LocalRhoByMuox2m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->RhoByMuox2p, &LocalRhoByMuox2p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->RhoByMuwx2m, &LocalRhoByMuwx2m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->RhoByMuwx2p, &LocalRhoByMuwx2p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->RhoByMugx2m, &LocalRhoByMugx2m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->RhoByMugx2p, &LocalRhoByMugx2p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->RhoByMuox3m, &LocalRhoByMuox3m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->RhoByMuox3p, &LocalRhoByMuox3p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->RhoByMuwx3m, &LocalRhoByMuwx3m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->RhoByMuwx3p, &LocalRhoByMuwx3p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->RhoByMugx3m, &LocalRhoByMugx3m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->RhoByMugx3p, &LocalRhoByMugx3p);CHKERRQ(ierr);

  for (k = zs; k < zs + zm; k++) {
    for (j = ys; j < ys + ym; j++) {
      for (i = xs; i < xs + xm; i++) {
        if (i == 0 || j == 0 || k == 0 || i == mx - 1 || j == my - 1 || k == mz
            - 1) {
        } else if (ABS(LocalFlowMask[k][j][i]-FLUID_FLOW) < EPSILON) {
          Beta = LocalPhi[k][j][i] / LocalPhi[k][j][i - 1];
          /* Density at Face */
          LocalRhoox1m[k][j][i] = Beta * LocalRhoo[k][j][i] + (1.0 - Beta)
              * LocalRhoo[k][j][i - 1];
          LocalRhowx1m[k][j][i] = Beta * LocalRhow[k][j][i] + (1.0 - Beta)
              * LocalRhow[k][j][i - 1];
          LocalRhogx1m[k][j][i] = Beta * LocalRhog[k][j][i] + (1.0 - Beta)
              * LocalRhog[k][j][i - 1];
          /* Viscosity at Face */
          LocalMuox1m[k][j][i] = 1.0 / (Beta / LocalMuo[k][j][i] + (1.0 - Beta)
              / LocalMuo[k][j][i - 1]);
          LocalMuwx1m[k][j][i] = 1.0 / (Beta / LocalMuw[k][j][i] + (1.0 - Beta)
              / LocalMuw[k][j][i - 1]);
          LocalMugx1m[k][j][i] = 1.0 / (Beta / LocalMug[k][j][i] + (1.0 - Beta)
              / LocalMug[k][j][i - 1]);
          /* Density by viscosity at the face */
          LocalRhoByMuox1m[k][j][i] = Beta * LocalRhoo[k][j][i]
              / LocalMuo[k][j][i] + (1.0 - Beta) * LocalRhoo[k][j][i - 1]
              / LocalMuo[k][j][i - 1];
          LocalRhoByMuwx1m[k][j][i] = Beta * LocalRhow[k][j][i]
              / LocalMuw[k][j][i] + (1.0 - Beta) * LocalRhow[k][j][i - 1]
              / LocalMuw[k][j][i - 1];
          LocalRhoByMugx1m[k][j][i] = Beta * LocalRhog[k][j][i]
              / LocalMug[k][j][i] + (1.0 - Beta) * LocalRhog[k][j][i - 1]
              / LocalMug[k][j][i - 1];

          Beta = LocalPhi[k][j][i] / LocalPhi[k][j][i + 1];
          /* Density at Face */
          LocalRhoox1p[k][j][i] = Beta * LocalRhoo[k][j][i] + (1.0 - Beta)
              * LocalRhoo[k][j][i + 1];
          LocalRhowx1p[k][j][i] = Beta * LocalRhow[k][j][i] + (1.0 - Beta)
              * LocalRhow[k][j][i + 1];
          LocalRhogx1p[k][j][i] = Beta * LocalRhog[k][j][i] + (1.0 - Beta)
              * LocalRhog[k][j][i + 1];
          /* Viscosity at Face */
          LocalMuox1p[k][j][i] = 1.0 / (Beta / LocalMuo[k][j][i] + (1.0 - Beta)
              / LocalMuo[k][j][i + 1]);
          LocalMuwx1p[k][j][i] = 1.0 / (Beta / LocalMuw[k][j][i] + (1.0 - Beta)
              / LocalMuw[k][j][i + 1]);
          LocalMugx1p[k][j][i] = 1.0 / (Beta / LocalMug[k][j][i] + (1.0 - Beta)
              / LocalMug[k][j][i + 1]);
          /* Density by viscosity at the face */
          LocalRhoByMuox1p[k][j][i] = Beta * LocalRhoo[k][j][i]
              / LocalMuo[k][j][i] + (1.0 - Beta) * LocalRhoo[k][j][i + 1]
              / LocalMuo[k][j][i + 1];
          LocalRhoByMuwx1p[k][j][i] = Beta * LocalRhow[k][j][i]
              / LocalMuw[k][j][i] + (1.0 - Beta) * LocalRhow[k][j][i + 1]
              / LocalMuw[k][j][i + 1];
          LocalRhoByMugx1p[k][j][i] = Beta * LocalRhog[k][j][i]
              / LocalMug[k][j][i] + (1.0 - Beta) * LocalRhog[k][j][i + 1]
              / LocalMug[k][j][i + 1];

          Beta = LocalPhi[k][j][i] / LocalPhi[k][j - 1][i];
          /* Density at Face */
          LocalRhoox2m[k][j][i] = Beta * LocalRhoo[k][j][i] + (1.0 - Beta)
              * LocalRhoo[k][j - 1][i];
          LocalRhowx2m[k][j][i] = Beta * LocalRhow[k][j][i] + (1.0 - Beta)
              * LocalRhow[k][j - 1][i];
          LocalRhogx2m[k][j][i] = Beta * LocalRhog[k][j][i] + (1.0 - Beta)
              * LocalRhog[k][j - 1][i];
          /* Viscosity at Face */
          LocalMuox2m[k][j][i] = 1.0 / (Beta / LocalMuo[k][j][i] + (1.0 - Beta)
              / LocalMuo[k][j - 1][i]);
          LocalMuwx2m[k][j][i] = 1.0 / (Beta / LocalMuw[k][j][i] + (1.0 - Beta)
              / LocalMuw[k][j - 1][i]);
          LocalMugx2m[k][j][i] = 1.0 / (Beta / LocalMug[k][j][i] + (1.0 - Beta)
              / LocalMug[k][j - 1][i]);
          /* Density by viscosity at the face */
          LocalRhoByMuox2m[k][j][i] = Beta * LocalRhoo[k][j][i]
              / LocalMuo[k][j][i] + (1.0 - Beta) * LocalRhoo[k][j - 1][i]
              / LocalMuo[k][j - 1][i];
          LocalRhoByMuwx2m[k][j][i] = Beta * LocalRhow[k][j][i]
              / LocalMuw[k][j][i] + (1.0 - Beta) * LocalRhow[k][j - 1][i]
              / LocalMuw[k][j - 1][i];
          LocalRhoByMugx2m[k][j][i] = Beta * LocalRhog[k][j][i]
              / LocalMug[k][j][i] + (1.0 - Beta) * LocalRhog[k][j - 1][i]
              / LocalMug[k][j - 1][i];

          Beta = LocalPhi[k][j][i] / LocalPhi[k][j + 1][i];
          /* Density at Face */
          LocalRhoox2p[k][j][i] = Beta * LocalRhoo[k][j][i] + (1.0 - Beta)
              * LocalRhoo[k][j + 1][i];
          LocalRhowx2p[k][j][i] = Beta * LocalRhow[k][j][i] + (1.0 - Beta)
              * LocalRhow[k][j + 1][i];
          LocalRhogx2p[k][j][i] = Beta * LocalRhog[k][j][i] + (1.0 - Beta)
              * LocalRhog[k][j + 1][i];
          /* Viscosity at Face */
          LocalMuox2p[k][j][i] = 1.0 / (Beta / LocalMuo[k][j][i] + (1.0 - Beta)
              / LocalMuo[k][j + 1][i]);
          LocalMuwx2p[k][j][i] = 1.0 / (Beta / LocalMuw[k][j][i] + (1.0 - Beta)
              / LocalMuw[k][j + 1][i]);
          LocalMugx2p[k][j][i] = 1.0 / (Beta / LocalMug[k][j][i] + (1.0 - Beta)
              / LocalMug[k][j + 1][i]);
          /* Density by viscosity at the face */
          LocalRhoByMuox2p[k][j][i] = Beta * LocalRhoo[k][j][i]
              / LocalMuo[k][j][i] + (1.0 - Beta) * LocalRhoo[k][j + 1][i]
              / LocalMuo[k][j + 1][i];
          LocalRhoByMuwx2p[k][j][i] = Beta * LocalRhow[k][j][i]
              / LocalMuw[k][j][i] + (1.0 - Beta) * LocalRhow[k][j + 1][i]
              / LocalMuw[k][j + 1][i];
          LocalRhoByMugx2p[k][j][i] = Beta * LocalRhog[k][j][i]
              / LocalMug[k][j][i] + (1.0 - Beta) * LocalRhog[k][j + 1][i]
              / LocalMug[k][j + 1][i];

          Beta = LocalPhi[k][j][i] / LocalPhi[k - 1][j][i];
          /* Density at Face */
          LocalRhoox3m[k][j][i] = Beta * LocalRhoo[k][j][i] + (1.0 - Beta)
              * LocalRhoo[k - 1][j][i];
          LocalRhowx3m[k][j][i] = Beta * LocalRhow[k][j][i] + (1.0 - Beta)
              * LocalRhow[k - 1][j][i];
          LocalRhogx3m[k][j][i] = Beta * LocalRhog[k][j][i] + (1.0 - Beta)
              * LocalRhog[k - 1][j][i];
          /* Viscosity at Face */
          LocalMuox3m[k][j][i] = 1.0 / (Beta / LocalMuo[k][j][i] + (1.0 - Beta)
              / LocalMuo[k - 1][j][i]);
          LocalMuwx3m[k][j][i] = 1.0 / (Beta / LocalMuw[k][j][i] + (1.0 - Beta)
              / LocalMuw[k - 1][j][i]);
          LocalMugx3m[k][j][i] = 1.0 / (Beta / LocalMug[k][j][i] + (1.0 - Beta)
              / LocalMug[k - 1][j][i]);
          /* Density by viscosity at the face */
          LocalRhoByMuox3m[k][j][i] = Beta * LocalRhoo[k][j][i]
              / LocalMuo[k][j][i] + (1.0 - Beta) * LocalRhoo[k - 1][j][i]
              / LocalMuo[k - 1][j][i];
          LocalRhoByMuwx3m[k][j][i] = Beta * LocalRhow[k][j][i]
              / LocalMuw[k][j][i] + (1.0 - Beta) * LocalRhow[k - 1][j][i]
              / LocalMuw[k - 1][j][i];
          LocalRhoByMugx3m[k][j][i] = Beta * LocalRhog[k][j][i]
              / LocalMug[k][j][i] + (1.0 - Beta) * LocalRhog[k - 1][j][i]
              / LocalMug[k - 1][j][i];

          Beta = LocalPhi[k][j][i] / LocalPhi[k + 1][j][i];
          /* Density at Face */
          LocalRhoox3p[k][j][i] = Beta * LocalRhoo[k][j][i] + (1.0 - Beta)
              * LocalRhoo[k + 1][j][i];
          LocalRhowx3p[k][j][i] = Beta * LocalRhow[k][j][i] + (1.0 - Beta)
              * LocalRhow[k + 1][j][i];
          LocalRhogx3p[k][j][i] = Beta * LocalRhog[k][j][i] + (1.0 - Beta)
              * LocalRhog[k + 1][j][i];
          /* Viscosity at Face */
          LocalMuox3p[k][j][i] = 1.0 / (Beta / LocalMuo[k][j][i] + (1.0 - Beta)
              / LocalMuo[k + 1][j][i]);
          LocalMuwx3p[k][j][i] = 1.0 / (Beta / LocalMuw[k][j][i] + (1.0 - Beta)
              / LocalMuw[k + 1][j][i]);
          LocalMugx3p[k][j][i] = 1.0 / (Beta / LocalMug[k][j][i] + (1.0 - Beta)
              / LocalMug[k + 1][j][i]);
          /* Density by viscosity at the face */
          LocalRhoByMuox3p[k][j][i] = Beta * LocalRhoo[k][j][i]
              / LocalMuo[k][j][i] + (1.0 - Beta) * LocalRhoo[k + 1][j][i]
              / LocalMuo[k + 1][j][i];
          LocalRhoByMuwx3p[k][j][i] = Beta * LocalRhow[k][j][i]
              / LocalMuw[k][j][i] + (1.0 - Beta) * LocalRhow[k + 1][j][i]
              / LocalMuw[k + 1][j][i];
          LocalRhoByMugx3p[k][j][i] = Beta * LocalRhog[k][j][i]
              / LocalMug[k][j][i] + (1.0 - Beta) * LocalRhog[k + 1][j][i]
              / LocalMug[k + 1][j][i];
        }
      }
    }
  }

  /* Restore the new arrays to their rightful place for face densities */
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Rhoox1m, &LocalRhoox1m);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Rhoox1p, &LocalRhoox1p);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Rhowx1m, &LocalRhowx1m);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Rhowx1p, &LocalRhowx1p);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Rhogx1m, &LocalRhogx1m);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Rhogx1p, &LocalRhogx1p);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Rhoox2m, &LocalRhoox2m);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Rhoox2p, &LocalRhoox2p);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Rhowx2m, &LocalRhowx2m);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Rhowx2p, &LocalRhowx2p);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Rhogx2m, &LocalRhogx2m);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Rhogx2p, &LocalRhogx2p);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Rhoox3m, &LocalRhoox3m);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Rhoox3p, &LocalRhoox3p);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Rhowx3m, &LocalRhowx3m);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Rhowx3p, &LocalRhowx3p);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Rhogx3m, &LocalRhogx3m);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Rhogx3p, &LocalRhogx3p);CHKERRQ(ierr);
  /* Restore the new arrays to their rightful place for face viscosities */
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Muox1m, &LocalMuox1m);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Muox1p, &LocalMuox1p);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Muwx1m, &LocalMuwx1m);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Muwx1p, &LocalMuwx1p);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Mugx1m, &LocalMugx1m);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Mugx1p, &LocalMugx1p);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Muox2m, &LocalMuox2m);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Muox2p, &LocalMuox2p);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Muwx2m, &LocalMuwx2m);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Muwx2p, &LocalMuwx2p);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Mugx2m, &LocalMugx2m);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Mugx2p, &LocalMugx2p);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Muox3m, &LocalMuox3m);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Muox3p, &LocalMuox3p);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Muwx3m, &LocalMuwx3m);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Muwx3p, &LocalMuwx3p);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Mugx3m, &LocalMugx3m);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Mugx3p, &LocalMugx3p);CHKERRQ(ierr);
  /* Restore density divided by viscosity at the faces */
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->RhoByMuox1m, &LocalRhoByMuox1m);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->RhoByMuox1p, &LocalRhoByMuox1p);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->RhoByMuwx1m, &LocalRhoByMuwx1m);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->RhoByMuwx1p, &LocalRhoByMuwx1p);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->RhoByMugx1m, &LocalRhoByMugx1m);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->RhoByMugx1p, &LocalRhoByMugx1p);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->RhoByMuox2m, &LocalRhoByMuox2m);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->RhoByMuox2p, &LocalRhoByMuox2p);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->RhoByMuwx2m, &LocalRhoByMuwx2m);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->RhoByMuwx2p, &LocalRhoByMuwx2p);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->RhoByMugx2m, &LocalRhoByMugx2m);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->RhoByMugx2p, &LocalRhoByMugx2p);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->RhoByMuox3m, &LocalRhoByMuox3m);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->RhoByMuox3p, &LocalRhoByMuox3p);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->RhoByMuwx3m, &LocalRhoByMuwx3m);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->RhoByMuwx3p, &LocalRhoByMuwx3p);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->RhoByMugx3m, &LocalRhoByMugx3m);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->RhoByMugx3p, &LocalRhoByMugx3p);CHKERRQ(ierr);

  /* Begin Assembly for vectors */
  ierr = VecAssemblyBegin(MySim->Rhoox1m);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->Rhoox1p);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->Rhowx1m);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->Rhowx1p);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->Rhogx1m);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->Rhogx1p);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->Rhoox2m);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->Rhoox2p);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->Rhowx2m);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->Rhowx2p);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->Rhogx2m);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->Rhogx2p);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->Rhoox3m);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->Rhoox3p);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->Rhowx3m);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->Rhowx3p);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->Rhogx3m);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->Rhogx3p);CHKERRQ(ierr);

  ierr = VecAssemblyBegin(MySim->Muox1m);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->Muox1p);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->Muwx1m);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->Muwx1p);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->Mugx1m);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->Mugx1p);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->Muox2m);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->Muox2p);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->Muwx2m);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->Muwx2p);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->Mugx2m);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->Mugx2p);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->Muox3m);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->Muox3p);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->Muwx3m);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->Muwx3p);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->Mugx3m);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->Mugx3p);CHKERRQ(ierr);

  ierr = VecAssemblyBegin(MySim->RhoByMuox1m);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->RhoByMuox1p);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->RhoByMuwx1m);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->RhoByMuwx1p);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->RhoByMugx1m);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->RhoByMugx1p);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->RhoByMuox2m);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->RhoByMuox2p);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->RhoByMuwx2m);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->RhoByMuwx2p);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->RhoByMugx2m);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->RhoByMugx2p);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->RhoByMuox3m);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->RhoByMuox3p);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->RhoByMuwx3m);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->RhoByMuwx3p);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->RhoByMugx3m);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->RhoByMugx3p);CHKERRQ(ierr);

  /* And end Assembly */
  ierr = VecAssemblyEnd(MySim->Rhoox1m);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->Rhoox1p);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->Rhowx1m);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->Rhowx1p);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->Rhogx1m);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->Rhogx1p);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->Rhoox2m);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->Rhoox2p);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->Rhowx2m);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->Rhowx2p);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->Rhogx2m);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->Rhogx2p);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->Rhoox3m);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->Rhoox3p);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->Rhowx3m);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->Rhowx3p);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->Rhogx3m);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->Rhogx3p);CHKERRQ(ierr);

  ierr = VecAssemblyEnd(MySim->Muox1m);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->Muox1p);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->Muwx1m);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->Muwx1p);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->Mugx1m);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->Mugx1p);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->Muox2m);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->Muox2p);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->Muwx2m);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->Muwx2p);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->Mugx2m);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->Mugx2p);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->Muox3m);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->Muox3p);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->Muwx3m);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->Muwx3p);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->Mugx3m);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->Mugx3p);CHKERRQ(ierr);

  ierr = VecAssemblyEnd(MySim->RhoByMuox1m);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->RhoByMuox1p);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->RhoByMuwx1m);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->RhoByMuwx1p);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->RhoByMugx1m);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->RhoByMugx1p);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->RhoByMuox2m);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->RhoByMuox2p);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->RhoByMuwx2m);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->RhoByMuwx2p);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->RhoByMugx2m);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->RhoByMugx2p);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->RhoByMuox3m);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->RhoByMuox3p);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->RhoByMuwx3m);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->RhoByMuwx3p);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->RhoByMugx3m);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->RhoByMugx3p);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DefiantComputeRelativePermsAtFaces"
PetscErrorCode DefiantComputeRelativePermsAtFaces(BlackOilReservoirSimulation* MySim) {
  PetscErrorCode ierr;
  PetscInt i, j, k, mx, my, mz, xm, ym, zm, xs, ys, zs;
  PetscScalar DeltaPot, Beta, BetaPrime; /* Difference in potential, weighting factors */
  PetscScalar ***LocalFlowMask;
  PetscScalar ***Localh1, ***Localh2, ***Localh3;
  PetscScalar ***Localx1, ***Localx2, ***Localx3; /* Geometry    */
  PetscScalar ***LocalPo, ***LocalPw, ***LocalPg; /* Pressure    */
  PetscScalar ***LocalRhoo, ***LocalRhow, ***LocalRhog; /* Densities   */
  /* Relative permeabilities at the cell center */
  PetscScalar ***LocalKro, ***LocalKrw, ***LocalKrg;
  /* Relative Permeabilities at the faces */
  PetscScalar ***LocalRelPermox1m, ***LocalRelPermox1p, ***LocalRelPermox2m,
      ***LocalRelPermox2p, ***LocalRelPermox3m, ***LocalRelPermox3p;
  PetscScalar ***LocalRelPermwx1m, ***LocalRelPermwx1p, ***LocalRelPermwx2m,
      ***LocalRelPermwx2p, ***LocalRelPermwx3m, ***LocalRelPermwx3p;
  PetscScalar ***LocalRelPermgx1m, ***LocalRelPermgx1p, ***LocalRelPermgx2m,
      ***LocalRelPermgx2p, ***LocalRelPermgx3m, ***LocalRelPermgx3p;
  /* Local Density at the faces */
  PetscScalar ***LocalRhoox1m, ***LocalRhoox1p, ***LocalRhowx1m,
      ***LocalRhowx1p, ***LocalRhogx1m, ***LocalRhogx1p;
  PetscScalar ***LocalRhoox2m, ***LocalRhoox2p, ***LocalRhowx2m,
      ***LocalRhowx2p, ***LocalRhogx2m, ***LocalRhogx2p;
  PetscScalar ***LocalRhoox3m, ***LocalRhoox3p, ***LocalRhowx3m,
      ***LocalRhowx3p, ***LocalRhogx3m, ***LocalRhogx3p;

  Vec vecLocalFlowMask;
  Vec vecLocalh1, vecLocalh2, vecLocalh3;
  Vec vecLocalx1, vecLocalx2, vecLocalx3; /* Geometry    */
  Vec vecLocalPo, vecLocalPw, vecLocalPg; /* Pressure    */
  Vec vecLocalRhoo, vecLocalRhow, vecLocalRhog; /* Densities   */
  /* Relative permeabilities at the cell center */
  Vec vecLocalKro, vecLocalKrw, vecLocalKrg;

  PetscFunctionBegin;

  /* Get dimensions and extents of the local vectors */
  ierr = DAGetInfo(MySim->SimDA, 0, &mx, &my, &mz, 0, 0, 0, 0, 0, 0, 0);CHKERRQ(ierr);
  ierr = DAGetCorners(MySim->SimDA, &xs, &ys, &zs, &xm, &ym, &zm);CHKERRQ(ierr);

  ierr = DAGetLocalVector(MySim->SimDA, &vecLocalFlowMask);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGetLocalVector(MySim->SimDA, &vecLocalh1);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGetLocalVector(MySim->SimDA, &vecLocalh2);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGetLocalVector(MySim->SimDA, &vecLocalh3);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGetLocalVector(MySim->SimDA, &vecLocalx1);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGetLocalVector(MySim->SimDA, &vecLocalx2);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGetLocalVector(MySim->SimDA, &vecLocalx3);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGetLocalVector(MySim->SimDA, &vecLocalPo);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGetLocalVector(MySim->SimDA, &vecLocalPw);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGetLocalVector(MySim->SimDA, &vecLocalPg);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGetLocalVector(MySim->SimDA, &vecLocalRhoo);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGetLocalVector(MySim->SimDA, &vecLocalRhow);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGetLocalVector(MySim->SimDA, &vecLocalRhog);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGetLocalVector(MySim->SimDA, &vecLocalKro);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGetLocalVector(MySim->SimDA, &vecLocalKrw);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGetLocalVector(MySim->SimDA, &vecLocalKrg);CHKERRQ(ierr);CHKMEMQ;

  ierr = DAGlobalToLocalBegin(MySim->SimDA,MySim->FlowMask,INSERT_VALUES,vecLocalFlowMask);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGlobalToLocalEnd(MySim->SimDA,MySim->FlowMask,INSERT_VALUES,vecLocalFlowMask);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGlobalToLocalBegin(MySim->SimDA,MySim->h1,INSERT_VALUES,vecLocalh1);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGlobalToLocalEnd(MySim->SimDA,MySim->h1,INSERT_VALUES,vecLocalh1);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGlobalToLocalBegin(MySim->SimDA,MySim->h2,INSERT_VALUES,vecLocalh2);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGlobalToLocalEnd(MySim->SimDA,MySim->h2,INSERT_VALUES,vecLocalh2);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGlobalToLocalBegin(MySim->SimDA,MySim->h3,INSERT_VALUES,vecLocalh3);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGlobalToLocalEnd(MySim->SimDA,MySim->h3,INSERT_VALUES,vecLocalh3);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGlobalToLocalBegin(MySim->SimDA,MySim->x1,INSERT_VALUES,vecLocalx1);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGlobalToLocalEnd(MySim->SimDA,MySim->x1,INSERT_VALUES,vecLocalx1);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGlobalToLocalBegin(MySim->SimDA,MySim->x2,INSERT_VALUES,vecLocalx2);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGlobalToLocalEnd(MySim->SimDA,MySim->x2,INSERT_VALUES,vecLocalx2);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGlobalToLocalBegin(MySim->SimDA,MySim->x3,INSERT_VALUES,vecLocalx3);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGlobalToLocalEnd(MySim->SimDA,MySim->x3,INSERT_VALUES,vecLocalx3);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGlobalToLocalBegin(MySim->SimDA,MySim->Po,INSERT_VALUES,vecLocalPo);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGlobalToLocalEnd(MySim->SimDA,MySim->Po,INSERT_VALUES,vecLocalPo);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGlobalToLocalBegin(MySim->SimDA,MySim->Pw,INSERT_VALUES,vecLocalPw);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGlobalToLocalEnd(MySim->SimDA,MySim->Pw,INSERT_VALUES,vecLocalPw);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGlobalToLocalBegin(MySim->SimDA,MySim->Pg,INSERT_VALUES,vecLocalPg);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGlobalToLocalEnd(MySim->SimDA,MySim->Pg,INSERT_VALUES,vecLocalPg);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGlobalToLocalBegin(MySim->SimDA,MySim->Rhoo,INSERT_VALUES,vecLocalRhoo);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGlobalToLocalEnd(MySim->SimDA,MySim->Rhoo,INSERT_VALUES,vecLocalRhoo);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGlobalToLocalBegin(MySim->SimDA,MySim->Rhow,INSERT_VALUES,vecLocalRhow);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGlobalToLocalEnd(MySim->SimDA,MySim->Rhow,INSERT_VALUES,vecLocalRhow);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGlobalToLocalBegin(MySim->SimDA,MySim->Rhog,INSERT_VALUES,vecLocalRhog);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGlobalToLocalEnd(MySim->SimDA,MySim->Rhog,INSERT_VALUES,vecLocalRhog);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGlobalToLocalBegin(MySim->SimDA,MySim->Kro,INSERT_VALUES,vecLocalKro);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGlobalToLocalEnd(MySim->SimDA,MySim->Kro,INSERT_VALUES,vecLocalKro);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGlobalToLocalBegin(MySim->SimDA,MySim->Krw,INSERT_VALUES,vecLocalKrw);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGlobalToLocalEnd(MySim->SimDA,MySim->Krw,INSERT_VALUES,vecLocalKrw);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGlobalToLocalBegin(MySim->SimDA,MySim->Krg,INSERT_VALUES,vecLocalKrg);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGlobalToLocalEnd(MySim->SimDA,MySim->Krg,INSERT_VALUES,vecLocalKrg);CHKERRQ(ierr);CHKMEMQ;

  /* Grab the data for the flow field */
  ierr = DAVecGetArray(MySim->SimDA, vecLocalFlowMask, &LocalFlowMask);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, vecLocalh1, &Localh1);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, vecLocalh2, &Localh2);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, vecLocalh3, &Localh3);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, vecLocalx1, &Localx1);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, vecLocalx2, &Localx2);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, vecLocalx3, &Localx3);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, vecLocalRhoo, &LocalRhoo);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, vecLocalRhow, &LocalRhow);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, vecLocalRhog, &LocalRhog);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, vecLocalPo, &LocalPo);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, vecLocalPw, &LocalPw);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, vecLocalPg, &LocalPg);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, vecLocalKro, &LocalKro);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, vecLocalKrw, &LocalKrw);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, vecLocalKrg, &LocalKrg);CHKERRQ(ierr);

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
  /* Grab the local data for the face densities */
  ierr = DAVecGetArray(MySim->SimDA, MySim->Rhoox1m, &LocalRhoox1m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Rhoox1p, &LocalRhoox1p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Rhowx1m, &LocalRhowx1m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Rhowx1p, &LocalRhowx1p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Rhogx1m, &LocalRhogx1m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Rhogx1p, &LocalRhogx1p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Rhoox2m, &LocalRhoox2m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Rhoox2p, &LocalRhoox2p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Rhowx2m, &LocalRhowx2m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Rhowx2p, &LocalRhowx2p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Rhogx2m, &LocalRhogx2m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Rhogx2p, &LocalRhogx2p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Rhoox3m, &LocalRhoox3m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Rhoox3p, &LocalRhoox3p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Rhowx3m, &LocalRhowx3m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Rhowx3p, &LocalRhowx3p);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Rhogx3m, &LocalRhogx3m);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Rhogx3p, &LocalRhogx3p);CHKERRQ(ierr);


  for (k = zs; k < zs + zm; k++) {
    for (j = ys; j < ys + ym; j++) {
      for (i = xs; i < xs + xm; i++) {
        if (i == 0 || j == 0 || k == 0 || i == mx - 1 || j == my - 1 || k == mz  - 1) {
        } else if (ABS(LocalFlowMask[k][j][i]-FLUID_FLOW) < EPSILON) {
          /* in the x1 direction */
          if ( i == 1 || mx -i == 2 ){  /* if I am the point right next to the physical boundary */
            if ( i  == 1 ) {
              LocalRelPermox1m[k][j][i] = 0.0;
              LocalRelPermwx1m[k][j][i] = 0.0;
              LocalRelPermgx1m[k][j][i] = 0.0;
              if ( xs + xm - i > 2 )
              {
                DeltaPot = LocalPo[k][j][i + 1] - LocalPo[k][j][i] - MySim->GravAcc * LocalRhoox1p[k][j][i] * (Localx3[k][j][i + 1] - Localx3[k][j][i]);
                if ( DeltaPot < 0.0 )
                  LocalRelPermox1p[k][j][i] = LocalKro[k][j][i];
                else
                  LocalRelPermox1p[k][j][i] = LocalKro[k][j][i + 1];

                DeltaPot = LocalPw[k][j][i + 1] - LocalPw[k][j][i] - MySim->GravAcc * LocalRhowx1p[k][j][i] * (Localx3[k][j][i + 1] - Localx3[k][j][i]);
                if ( DeltaPot < 0.0 )
                  LocalRelPermwx1p[k][j][i] = LocalKrw[k][j][i];
                else
                  LocalRelPermwx1p[k][j][i] = LocalKrw[k][j][i + 1];

                DeltaPot = LocalPg[k][j][i + 1] - LocalPg[k][j][i] - MySim->GravAcc * LocalRhogx1p[k][j][i] * (Localx3[k][j][i + 1] - Localx3[k][j][i]);
                if ( DeltaPot < 0.0 )
                  LocalRelPermgx1p[k][j][i] = LocalKrg[k][j][i];
                else
                  LocalRelPermgx1p[k][j][i] = LocalKrg[k][j][i + 1];
              }
            }
            if ( mx - i == 2 ) {
              LocalRelPermox1p[k][j][i] = 0.0;
              LocalRelPermwx1p[k][j][i] = 0.0;
              LocalRelPermgx1p[k][j][i] = 0.0;
              if ( i - xs > 1 ) {
                DeltaPot = LocalPo[k][j][i] - LocalPo[k][j][i - 1] - MySim->GravAcc * LocalRhoox1m[k][j][i] * (Localx3[k][j][i] - Localx3[k][j][i - 1]);
                if ( DeltaPot < 0.0 )
                  LocalRelPermox1m[k][j][i] = LocalKro[k][j][i - 1];
                else
                  LocalRelPermox1m[k][j][i] = LocalKro[k][j][i];

                DeltaPot = LocalPw[k][j][i] - LocalPw[k][j][i - 1] - MySim->GravAcc * LocalRhowx1m[k][j][i] * (Localx3[k][j][i] - Localx3[k][j][i - 1]);
                if ( DeltaPot < 0.0 )
                  LocalRelPermwx1m[k][j][i] = LocalKrw[k][j][i - 1];
                else
                  LocalRelPermwx1m[k][j][i] = LocalKrw[k][j][i];

                DeltaPot = LocalPg[k][j][i] - LocalPg[k][j][i - 1] - MySim->GravAcc * LocalRhogx1m[k][j][i] * (Localx3[k][j][i] - Localx3[k][j][i - 1]);
                if ( DeltaPot < 0.0 )
                  LocalRelPermgx1m[k][j][i] = LocalKrg[k][j][i - 1];
                else
                  LocalRelPermgx1m[k][j][i] = LocalKrg[k][j][i];
              }
            }
          }
          else if ( i == 2 || i == mx - 3 ) { /* else if I am the point next to next to a physical boundary */
            DeltaPot = LocalPo[k][j][i] - LocalPo[k][j][i - 1] - MySim->GravAcc * LocalRhoox1m[k][j][i] * (Localx3[k][j][i] - Localx3[k][j][i - 1]);
            if ( DeltaPot < 0.0 )
              LocalRelPermox1m[k][j][i] = LocalKro[k][j][i - 1];
            else
              LocalRelPermox1m[k][j][i] = LocalKro[k][j][i];

            DeltaPot = LocalPw[k][j][i] - LocalPw[k][j][i - 1] - MySim->GravAcc * LocalRhowx1m[k][j][i] * (Localx3[k][j][i] - Localx3[k][j][i - 1]);
            if ( DeltaPot < 0.0 )
              LocalRelPermwx1m[k][j][i] = LocalKrw[k][j][i - 1];
            else
              LocalRelPermwx1m[k][j][i] = LocalKrw[k][j][i];

            DeltaPot = LocalPg[k][j][i] - LocalPg[k][j][i - 1] - MySim->GravAcc * LocalRhogx1m[k][j][i] * (Localx3[k][j][i] - Localx3[k][j][i - 1]);
            if ( DeltaPot < 0.0 )
              LocalRelPermgx1m[k][j][i] = LocalKrg[k][j][i - 1];
            else
              LocalRelPermgx1m[k][j][i] = LocalKrg[k][j][i];

            DeltaPot = LocalPo[k][j][i + 1] - LocalPo[k][j][i] - MySim->GravAcc * LocalRhoox1p[k][j][i] * (Localx3[k][j][i + 1] - Localx3[k][j][i]);
            if ( DeltaPot < 0.0 )
              LocalRelPermox1p[k][j][i] = LocalKro[k][j][i];
            else
              LocalRelPermox1p[k][j][i] = LocalKro[k][j][i + 1];

            DeltaPot = LocalPw[k][j][i + 1] - LocalPw[k][j][i] - MySim->GravAcc * LocalRhowx1p[k][j][i] * (Localx3[k][j][i + 1] - Localx3[k][j][i]);
            if ( DeltaPot < 0.0 )
              LocalRelPermwx1p[k][j][i] = LocalKrw[k][j][i];
            else
              LocalRelPermwx1p[k][j][i] = LocalKrw[k][j][i + 1];

            DeltaPot = LocalPg[k][j][i + 1] - LocalPg[k][j][i] - MySim->GravAcc * LocalRhogx1p[k][j][i] * (Localx3[k][j][i + 1] - Localx3[k][j][i]);
            if ( DeltaPot < 0.0 )
              LocalRelPermgx1p[k][j][i] = LocalKrg[k][j][i];
            else
              LocalRelPermgx1p[k][j][i] = LocalKrg[k][j][i + 1];
          }
          else if (( i == xs && xs != 0 ) || ( i == xs + xm - 1 && xs + xm - 1 != mx - 1 )){ /* else if I am the first and last point next to a ghost zone */
            DeltaPot = LocalPo[k][j][i] - LocalPo[k][j][i - 1] - MySim->GravAcc * LocalRhoox1m[k][j][i] * (Localx3[k][j][i] - Localx3[k][j][i - 1]);
            if ( DeltaPot < 0.0 )
              LocalRelPermox1m[k][j][i] = LocalKro[k][j][i - 1];
            else
              LocalRelPermox1m[k][j][i] = LocalKro[k][j][i];

            DeltaPot = LocalPw[k][j][i] - LocalPw[k][j][i - 1] - MySim->GravAcc * LocalRhowx1m[k][j][i] * (Localx3[k][j][i] - Localx3[k][j][i - 1]);
            if ( DeltaPot < 0.0 )
              LocalRelPermwx1m[k][j][i] = LocalKrw[k][j][i - 1];
            else
              LocalRelPermwx1m[k][j][i] = LocalKrw[k][j][i];

            DeltaPot = LocalPg[k][j][i] - LocalPg[k][j][i - 1] - MySim->GravAcc * LocalRhogx1m[k][j][i] * (Localx3[k][j][i] - Localx3[k][j][i - 1]);
            if ( DeltaPot < 0.0 )
              LocalRelPermgx1m[k][j][i] = LocalKrg[k][j][i - 1];
            else
              LocalRelPermgx1m[k][j][i] = LocalKrg[k][j][i];

            DeltaPot = LocalPo[k][j][i + 1] - LocalPo[k][j][i] - MySim->GravAcc * LocalRhoox1p[k][j][i] * (Localx3[k][j][i + 1] - Localx3[k][j][i]);
            if ( DeltaPot < 0.0 )
              LocalRelPermox1p[k][j][i] = LocalKro[k][j][i];
            else
              LocalRelPermox1p[k][j][i] = LocalKro[k][j][i + 1];

            DeltaPot = LocalPw[k][j][i + 1] - LocalPw[k][j][i] - MySim->GravAcc * LocalRhowx1p[k][j][i] * (Localx3[k][j][i + 1] - Localx3[k][j][i]);
            if ( DeltaPot < 0.0 )
              LocalRelPermwx1p[k][j][i] = LocalKrw[k][j][i];
            else
              LocalRelPermwx1p[k][j][i] = LocalKrw[k][j][i + 1];

            DeltaPot = LocalPg[k][j][i + 1] - LocalPg[k][j][i] - MySim->GravAcc * LocalRhogx1p[k][j][i] * (Localx3[k][j][i + 1] - Localx3[k][j][i]);
            if ( DeltaPot < 0.0 )
              LocalRelPermgx1p[k][j][i] = LocalKrg[k][j][i];
            else
              LocalRelPermgx1p[k][j][i] = LocalKrg[k][j][i + 1];
          }
          else if ( i - xs > 1 && xs + xm - i > 2 ) {
            Beta = 0.5 * Localh1[k][j][i - 1] / (Localx1[k][j][i - 1] - Localx1[k][j][i - 2]);
            BetaPrime = 0.5 * Localh1[k][j][i] / (Localx1[k][j][i + 1] - Localx1[k][j][i]);

            DeltaPot = LocalPo[k][j][i] - LocalPo[k][j][i - 1] - MySim->GravAcc * LocalRhoox1m[k][j][i] * (Localx3[k][j][i] - Localx3[k][j][i - 1]);
            if ( DeltaPot < 0.0 )
              LocalRelPermox1m[k][j][i] = (1 + Beta) * LocalKro[k][j][i - 1] - Beta * LocalKro[k][j][i - 2];
            else
              LocalRelPermox1m[k][j][i] = (1 + BetaPrime) * LocalKro[k][j][i] - BetaPrime * LocalKro[k][j][i + 1];

            DeltaPot = LocalPw[k][j][i] - LocalPw[k][j][i - 1] - MySim->GravAcc * LocalRhowx1m[k][j][i] * (Localx3[k][j][i] - Localx3[k][j][i - 1]);
            if ( DeltaPot < 0.0 )
              LocalRelPermwx1m[k][j][i] = (1 + Beta) * LocalKrw[k][j][i - 1] - Beta * LocalKrw[k][j][i - 2];
            else
              LocalRelPermwx1m[k][j][i] = (1 + BetaPrime) * LocalKrw[k][j][i] - BetaPrime * LocalKrw[k][j][i + 1];

            DeltaPot = LocalPg[k][j][i] - LocalPg[k][j][i - 1] - MySim->GravAcc * LocalRhogx1m[k][j][i] * (Localx3[k][j][i] - Localx3[k][j][i - 1]);
            if ( DeltaPot < 0.0 )
              LocalRelPermgx1m[k][j][i] = (1 + Beta) * LocalKrg[k][j][i - 1] - Beta * LocalKrg[k][j][i - 2];
            else
              LocalRelPermgx1m[k][j][i] = (1 + BetaPrime) * LocalKrg[k][j][i] - BetaPrime * LocalKrg[k][j][i + 1];

            Beta = 0.5 * Localh1[k][j][i] / (Localx1[k][j][i] - Localx1[k][j][i - 1]);
            BetaPrime = 0.5 * Localh1[k][j][i + 1] / (Localx1[k][j][i + 2] - Localx1[k][j][i + 1]);

            DeltaPot = LocalPo[k][j][i + 1] - LocalPo[k][j][i] - MySim->GravAcc * LocalRhoox1p[k][j][i] * (Localx3[k][j][i + 1] - Localx3[k][j][i]);
            if ( DeltaPot < 0.0 )
              LocalRelPermox1p[k][j][i] = (1 + Beta) * LocalKro[k][j][i] - Beta * LocalKro[k][j][i - 1];
            else
              LocalRelPermox1p[k][j][i] = (1 + BetaPrime) * LocalKro[k][j][i + 1] - BetaPrime * LocalKro[k][j][i + 2];

            DeltaPot = LocalPw[k][j][i + 1] - LocalPw[k][j][i] - MySim->GravAcc * LocalRhowx1p[k][j][i] * (Localx3[k][j][i + 1] - Localx3[k][j][i]);
            if ( DeltaPot < 0.0 )
              LocalRelPermwx1p[k][j][i] = (1 + Beta) * LocalKrw[k][j][i] - Beta * LocalKrw[k][j][i - 1];
            else
              LocalRelPermwx1p[k][j][i] = (1 + BetaPrime) * LocalKrw[k][j][i + 1] - BetaPrime * LocalKrw[k][j][i + 2];

            DeltaPot = LocalPg[k][j][i + 1] - LocalPg[k][j][i] - MySim->GravAcc * LocalRhogx1p[k][j][i] * (Localx3[k][j][i + 1] - Localx3[k][j][i]);
            if ( DeltaPot < 0.0 )
              LocalRelPermgx1p[k][j][i] = (1 + Beta) * LocalKrg[k][j][i] - Beta * LocalKrg[k][j][i - 1];
            else
              LocalRelPermgx1p[k][j][i] = (1 + BetaPrime) * LocalKrg[k][j][i + 1] - BetaPrime * LocalKrg[k][j][i + 2];
          }
          else if (( i - xs > 0 && xs != 0 ) || ( xs + xm - i > 1 && xs + xm != mx )) {
            Beta = 0.5 * Localh1[k][j][i - 1] / (Localx1[k][j][i - 1] - Localx1[k][j][i - 2]);
            BetaPrime = 0.5 * Localh1[k][j][i] / (Localx1[k][j][i + 1] - Localx1[k][j][i]);

            DeltaPot = LocalPo[k][j][i] - LocalPo[k][j][i - 1] - MySim->GravAcc * LocalRhoox1m[k][j][i] * (Localx3[k][j][i] - Localx3[k][j][i - 1]);
            if ( DeltaPot < 0.0 )
              LocalRelPermox1m[k][j][i] = (1 + Beta) * LocalKro[k][j][i - 1] - Beta * LocalKro[k][j][i - 2];
            else
              LocalRelPermox1m[k][j][i] = (1 + BetaPrime) * LocalKro[k][j][i] - BetaPrime * LocalKro[k][j][i + 1];

            DeltaPot = LocalPw[k][j][i] - LocalPw[k][j][i - 1] - MySim->GravAcc * LocalRhowx1m[k][j][i] * (Localx3[k][j][i] - Localx3[k][j][i - 1]);
            if ( DeltaPot < 0.0 )
              LocalRelPermwx1m[k][j][i] = (1 + Beta) * LocalKrw[k][j][i - 1] - Beta * LocalKrw[k][j][i - 2];
            else
              LocalRelPermwx1m[k][j][i] = (1 + BetaPrime) * LocalKrw[k][j][i] - BetaPrime * LocalKrw[k][j][i + 1];

            DeltaPot = LocalPg[k][j][i] - LocalPg[k][j][i - 1] - MySim->GravAcc * LocalRhogx1m[k][j][i] * (Localx3[k][j][i] - Localx3[k][j][i - 1]);
            if ( DeltaPot < 0.0 )
              LocalRelPermgx1m[k][j][i] = (1 + Beta) * LocalKrg[k][j][i - 1] - Beta * LocalKrg[k][j][i - 2];
            else
              LocalRelPermgx1m[k][j][i] = (1 + BetaPrime) * LocalKrg[k][j][i] - BetaPrime * LocalKrg[k][j][i + 1];

            Beta = 0.5 * Localh1[k][j][i] / (Localx1[k][j][i] - Localx1[k][j][i - 1]);
            BetaPrime = 0.5 * Localh1[k][j][i + 1] / (Localx1[k][j][i + 2] - Localx1[k][j][i + 1]);

            DeltaPot = LocalPo[k][j][i + 1] - LocalPo[k][j][i] - MySim->GravAcc * LocalRhoox1p[k][j][i] * (Localx3[k][j][i + 1] - Localx3[k][j][i]);
            if ( DeltaPot < 0.0 )
              LocalRelPermox1p[k][j][i] = (1 + Beta) * LocalKro[k][j][i] - Beta * LocalKro[k][j][i - 1];
            else
              LocalRelPermox1p[k][j][i] = (1 + BetaPrime) * LocalKro[k][j][i + 1] - BetaPrime * LocalKro[k][j][i + 2];

            DeltaPot = LocalPw[k][j][i + 1] - LocalPw[k][j][i] - MySim->GravAcc * LocalRhowx1p[k][j][i] * (Localx3[k][j][i + 1] - Localx3[k][j][i]);
            if ( DeltaPot < 0.0 )
              LocalRelPermwx1p[k][j][i] = (1 + Beta) * LocalKrw[k][j][i] - Beta * LocalKrw[k][j][i - 1];
            else
              LocalRelPermwx1p[k][j][i] = (1 + BetaPrime) * LocalKrw[k][j][i + 1] - BetaPrime * LocalKrw[k][j][i + 2];

            DeltaPot = LocalPg[k][j][i + 1] - LocalPg[k][j][i] - MySim->GravAcc * LocalRhogx1p[k][j][i] * (Localx3[k][j][i + 1] - Localx3[k][j][i]);
            if ( DeltaPot < 0.0 )
              LocalRelPermgx1p[k][j][i] = (1 + Beta) * LocalKrg[k][j][i] - Beta * LocalKrg[k][j][i - 1];
            else
              LocalRelPermgx1p[k][j][i] = (1 + BetaPrime) * LocalKrg[k][j][i + 1] - BetaPrime * LocalKrg[k][j][i + 2];
          }

          /* in the x2 direction */
          if ( j == 1 || my -j == 2 ){  /* if I am the point right next to the physical boundary */
            if ( j  == 1 ) {
              LocalRelPermox2m[k][j][i] = 0.0;
              LocalRelPermwx2m[k][j][i] = 0.0;
              LocalRelPermgx2m[k][j][i] = 0.0;
              if ( ys + ym - j > 2 )
              {
                DeltaPot = LocalPo[k][j + 1][i] - LocalPo[k][j][i] - MySim->GravAcc * LocalRhoox2p[k][j][i] * (Localx3[k][j + 1][i] - Localx3[k][j][i]);
                if ( DeltaPot < 0.0 )
                  LocalRelPermox2p[k][j][i] = LocalKro[k][j][i];
                else
                  LocalRelPermox2p[k][j][i] = LocalKro[k][j + 1][i];

                DeltaPot = LocalPw[k][j + 1][i] - LocalPw[k][j][i] - MySim->GravAcc * LocalRhowx2p[k][j][i] * (Localx3[k][j + 1][i] - Localx3[k][j][i]);
                if ( DeltaPot < 0.0 )
                  LocalRelPermwx2p[k][j][i] = LocalKrw[k][j][i];
                else
                  LocalRelPermwx2p[k][j][i] = LocalKrw[k][j + 1][i];

                DeltaPot = LocalPg[k][j + 1][i] - LocalPg[k][j][i] - MySim->GravAcc * LocalRhogx2p[k][j][i] * (Localx3[k][j + 1][i] - Localx3[k][j][i]);
                if ( DeltaPot < 0.0 )
                  LocalRelPermgx2p[k][j][i] = LocalKrg[k][j][i];
                else
                  LocalRelPermgx2p[k][j][i] = LocalKrg[k][j + 1][i];
              }
            }
            if ( my - j == 2 ) {
              LocalRelPermox2p[k][j][i] = 0.0;
              LocalRelPermwx2p[k][j][i] = 0.0;
              LocalRelPermgx2p[k][j][i] = 0.0;
              if ( j - ys > 1 ) {
                DeltaPot = LocalPo[k][j][i] - LocalPo[k][j - 1][i] - MySim->GravAcc * LocalRhoox2m[k][j][i] * (Localx3[k][j][i] - Localx3[k][j - 1][i]);
                if ( DeltaPot < 0.0 )
                  LocalRelPermox2m[k][j][i] = LocalKro[k][j - 1][i];
                else
                  LocalRelPermox2m[k][j][i] = LocalKro[k][j][i];

                DeltaPot = LocalPw[k][j][i] - LocalPw[k][j - 1][i] - MySim->GravAcc * LocalRhowx2m[k][j][i] * (Localx3[k][j][i] - Localx3[k][j - 1][i]);
                if ( DeltaPot < 0.0 )
                  LocalRelPermwx2m[k][j][i] = LocalKrw[k][j - 1][i];
                else
                  LocalRelPermwx2m[k][j][i] = LocalKrw[k][j][i];

                DeltaPot = LocalPg[k][j][i] - LocalPg[k][j - 1][i] - MySim->GravAcc * LocalRhogx2m[k][j][i] * (Localx3[k][j][i] - Localx3[k][j - 1][i]);
                if ( DeltaPot < 0.0 )
                  LocalRelPermgx2m[k][j][i] = LocalKrg[k][j - 1][i];
                else
                  LocalRelPermgx2m[k][j][i] = LocalKrg[k][j][i];
              }
            }
          }
          else if ( j == 2 || j == my - 3 ) { /* else if I am the point next to next to a physical boundary */
            DeltaPot = LocalPo[k][j][i] - LocalPo[k][j - 1][i] - MySim->GravAcc * LocalRhoox2m[k][j][i] * (Localx3[k][j][i] - Localx3[k][j - 1][i]);
            if ( DeltaPot < 0.0 )
              LocalRelPermox2m[k][j][i] = LocalKro[k][j - 1][i];
            else
              LocalRelPermox2m[k][j][i] = LocalKro[k][j][i];

            DeltaPot = LocalPw[k][j][i] - LocalPw[k][j - 1][i] - MySim->GravAcc * LocalRhowx2m[k][j][i] * (Localx3[k][j][i] - Localx3[k][j - 1][i]);
            if ( DeltaPot < 0.0 )
              LocalRelPermwx2m[k][j][i] = LocalKrw[k][j - 1][i];
            else
              LocalRelPermwx2m[k][j][i] = LocalKrw[k][j][i];

            DeltaPot = LocalPg[k][j][i] - LocalPg[k][j - 1][i] - MySim->GravAcc * LocalRhogx2m[k][j][i] * (Localx3[k][j][i] - Localx3[k][j - 1][i]);
            if ( DeltaPot < 0.0 )
              LocalRelPermgx2m[k][j][i] = LocalKrg[k][j - 1][i];
            else
              LocalRelPermgx2m[k][j][i] = LocalKrg[k][j][i];

            DeltaPot = LocalPo[k][j + 1][i] - LocalPo[k][j][i] - MySim->GravAcc * LocalRhoox2p[k][j][i] * (Localx3[k][j + 1][i] - Localx3[k][j][i]);
            if ( DeltaPot < 0.0 )
              LocalRelPermox2p[k][j][i] = LocalKro[k][j][i];
            else
              LocalRelPermox2p[k][j][i] = LocalKro[k][j + 1][i];

            DeltaPot = LocalPw[k][j + 1][i] - LocalPw[k][j][i] - MySim->GravAcc * LocalRhowx2p[k][j][i] * (Localx3[k][j + 1][i] - Localx3[k][j][i]);
            if ( DeltaPot < 0.0 )
              LocalRelPermwx2p[k][j][i] = LocalKrw[k][j][i];
            else
              LocalRelPermwx2p[k][j][i] = LocalKrw[k][j + 1][i];

            DeltaPot = LocalPg[k][j + 1][i] - LocalPg[k][j][i] - MySim->GravAcc * LocalRhogx2p[k][j][i] * (Localx3[k][j + 1][i] - Localx3[k][j][i]);
            if ( DeltaPot < 0.0 )
              LocalRelPermgx2p[k][j][i] = LocalKrg[k][j][i];
            else
              LocalRelPermgx2p[k][j][i] = LocalKrg[k][j + 1][i];
          }
          else if (( j == ys && ys != 0 ) || ( j == ys + ym - 1 && ys + ym - 1 != my - 1 )){ /* else if I am the first and last point next to a ghost zone */
            DeltaPot = LocalPo[k][j][i] - LocalPo[k][j - 1][i] - MySim->GravAcc * LocalRhoox2m[k][j][i] * (Localx3[k][j][i] - Localx3[k][j - 1][i]);
            if ( DeltaPot < 0.0 )
              LocalRelPermox2m[k][j][i] = LocalKro[k][j - 1][i];
            else
              LocalRelPermox2m[k][j][i] = LocalKro[k][j][i];

            DeltaPot = LocalPw[k][j][i] - LocalPw[k][j - 1][i] - MySim->GravAcc * LocalRhowx2m[k][j][i] * (Localx3[k][j][i] - Localx3[k][j - 1][i]);
            if ( DeltaPot < 0.0 )
              LocalRelPermwx2m[k][j][i] = LocalKrw[k][j - 1][i];
            else
              LocalRelPermwx2m[k][j][i] = LocalKrw[k][j][i];

            DeltaPot = LocalPg[k][j][i] - LocalPg[k][j - 1][i] - MySim->GravAcc * LocalRhogx2m[k][j][i] * (Localx3[k][j][i] - Localx3[k][j - 1][i]);
            if ( DeltaPot < 0.0 )
              LocalRelPermgx2m[k][j][i] = LocalKrg[k][j - 1][i];
            else
              LocalRelPermgx2m[k][j][i] = LocalKrg[k][j][i];

            DeltaPot = LocalPo[k][j + 1][i] - LocalPo[k][j][i] - MySim->GravAcc * LocalRhoox2p[k][j][i] * (Localx3[k][j + 1][i] - Localx3[k][j][i]);
            if ( DeltaPot < 0.0 )
              LocalRelPermox2p[k][j][i] = LocalKro[k][j][i];
            else
              LocalRelPermox2p[k][j][i] = LocalKro[k][j + 1][i];

            DeltaPot = LocalPw[k][j + 1][i] - LocalPw[k][j][i] - MySim->GravAcc * LocalRhowx2p[k][j][i] * (Localx3[k][j + 1][i] - Localx3[k][j][i]);
            if ( DeltaPot < 0.0 )
              LocalRelPermwx2p[k][j][i] = LocalKrw[k][j][i];
            else
              LocalRelPermwx2p[k][j][i] = LocalKrw[k][j + 1][i];

            DeltaPot = LocalPg[k][j + 1][i] - LocalPg[k][j][i] - MySim->GravAcc * LocalRhogx2p[k][j][i] * (Localx3[k][j + 1][i] - Localx3[k][j][i]);
            if ( DeltaPot < 0.0 )
              LocalRelPermgx2p[k][j][i] = LocalKrg[k][j][i];
            else
              LocalRelPermgx2p[k][j][i] = LocalKrg[k][j + 1][i];
          }
          else if ( j - ys > 1 && ys + ym - j > 2 ) {
            Beta = 0.5 * Localh2[k][j - 1][i] / (Localx2[k][j - 1][i] - Localx2[k][j - 2][i]);
            BetaPrime = 0.5 * Localh2[k][j][i] / (Localx2[k][j + 1][i] - Localx2[k][j][i]);

            DeltaPot = LocalPo[k][j][i] - LocalPo[k][j - 1][i] - MySim->GravAcc * LocalRhoox2m[k][j][i] * (Localx3[k][j][i] - Localx3[k][j - 1][i]);
            if ( DeltaPot < 0.0 )
              LocalRelPermox2m[k][j][i] = (1 + Beta) * LocalKro[k][j - 1][i] - Beta * LocalKro[k][j - 2][i];
            else
              LocalRelPermox2m[k][j][i] = (1 + BetaPrime) * LocalKro[k][j][i] - BetaPrime * LocalKro[k][j + 1][i];

            DeltaPot = LocalPw[k][j][i] - LocalPw[k][j - 1][i] - MySim->GravAcc * LocalRhowx2m[k][j][i] * (Localx3[k][j][i] - Localx3[k][j - 1][i]);
            if ( DeltaPot < 0.0 )
              LocalRelPermwx2m[k][j][i] = (1 + Beta) * LocalKrw[k][j - 1][i] - Beta * LocalKrw[k][j - 2][i];
            else
              LocalRelPermwx2m[k][j][i] = (1 + BetaPrime) * LocalKrw[k][j][i] - BetaPrime * LocalKrw[k][j + 1][i];

            DeltaPot = LocalPg[k][j][i] - LocalPg[k][j - 1][i] - MySim->GravAcc * LocalRhogx2m[k][j][i] * (Localx3[k][j][i] - Localx3[k][j - 1][i]);
            if ( DeltaPot < 0.0 )
              LocalRelPermgx2m[k][j][i] = (1 + Beta) * LocalKrg[k][j - 1][i] - Beta * LocalKrg[k][j - 2][i];
            else
              LocalRelPermgx2m[k][j][i] = (1 + BetaPrime) * LocalKrg[k][j][i] - BetaPrime * LocalKrg[k][j + 1][i];

            Beta = 0.5 * Localh2[k][j][i] / (Localx2[k][j][i] - Localx2[k][j - 1][i]);
            BetaPrime = 0.5 * Localh2[k][j + 1][i] / (Localx2[k][j + 2][i] - Localx2[k][j + 1][i]);

            DeltaPot = LocalPo[k][j + 1][i] - LocalPo[k][j][i] - MySim->GravAcc * LocalRhoox2p[k][j][i] * (Localx3[k][j + 1][i] - Localx3[k][j][i]);
            if ( DeltaPot < 0.0 )
              LocalRelPermox2p[k][j][i] = (1 + Beta) * LocalKro[k][j][i] - Beta * LocalKro[k][j - 1][i];
            else
              LocalRelPermox2p[k][j][i] = (1 + BetaPrime) * LocalKro[k][j + 1][i] - BetaPrime * LocalKro[k][j + 2][i];

            DeltaPot = LocalPw[k][j + 1][i] - LocalPw[k][j][i] - MySim->GravAcc * LocalRhowx2p[k][j][i] * (Localx3[k][j + 1][i] - Localx3[k][j][i]);
            if ( DeltaPot < 0.0 )
              LocalRelPermwx2p[k][j][i] = (1 + Beta) * LocalKrw[k][j][i] - Beta * LocalKrw[k][j - 1][i];
            else
              LocalRelPermwx2p[k][j][i] = (1 + BetaPrime) * LocalKrw[k][j + 1][i] - BetaPrime * LocalKrw[k][j + 2][i];

            DeltaPot = LocalPg[k][j + 1][i] - LocalPg[k][j][i] - MySim->GravAcc * LocalRhogx2p[k][j][i] * (Localx3[k][j + 1][i] - Localx3[k][j][i]);
            if ( DeltaPot < 0.0 )
              LocalRelPermgx2p[k][j][i] = (1 + Beta) * LocalKrg[k][j][i] - Beta * LocalKrg[k][j - 1][i];
            else
              LocalRelPermgx2p[k][j][i] = (1 + BetaPrime) * LocalKrg[k][j + 1][i] - BetaPrime * LocalKrg[k][j + 2][i];
          }
          else if (( j - ys > 0 && ys != 0 ) || ( ys + ym - j > 1 && ys + ym != my )) {
            Beta = 0.5 * Localh2[k][j - 1][i] / (Localx2[k][j - 1][i] - Localx2[k][j - 2][i]);
            BetaPrime = 0.5 * Localh2[k][j][i] / (Localx2[k][j + 1][i] - Localx2[k][j][i]);

            DeltaPot = LocalPo[k][j][i] - LocalPo[k][j - 1][i] - MySim->GravAcc * LocalRhoox2m[k][j][i] * (Localx3[k][j][i] - Localx3[k][j - 1][i]);
            if ( DeltaPot < 0.0 )
              LocalRelPermox2m[k][j][i] = (1 + Beta) * LocalKro[k][j - 1][i] - Beta * LocalKro[k][j - 2][i];
            else
              LocalRelPermox2m[k][j][i] = (1 + BetaPrime) * LocalKro[k][j][i] - BetaPrime * LocalKro[k][j + 1][i];

            DeltaPot = LocalPw[k][j][i] - LocalPw[k][j - 1][i] - MySim->GravAcc * LocalRhowx2m[k][j][i] * (Localx3[k][j][i] - Localx3[k][j - 1][i]);
            if ( DeltaPot < 0.0 )
              LocalRelPermwx2m[k][j][i] = (1 + Beta) * LocalKrw[k][j - 1][i] - Beta * LocalKrw[k][j - 2][i];
            else
              LocalRelPermwx2m[k][j][i] = (1 + BetaPrime) * LocalKrw[k][j][i] - BetaPrime * LocalKrw[k][j + 1][i];

            DeltaPot = LocalPg[k][j][i] - LocalPg[k][j - 1][i] - MySim->GravAcc * LocalRhogx2m[k][j][i] * (Localx3[k][j][i] - Localx3[k][j - 1][i]);
            if ( DeltaPot < 0.0 )
              LocalRelPermgx2m[k][j][i] = (1 + Beta) * LocalKrg[k][j - 1][i] - Beta * LocalKrg[k][j - 2][i];
            else
              LocalRelPermgx2m[k][j][i] = (1 + BetaPrime) * LocalKrg[k][j][i] - BetaPrime * LocalKrg[k][j + 1][i];

            Beta = 0.5 * Localh2[k][j][i] / (Localx2[k][j][i] - Localx2[k][j - 1][i]);
            BetaPrime = 0.5 * Localh2[k][j + 1][i] / (Localx2[k][j + 2][i] - Localx2[k][j + 1][i]);

            DeltaPot = LocalPo[k][j + 1][i] - LocalPo[k][j][i] - MySim->GravAcc * LocalRhoox2p[k][j][i] * (Localx3[k][j + 1][i] - Localx3[k][j][i]);
            if ( DeltaPot < 0.0 )
              LocalRelPermox2p[k][j][i] = (1 + Beta) * LocalKro[k][j][i] - Beta * LocalKro[k][j - 1][i];
            else
              LocalRelPermox2p[k][j][i] = (1 + BetaPrime) * LocalKro[k][j + 1][i] - BetaPrime * LocalKro[k][j + 2][i];

            DeltaPot = LocalPw[k][j + 1][i] - LocalPw[k][j][i] - MySim->GravAcc * LocalRhowx2p[k][j][i] * (Localx3[k][j + 1][i] - Localx3[k][j][i]);
            if ( DeltaPot < 0.0 )
              LocalRelPermwx2p[k][j][i] = (1 + Beta) * LocalKrw[k][j][i] - Beta * LocalKrw[k][j - 1][i];
            else
              LocalRelPermwx2p[k][j][i] = (1 + BetaPrime) * LocalKrw[k][j + 1][i] - BetaPrime * LocalKrw[k][j + 2][i];

            DeltaPot = LocalPg[k][j + 1][i] - LocalPg[k][j][i] - MySim->GravAcc * LocalRhogx2p[k][j][i] * (Localx3[k][j + 1][i] - Localx3[k][j][i]);
            if ( DeltaPot < 0.0 )
              LocalRelPermgx2p[k][j][i] = (1 + Beta) * LocalKrg[k][j][i] - Beta * LocalKrg[k][j - 1][i];
            else
              LocalRelPermgx2p[k][j][i] = (1 + BetaPrime) * LocalKrg[k][j + 1][i] - BetaPrime * LocalKrg[k][j + 2][i];
          }

          /* in the x3 direction */
          if ( k == 1 || mz -k == 2 ){  /* if I am the point right next to the physical boundary */
            if ( k  == 1 ) {
              LocalRelPermox3m[k][j][i] = 0.0;
              LocalRelPermwx3m[k][j][i] = 0.0;
              LocalRelPermgx3m[k][j][i] = 0.0;
              if ( zs + zm - k > 2 )
              {
                DeltaPot = LocalPo[k + 1][j][i] - LocalPo[k][j][i] - MySim->GravAcc * LocalRhoox3p[k][j][i] * (Localx3[k + 1][j][i] - Localx3[k][j][i]);
                if ( DeltaPot < 0.0 )
                  LocalRelPermox3p[k][j][i] = LocalKro[k][j][i];
                else
                  LocalRelPermox3p[k][j][i] = LocalKro[k + 1][j][i];

                DeltaPot = LocalPw[k + 1][j][i] - LocalPw[k][j][i] - MySim->GravAcc * LocalRhowx3p[k][j][i] * (Localx3[k + 1][j][i] - Localx3[k][j][i]);
                if ( DeltaPot < 0.0 )
                  LocalRelPermwx3p[k][j][i] = LocalKrw[k][j][i];
                else
                  LocalRelPermwx3p[k][j][i] = LocalKrw[k + 1][j][i];

                DeltaPot = LocalPg[k + 1][j][i] - LocalPg[k][j][i] - MySim->GravAcc * LocalRhogx3p[k][j][i] * (Localx3[k + 1][j][i] - Localx3[k][j][i]);
                if ( DeltaPot < 0.0 )
                  LocalRelPermgx3p[k][j][i] = LocalKrg[k][j][i];
                else
                  LocalRelPermgx3p[k][j][i] = LocalKrg[k + 1][j][i];
              }
            }
            if ( mz - k == 2 ) {
              LocalRelPermox3p[k][j][i] = 0.0;
              LocalRelPermwx3p[k][j][i] = 0.0;
              LocalRelPermgx3p[k][j][i] = 0.0;
              if ( k - zs > 1 ) {
                DeltaPot = LocalPo[k][j][i] - LocalPo[k - 1][j][i] - MySim->GravAcc * LocalRhoox3m[k][j][i] * (Localx3[k][j][i] - Localx3[k - 1][j][i]);
                if ( DeltaPot < 0.0 )
                  LocalRelPermox3m[k][j][i] = LocalKro[k - 1][j][i];
                else
                  LocalRelPermox3m[k][j][i] = LocalKro[k][j][i];

                DeltaPot = LocalPw[k][j][i] - LocalPw[k - 1][j][i] - MySim->GravAcc * LocalRhowx3m[k][j][i] * (Localx3[k][j][i] - Localx3[k - 1][j][i]);
                if ( DeltaPot < 0.0 )
                  LocalRelPermwx3m[k][j][i] = LocalKrw[k - 1][j][i];
                else
                  LocalRelPermwx3m[k][j][i] = LocalKrw[k][j][i];

                DeltaPot = LocalPg[k][j][i] - LocalPg[k - 1][j][i] - MySim->GravAcc * LocalRhogx3m[k][j][i] * (Localx3[k][j][i] - Localx3[k - 1][j][i]);
                if ( DeltaPot < 0.0 )
                  LocalRelPermgx3m[k][j][i] = LocalKrg[k - 1][j][i];
                else
                  LocalRelPermgx3m[k][j][i] = LocalKrg[k][j][i];
              }
            }
          }
          else if ( k == 2 || k == mz - 3 ) { /* else if I am the point next to next to a physical boundary */
            DeltaPot = LocalPo[k][j][i] - LocalPo[k - 1][j][i] - MySim->GravAcc * LocalRhoox3m[k][j][i] * (Localx3[k][j][i] - Localx3[k - 1][j][i]);
            if ( DeltaPot < 0.0 )
              LocalRelPermox3m[k][j][i] = LocalKro[k - 1][j][i];
            else
              LocalRelPermox3m[k][j][i] = LocalKro[k][j][i];

            DeltaPot = LocalPw[k][j][i] - LocalPw[k - 1][j][i] - MySim->GravAcc * LocalRhowx3m[k][j][i] * (Localx3[k][j][i] - Localx3[k - 1][j][i]);
            if ( DeltaPot < 0.0 )
              LocalRelPermwx3m[k][j][i] = LocalKrw[k - 1][j][i];
            else
              LocalRelPermwx3m[k][j][i] = LocalKrw[k][j][i];

            DeltaPot = LocalPg[k][j][i] - LocalPg[k - 1][j][i] - MySim->GravAcc * LocalRhogx3m[k][j][i] * (Localx3[k][j][i] - Localx3[k - 1][j][i]);
            if ( DeltaPot < 0.0 )
              LocalRelPermgx3m[k][j][i] = LocalKrg[k - 1][j][i];
            else
              LocalRelPermgx3m[k][j][i] = LocalKrg[k][j][i];

            DeltaPot = LocalPo[k + 1][j][i] - LocalPo[k][j][i] - MySim->GravAcc * LocalRhoox3p[k][j][i] * (Localx3[k + 1][j][i] - Localx3[k][j][i]);
            if ( DeltaPot < 0.0 )
              LocalRelPermox3p[k][j][i] = LocalKro[k][j][i];
            else
              LocalRelPermox3p[k][j][i] = LocalKro[k + 1][j][i];

            DeltaPot = LocalPw[k + 1][j][i] - LocalPw[k][j][i] - MySim->GravAcc * LocalRhowx3p[k][j][i] * (Localx3[k + 1][j][i] - Localx3[k][j][i]);
            if ( DeltaPot < 0.0 )
              LocalRelPermwx3p[k][j][i] = LocalKrw[k][j][i];
            else
              LocalRelPermwx3p[k][j][i] = LocalKrw[k + 1][j][i];

            DeltaPot = LocalPg[k + 1][j][i] - LocalPg[k][j][i] - MySim->GravAcc * LocalRhogx3p[k][j][i] * (Localx3[k + 1][j][i] - Localx3[k][j][i]);
            if ( DeltaPot < 0.0 )
              LocalRelPermgx3p[k][j][i] = LocalKrg[k][j][i];
            else
              LocalRelPermgx3p[k][j][i] = LocalKrg[k + 1][j][i];
          }
          else if (( k == zs && zs != 0 ) || ( k == zs + zm - 1 && zs + zm - 1 != mz - 1 )){ /* else if I am the first and last point next to a ghost zone */
            DeltaPot = LocalPo[k][j][i] - LocalPo[k - 1][j][i] - MySim->GravAcc * LocalRhoox3m[k][j][i] * (Localx3[k][j][i] - Localx3[k - 1][j][i]);
            if ( DeltaPot < 0.0 )
              LocalRelPermox3m[k][j][i] = LocalKro[k - 1][j][i];
            else
              LocalRelPermox3m[k][j][i] = LocalKro[k][j][i];

            DeltaPot = LocalPw[k][j][i] - LocalPw[k - 1][j][i] - MySim->GravAcc * LocalRhowx3m[k][j][i] * (Localx3[k][j][i] - Localx3[k - 1][j][i]);
            if ( DeltaPot < 0.0 )
              LocalRelPermwx3m[k][j][i] = LocalKrw[k - 1][j][i];
            else
              LocalRelPermwx3m[k][j][i] = LocalKrw[k][j][i];

            DeltaPot = LocalPg[k][j][i] - LocalPg[k - 1][j][i] - MySim->GravAcc * LocalRhogx3m[k][j][i] * (Localx3[k][j][i] - Localx3[k - 1][j][i]);
            if ( DeltaPot < 0.0 )
              LocalRelPermgx3m[k][j][i] = LocalKrg[k - 1][j][i];
            else
              LocalRelPermgx3m[k][j][i] = LocalKrg[k][j][i];

            DeltaPot = LocalPo[k + 1][j][i] - LocalPo[k][j][i] - MySim->GravAcc * LocalRhoox3p[k][j][i] * (Localx3[k + 1][j][i] - Localx3[k][j][i]);
            if ( DeltaPot < 0.0 )
              LocalRelPermox3p[k][j][i] = LocalKro[k][j][i];
            else
              LocalRelPermox3p[k][j][i] = LocalKro[k + 1][j][i];

            DeltaPot = LocalPw[k + 1][j][i] - LocalPw[k][j][i] - MySim->GravAcc * LocalRhowx3p[k][j][i] * (Localx3[k + 1][j][i] - Localx3[k][j][i]);
            if ( DeltaPot < 0.0 )
              LocalRelPermwx3p[k][j][i] = LocalKrw[k][j][i];
            else
              LocalRelPermwx3p[k][j][i] = LocalKrw[k + 1][j][i];

            DeltaPot = LocalPg[k + 1][j][i] - LocalPg[k][j][i] - MySim->GravAcc * LocalRhogx3p[k][j][i] * (Localx3[k + 1][j][i] - Localx3[k][j][i]);
            if ( DeltaPot < 0.0 )
              LocalRelPermgx3p[k][j][i] = LocalKrg[k][j][i];
            else
              LocalRelPermgx3p[k][j][i] = LocalKrg[k + 1][j][i];
          }
          else if ( k - zs > 1 && zs + zm - k > 2 ) {
            Beta = 0.5 * Localh3[k - 1][j][i] / (Localx3[k - 1][j][i] - Localx3[k - 2][j][i]);
            BetaPrime = 0.5 * Localh3[k][j][i] / (Localx3[k + 1][j][i] - Localx3[k][j][i]);

            DeltaPot = LocalPo[k][j][i] - LocalPo[k - 1][j][i] - MySim->GravAcc * LocalRhoox3m[k][j][i] * (Localx3[k][j][i] - Localx3[k - 1][j][i]);
            if ( DeltaPot < 0.0 )
              LocalRelPermox3m[k][j][i] = (1 + Beta) * LocalKro[k - 1][j][i] - Beta * LocalKro[k - 2][j][i];
            else
              LocalRelPermox3m[k][j][i] = (1 + BetaPrime) * LocalKro[k][j][i] - BetaPrime * LocalKro[k + 1][j][i];

            DeltaPot = LocalPw[k][j][i] - LocalPw[k - 1][j][i] - MySim->GravAcc * LocalRhowx3m[k][j][i] * (Localx3[k][j][i] - Localx3[k - 1][j][i]);
            if ( DeltaPot < 0.0 )
              LocalRelPermwx3m[k][j][i] = (1 + Beta) * LocalKrw[k - 1][j][i] - Beta * LocalKrw[k - 2][j][i];
            else
              LocalRelPermwx3m[k][j][i] = (1 + BetaPrime) * LocalKrw[k][j][i] - BetaPrime * LocalKrw[k + 1][j][i];

            DeltaPot = LocalPg[k][j][i] - LocalPg[k - 1][j][i] - MySim->GravAcc * LocalRhogx3m[k][j][i] * (Localx3[k][j][i] - Localx3[k - 1][j][i]);
            if ( DeltaPot < 0.0 )
              LocalRelPermgx3m[k][j][i] = (1 + Beta) * LocalKrg[k - 1][j][i] - Beta * LocalKrg[k - 2][j][i];
            else
              LocalRelPermgx3m[k][j][i] = (1 + BetaPrime) * LocalKrg[k][j][i] - BetaPrime * LocalKrg[k + 1][j][i];

            Beta = 0.5 * Localh3[k][j][i] / (Localx3[k][j][i] - Localx3[k - 1][j][i]);
            BetaPrime = 0.5 * Localh3[k + 1][j][i] / (Localx3[k + 2][j][i] - Localx3[k + 1][j][i]);

            DeltaPot = LocalPo[k + 1][j][i] - LocalPo[k][j][i] - MySim->GravAcc * LocalRhoox3p[k][j][i] * (Localx3[k + 1][j][i] - Localx3[k][j][i]);
            if ( DeltaPot < 0.0 )
              LocalRelPermox3p[k][j][i] = (1 + Beta) * LocalKro[k][j][i] - Beta * LocalKro[k - 1][j][i];
            else
              LocalRelPermox3p[k][j][i] = (1 + BetaPrime) * LocalKro[k + 1][j][i] - BetaPrime * LocalKro[k + 2][j][i];

            DeltaPot = LocalPw[k + 1][j][i] - LocalPw[k][j][i] - MySim->GravAcc * LocalRhowx3p[k][j][i] * (Localx3[k + 1][j][i] - Localx3[k][j][i]);
            if ( DeltaPot < 0.0 )
              LocalRelPermwx3p[k][j][i] = (1 + Beta) * LocalKrw[k][j][i] - Beta * LocalKrw[k - 1][j][i];
            else
              LocalRelPermwx3p[k][j][i] = (1 + BetaPrime) * LocalKrw[k + 1][j][i] - BetaPrime * LocalKrw[k + 2][j][i];

            DeltaPot = LocalPg[k + 1][j][i] - LocalPg[k][j][i] - MySim->GravAcc * LocalRhogx3p[k][j][i] * (Localx3[k + 1][j][i] - Localx3[k][j][i]);
            if ( DeltaPot < 0.0 )
              LocalRelPermgx3p[k][j][i] = (1 + Beta) * LocalKrg[k][j][i] - Beta * LocalKrg[k - 1][j][i];
            else
              LocalRelPermgx3p[k][j][i] = (1 + BetaPrime) * LocalKrg[k + 1][j][i] - BetaPrime * LocalKrg[k + 2][j][i];
          }
          else if (( k - zs > 0 && zs != 0 ) || ( zs + zm - k > 1 && zs + zm != mz )) {
            Beta = 0.5 * Localh3[k - 1][j][i] / (Localx3[k - 1][j][i] - Localx3[k - 2][j][i]);
            BetaPrime = 0.5 * Localh3[k][j][i] / (Localx3[k + 1][j][i] - Localx3[k][j][i]);

            DeltaPot = LocalPo[k][j][i] - LocalPo[k - 1][j][i] - MySim->GravAcc * LocalRhoox3m[k][j][i] * (Localx3[k][j][i] - Localx3[k - 1][j][i]);
            if ( DeltaPot < 0.0 )
              LocalRelPermox3m[k][j][i] = (1 + Beta) * LocalKro[k - 1][j][i] - Beta * LocalKro[k - 2][j][i];
            else
              LocalRelPermox3m[k][j][i] = (1 + BetaPrime) * LocalKro[k][j][i] - BetaPrime * LocalKro[k + 1][j][i];

            DeltaPot = LocalPw[k][j][i] - LocalPw[k - 1][j][i] - MySim->GravAcc * LocalRhowx3m[k][j][i] * (Localx3[k][j][i] - Localx3[k - 1][j][i]);
            if ( DeltaPot < 0.0 )
              LocalRelPermwx3m[k][j][i] = (1 + Beta) * LocalKrw[k - 1][j][i] - Beta * LocalKrw[k - 2][j][i];
            else
              LocalRelPermwx3m[k][j][i] = (1 + BetaPrime) * LocalKrw[k][j][i] - BetaPrime * LocalKrw[k + 1][j][i];

            DeltaPot = LocalPg[k][j][i] - LocalPg[k - 1][j][i] - MySim->GravAcc * LocalRhogx3m[k][j][i] * (Localx3[k][j][i] - Localx3[k - 1][j][i]);
            if ( DeltaPot < 0.0 )
              LocalRelPermgx3m[k][j][i] = (1 + Beta) * LocalKrg[k - 1][j][i] - Beta * LocalKrg[k - 2][j][i];
            else
              LocalRelPermgx3m[k][j][i] = (1 + BetaPrime) * LocalKrg[k][j][i] - BetaPrime * LocalKrg[k + 1][j][i];

            Beta = 0.5 * Localh3[k][j][i] / (Localx3[k][j][i] - Localx3[k - 1][j][i]);
            BetaPrime = 0.5 * Localh3[k + 1][j][i] / (Localx3[k + 2][j][i] - Localx3[k + 1][j][i]);

            DeltaPot = LocalPo[k + 1][j][i] - LocalPo[k][j][i] - MySim->GravAcc * LocalRhoox3p[k][j][i] * (Localx3[k + 1][j][i] - Localx3[k][j][i]);
            if ( DeltaPot < 0.0 )
              LocalRelPermox3p[k][j][i] = (1 + Beta) * LocalKro[k][j][i] - Beta * LocalKro[k - 1][j][i];
            else
              LocalRelPermox3p[k][j][i] = (1 + BetaPrime) * LocalKro[k + 1][j][i] - BetaPrime * LocalKro[k + 2][j][i];

            DeltaPot = LocalPw[k + 1][j][i] - LocalPw[k][j][i] - MySim->GravAcc * LocalRhowx3p[k][j][i] * (Localx3[k + 1][j][i] - Localx3[k][j][i]);
            if ( DeltaPot < 0.0 )
              LocalRelPermwx3p[k][j][i] = (1 + Beta) * LocalKrw[k][j][i] - Beta * LocalKrw[k - 1][j][i];
            else
              LocalRelPermwx3p[k][j][i] = (1 + BetaPrime) * LocalKrw[k + 1][j][i] - BetaPrime * LocalKrw[k + 2][j][i];

            DeltaPot = LocalPg[k + 1][j][i] - LocalPg[k][j][i] - MySim->GravAcc * LocalRhogx3p[k][j][i] * (Localx3[k + 1][j][i] - Localx3[k][j][i]);
            if ( DeltaPot < 0.0 )
              LocalRelPermgx3p[k][j][i] = (1 + Beta) * LocalKrg[k][j][i] - Beta * LocalKrg[k - 1][j][i];
            else
              LocalRelPermgx3p[k][j][i] = (1 + BetaPrime) * LocalKrg[k + 1][j][i] - BetaPrime * LocalKrg[k + 2][j][i];
          }
          /* Now we override based on the flow mask */
          if (LocalFlowMask[k][j][i-1] == NO_FLUID_FLOW) {LocalRelPermox1m[i][j][k]=0.0;LocalRelPermwx1m[i][j][k]=0.0;LocalRelPermgx1m[i][j][k]=0.0;}
          if (LocalFlowMask[k][j][i+1] == NO_FLUID_FLOW) {LocalRelPermox1p[i][j][k]=0.0;LocalRelPermwx1p[i][j][k]=0.0;LocalRelPermgx1p[i][j][k]=0.0;}
          if (LocalFlowMask[k][j-1][i] == NO_FLUID_FLOW) {LocalRelPermox2m[i][j][k]=0.0;LocalRelPermwx2m[i][j][k]=0.0;LocalRelPermgx2m[i][j][k]=0.0;}
          if (LocalFlowMask[k][j+1][i] == NO_FLUID_FLOW) {LocalRelPermox2p[i][j][k]=0.0;LocalRelPermwx2p[i][j][k]=0.0;LocalRelPermgx2p[i][j][k]=0.0;}
          if (LocalFlowMask[k-1][j][i] == NO_FLUID_FLOW) {LocalRelPermox3m[i][j][k]=0.0;LocalRelPermwx3m[i][j][k]=0.0;LocalRelPermgx3m[i][j][k]=0.0;}
          if (LocalFlowMask[k+1][j][i] == NO_FLUID_FLOW) {LocalRelPermox3p[i][j][k]=0.0;LocalRelPermwx3p[i][j][k]=0.0;LocalRelPermgx3p[i][j][k]=0.0;}
        }

        if (LocalRelPermox1m[k][j][i] < 0.0) LocalRelPermox1m[k][j][i] = 0.0;
        if (LocalRelPermox2m[k][j][i] < 0.0) LocalRelPermox2m[k][j][i] = 0.0;
        if (LocalRelPermox3m[k][j][i] < 0.0) LocalRelPermox3m[k][j][i] = 0.0;
        if (LocalRelPermox1p[k][j][i] < 0.0) LocalRelPermox1p[k][j][i] = 0.0;
        if (LocalRelPermox2p[k][j][i] < 0.0) LocalRelPermox2p[k][j][i] = 0.0;
        if (LocalRelPermox3p[k][j][i] < 0.0) LocalRelPermox3p[k][j][i] = 0.0;

        if (LocalRelPermwx1m[k][j][i] < 0.0) LocalRelPermwx1m[k][j][i] = 0.0;
        if (LocalRelPermwx2m[k][j][i] < 0.0) LocalRelPermwx2m[k][j][i] = 0.0;
        if (LocalRelPermwx3m[k][j][i] < 0.0) LocalRelPermwx3m[k][j][i] = 0.0;
        if (LocalRelPermwx1p[k][j][i] < 0.0) LocalRelPermwx1p[k][j][i] = 0.0;
        if (LocalRelPermwx2p[k][j][i] < 0.0) LocalRelPermwx2p[k][j][i] = 0.0;
        if (LocalRelPermwx3p[k][j][i] < 0.0) LocalRelPermwx3p[k][j][i] = 0.0;

        if (LocalRelPermgx1m[k][j][i] < 0.0) LocalRelPermgx1m[k][j][i] = 0.0;
        if (LocalRelPermgx2m[k][j][i] < 0.0) LocalRelPermgx2m[k][j][i] = 0.0;
        if (LocalRelPermgx3m[k][j][i] < 0.0) LocalRelPermgx3m[k][j][i] = 0.0;
        if (LocalRelPermgx1p[k][j][i] < 0.0) LocalRelPermgx1p[k][j][i] = 0.0;
        if (LocalRelPermgx2p[k][j][i] < 0.0) LocalRelPermgx2p[k][j][i] = 0.0;
        if (LocalRelPermgx3p[k][j][i] < 0.0) LocalRelPermgx3p[k][j][i] = 0.0;
      }
    }
  }

  /* Restore the new arrays to their rightful place */
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->RelPermox1m, &LocalRelPermox1m);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->RelPermox1p, &LocalRelPermox1p);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->RelPermox2m, &LocalRelPermox2m);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->RelPermox2p, &LocalRelPermox2p);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->RelPermox3m, &LocalRelPermox3m);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->RelPermox3p, &LocalRelPermox3p);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->RelPermwx1m, &LocalRelPermwx1m);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->RelPermwx1p, &LocalRelPermwx1p);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->RelPermwx2m, &LocalRelPermwx2m);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->RelPermwx2p, &LocalRelPermwx2p);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->RelPermwx3m, &LocalRelPermwx3m);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->RelPermwx3p, &LocalRelPermwx3p);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->RelPermgx1m, &LocalRelPermgx1m);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->RelPermgx1p, &LocalRelPermgx1p);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->RelPermgx2m, &LocalRelPermgx2m);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->RelPermgx2p, &LocalRelPermgx2p);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->RelPermgx3m, &LocalRelPermgx3m);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->RelPermgx3p, &LocalRelPermgx3p);CHKERRQ(ierr);

  /* Begin Assembly for vectors */
  ierr = VecAssemblyBegin(MySim->RelPermox1m);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->RelPermox1p);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->RelPermox2m);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->RelPermox2p);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->RelPermox3m);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->RelPermox3p);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->RelPermwx1m);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->RelPermwx1p);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->RelPermwx2m);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->RelPermwx2p);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->RelPermwx3m);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->RelPermwx3p);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->RelPermgx1m);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->RelPermgx1p);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->RelPermgx2m);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->RelPermgx2p);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->RelPermgx3m);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->RelPermgx3p);CHKERRQ(ierr);

  /* And end Assembly */
  ierr = VecAssemblyEnd(MySim->RelPermox1m);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->RelPermox1p);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->RelPermox2m);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->RelPermox2p);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->RelPermox3m);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->RelPermox3p);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->RelPermwx1m);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->RelPermwx1p);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->RelPermwx2m);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->RelPermwx2p);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->RelPermwx3m);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->RelPermwx3p);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->RelPermgx1m);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->RelPermgx1p);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->RelPermgx2m);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->RelPermgx2p);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->RelPermgx3m);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->RelPermgx3p);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DefiantComputeVolumeFactorsAtFaces"
PetscErrorCode DefiantComputeVolumeFactorsAtFaces(BlackOilReservoirSimulation* MySim)
{
  PetscErrorCode ierr;
  PetscInt i, j, k, mx, my, mz, xm, ym, zm, xs, ys, zs;
  PetscReal Beta;
  PetscScalar ***LocalFlowMask, ***LocalPhi;
  PetscScalar ***LocalBo, ***LocalBw, ***LocalBg;
  PetscScalar ***LocalBox1p, ***LocalBox2p, ***LocalBox3p;
  PetscScalar ***LocalBox1m, ***LocalBox2m, ***LocalBox3m;
  PetscScalar ***LocalBwx1p, ***LocalBwx2p, ***LocalBwx3p;
  PetscScalar ***LocalBwx1m, ***LocalBwx2m, ***LocalBwx3m;
  PetscScalar ***LocalBgx1p, ***LocalBgx2p, ***LocalBgx3p;
  PetscScalar ***LocalBgx1m, ***LocalBgx2m, ***LocalBgx3m;

  Vec vecLocalFlowMask, vecLocalPhi;
  Vec vecLocalBo, vecLocalBw, vecLocalBg;

  PetscFunctionBegin;
  /* Get dimensions and extents of the local vectors */
  ierr = DAGetInfo(MySim->SimDA, 0, &mx, &my, &mz, 0, 0, 0, 0, 0, 0, 0);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGetCorners(MySim->SimDA, &xs, &ys, &zs, &xm, &ym, &zm);CHKERRQ(ierr);CHKMEMQ;


  ierr = DAGetLocalVector(MySim->SimDA, &vecLocalFlowMask);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGetLocalVector(MySim->SimDA, &vecLocalPhi);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGetLocalVector(MySim->SimDA, &vecLocalBo);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGetLocalVector(MySim->SimDA, &vecLocalBw);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGetLocalVector(MySim->SimDA, &vecLocalBg);CHKERRQ(ierr);CHKMEMQ;

  ierr = DAGlobalToLocalBegin(MySim->SimDA, MySim-> FlowMask  ,INSERT_VALUES, vecLocalFlowMask  );CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGlobalToLocalEnd(MySim->SimDA, MySim-> FlowMask  ,INSERT_VALUES, vecLocalFlowMask  );CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGlobalToLocalBegin(MySim->SimDA, MySim-> Phi  ,INSERT_VALUES, vecLocalPhi  );CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGlobalToLocalEnd(MySim->SimDA, MySim-> Phi  ,INSERT_VALUES, vecLocalPhi  );CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGlobalToLocalBegin(MySim->SimDA, MySim-> Bo  ,INSERT_VALUES, vecLocalBo  );CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGlobalToLocalEnd(MySim->SimDA, MySim-> Bo  ,INSERT_VALUES, vecLocalBo  );CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGlobalToLocalBegin(MySim->SimDA, MySim-> Bw  ,INSERT_VALUES, vecLocalBw  );CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGlobalToLocalEnd(MySim->SimDA, MySim-> Bw  ,INSERT_VALUES, vecLocalBw  );CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGlobalToLocalBegin(MySim->SimDA, MySim-> Bg  ,INSERT_VALUES, vecLocalBg  );CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGlobalToLocalEnd(MySim->SimDA, MySim-> Bg  ,INSERT_VALUES, vecLocalBg  );CHKERRQ(ierr);CHKMEMQ;

  /* Grab the data for the flow field */
  ierr = DAVecGetArray(MySim->SimDA, vecLocalFlowMask, &LocalFlowMask);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAVecGetArray(MySim->SimDA, vecLocalPhi, &LocalPhi);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAVecGetArray(MySim->SimDA, vecLocalBo, &LocalBo);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAVecGetArray(MySim->SimDA, vecLocalBw, &LocalBw);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAVecGetArray(MySim->SimDA, vecLocalBg, &LocalBg);CHKERRQ(ierr);CHKMEMQ;

  /* Grab the local data for Volume factors at the faces */
  ierr = DAVecGetArray(MySim->SimDA, MySim->Box1p, &LocalBox1p);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAVecGetArray(MySim->SimDA, MySim->Box2p, &LocalBox2p);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAVecGetArray(MySim->SimDA, MySim->Box3p, &LocalBox3p);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAVecGetArray(MySim->SimDA, MySim->Box1m, &LocalBox1m);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAVecGetArray(MySim->SimDA, MySim->Box2m, &LocalBox2m);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAVecGetArray(MySim->SimDA, MySim->Box3m, &LocalBox3m);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAVecGetArray(MySim->SimDA, MySim->Bwx1p, &LocalBwx1p);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAVecGetArray(MySim->SimDA, MySim->Bwx2p, &LocalBwx2p);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAVecGetArray(MySim->SimDA, MySim->Bwx3p, &LocalBwx3p);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAVecGetArray(MySim->SimDA, MySim->Bwx1m, &LocalBwx1m);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAVecGetArray(MySim->SimDA, MySim->Bwx2m, &LocalBwx2m);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAVecGetArray(MySim->SimDA, MySim->Bwx3m, &LocalBwx3m);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAVecGetArray(MySim->SimDA, MySim->Bgx1p, &LocalBgx1p);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAVecGetArray(MySim->SimDA, MySim->Bgx2p, &LocalBgx2p);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAVecGetArray(MySim->SimDA, MySim->Bgx3p, &LocalBgx3p);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAVecGetArray(MySim->SimDA, MySim->Bgx1m, &LocalBgx1m);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAVecGetArray(MySim->SimDA, MySim->Bgx2m, &LocalBgx2m);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAVecGetArray(MySim->SimDA, MySim->Bgx3m, &LocalBgx3m);CHKERRQ(ierr);CHKMEMQ;

  /* Set the gas oil solubility */
  for (k = zs; k < zs + zm; k++) {
    for (j = ys; j < ys + ym; j++) {
      for (i = xs; i < xs + xm; i++) {
        if (i == 0 || j == 0 || k == 0 || i == mx - 1 || j == my - 1 || k == mz
            - 1) {
        } else if (ABS(LocalFlowMask[k][j][i]-FLUID_FLOW) < EPSILON) {
          Beta = LocalPhi[k][j][i] / LocalPhi[k][j][i - 1];
          /* Density at Face */
          LocalBox1m[k][j][i] = Beta * LocalBo[k][j][i] + (1.0 - Beta)
              * LocalBo[k][j][i - 1];
          LocalBwx1m[k][j][i] = Beta * LocalBw[k][j][i] + (1.0 - Beta)
              * LocalBw[k][j][i - 1];
          LocalBgx1m[k][j][i] = Beta * LocalBg[k][j][i] + (1.0 - Beta)
              * LocalBg[k][j][i - 1];

          Beta = LocalPhi[k][j][i] / LocalPhi[k][j][i + 1];
          /* Density at Face */
          LocalBox1p[k][j][i] = Beta * LocalBo[k][j][i] + (1.0 - Beta)
              * LocalBo[k][j][i + 1];
          LocalBwx1p[k][j][i] = Beta * LocalBw[k][j][i] + (1.0 - Beta)
              * LocalBw[k][j][i + 1];
          LocalBgx1p[k][j][i] = Beta * LocalBg[k][j][i] + (1.0 - Beta)
              * LocalBg[k][j][i + 1];

          Beta = LocalPhi[k][j][i] / LocalPhi[k][j - 1][i];
          /* Density at Face */
          LocalBox2m[k][j][i] = Beta * LocalBo[k][j][i] + (1.0 - Beta)
              * LocalBo[k][j - 1][i];
          LocalBwx2m[k][j][i] = Beta * LocalBw[k][j][i] + (1.0 - Beta)
              * LocalBw[k][j - 1][i];
          LocalBgx2m[k][j][i] = Beta * LocalBg[k][j][i] + (1.0 - Beta)
              * LocalBg[k][j - 1][i];

          Beta = LocalPhi[k][j][i] / LocalPhi[k][j + 1][i];
          /* Density at Face */
          LocalBox2p[k][j][i] = Beta * LocalBo[k][j][i] + (1.0 - Beta)
              * LocalBo[k][j + 1][i];
          LocalBwx2p[k][j][i] = Beta * LocalBw[k][j][i] + (1.0 - Beta)
              * LocalBw[k][j + 1][i];
          LocalBgx2p[k][j][i] = Beta * LocalBg[k][j][i] + (1.0 - Beta)
              * LocalBg[k][j + 1][i];

          Beta = LocalPhi[k][j][i] / LocalPhi[k - 1][j][i];
          /* Density at Face */
          LocalBox3m[k][j][i] = Beta * LocalBo[k][j][i] + (1.0 - Beta)
              * LocalBo[k - 1][j][i];
          LocalBwx3m[k][j][i] = Beta * LocalBw[k][j][i] + (1.0 - Beta)
              * LocalBw[k - 1][j][i];
          LocalBgx3m[k][j][i] = Beta * LocalBg[k][j][i] + (1.0 - Beta)
              * LocalBg[k - 1][j][i];

          Beta = LocalPhi[k][j][i] / LocalPhi[k + 1][j][i];
          /* Density at Face */
          LocalBox3p[k][j][i] = Beta * LocalBo[k][j][i] + (1.0 - Beta)
              * LocalBo[k + 1][j][i];
          LocalBwx3p[k][j][i] = Beta * LocalBw[k][j][i] + (1.0 - Beta)
              * LocalBw[k + 1][j][i];
          LocalBgx3p[k][j][i] = Beta * LocalBg[k][j][i] + (1.0 - Beta)
              * LocalBg[k + 1][j][i];

        }
      }
    }
  }

  /* Restore data for Volume factors at the faces */
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Box1p, &LocalBox1p);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Box2p, &LocalBox2p);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Box3p, &LocalBox3p);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Box1m, &LocalBox1m);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Box2m, &LocalBox2m);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Box3m, &LocalBox3m);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Bwx1p, &LocalBwx1p);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Bwx2p, &LocalBwx2p);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Bwx3p, &LocalBwx3p);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Bwx1m, &LocalBwx1m);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Bwx2m, &LocalBwx2m);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Bwx3m, &LocalBwx3m);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Bgx1p, &LocalBgx1p);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Bgx2p, &LocalBgx2p);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Bgx3p, &LocalBgx3p);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Bgx1m, &LocalBgx1m);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Bgx2m, &LocalBgx2m);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Bgx3m, &LocalBgx3m);CHKERRQ(ierr);CHKMEMQ;

  /* Assemble the vectors */
  ierr = VecAssemblyBegin(MySim->Box1p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecAssemblyBegin(MySim->Box2p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecAssemblyBegin(MySim->Box3p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecAssemblyBegin(MySim->Box1m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecAssemblyBegin(MySim->Box2m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecAssemblyBegin(MySim->Box3m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecAssemblyBegin(MySim->Bwx1p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecAssemblyBegin(MySim->Bwx2p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecAssemblyBegin(MySim->Bwx3p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecAssemblyBegin(MySim->Bwx1m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecAssemblyBegin(MySim->Bwx2m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecAssemblyBegin(MySim->Bwx3m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecAssemblyBegin(MySim->Bgx1p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecAssemblyBegin(MySim->Bgx2p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecAssemblyBegin(MySim->Bgx3p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecAssemblyBegin(MySim->Bgx1m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecAssemblyBegin(MySim->Bgx2m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecAssemblyBegin(MySim->Bgx3m);CHKERRQ(ierr);CHKMEMQ;

  /* finish the assembly */
  ierr = VecAssemblyEnd(MySim->Box1p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecAssemblyEnd(MySim->Box2p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecAssemblyEnd(MySim->Box3p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecAssemblyEnd(MySim->Box1m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecAssemblyEnd(MySim->Box2m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecAssemblyEnd(MySim->Box3m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecAssemblyEnd(MySim->Bwx1p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecAssemblyEnd(MySim->Bwx2p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecAssemblyEnd(MySim->Bwx3p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecAssemblyEnd(MySim->Bwx1m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecAssemblyEnd(MySim->Bwx2m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecAssemblyEnd(MySim->Bwx3m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecAssemblyEnd(MySim->Bgx1p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecAssemblyEnd(MySim->Bgx2p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecAssemblyEnd(MySim->Bgx3p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecAssemblyEnd(MySim->Bgx1m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecAssemblyEnd(MySim->Bgx2m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecAssemblyEnd(MySim->Bgx3m);CHKERRQ(ierr);CHKMEMQ;

  PetscFunctionReturn(0);
}
