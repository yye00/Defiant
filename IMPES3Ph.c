/*
 * IMPES3Ph.c
 *
 *  Created on: Sep 6, 2009
 *      Author: yye00
 */

#include "Defiant.h"

#undef __FUNCT__
#define __FUNCT__ "DefiantIMPES3PhHandleWells"
extern PetscErrorCode DefiantIMPES3PhHandleWells(BlackOilReservoirSimulation* MySim) {
  PetscErrorCode ierr;
  PetscInt mx, my, mz, xm, ym, zm, xs, ys, zs;
  PetscScalar re;
  PetscScalar v[1];
  MatStencil row, col;
  /* Well handling variables */
  PetscInt PerfIDMine, PerfIDOther, WellID;
  PetscScalar apo, apw, apg;
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
  PetscScalar ***LocalPcow, ***LocalPcog;
  PetscScalar ***LocalRHS;
  PetscScalar ***LocalBo, ***LocalBw, ***LocalBg;

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
  /* Grab the local capillary pressure */
  ierr = DAVecGetArray(MySim->SimDA, MySim->Pcow, &LocalPcow);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Pcog, &LocalPcog);CHKERRQ(ierr);
  /* Grab the local data for volume factors at the cell centers */
  ierr = DAVecGetArray(MySim->SimDA, MySim->Bo, &LocalBo);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Bw, &LocalBw);CHKERRQ(ierr);
  ierr = DAVecGetArray(MySim->SimDA, MySim->Bg, &LocalBg);CHKERRQ(ierr);
  /* Grab the local RHS */
  ierr = DAVecGetArray(MySim->SimDA, MySim->RHS, &LocalRHS);CHKERRQ(ierr);

  for (WellID = 0; WellID < MySim->NumberOfWells; WellID++) {
    for (PerfIDMine = 0; PerfIDMine
        < (MySim->Wells[WellID]).NumberOfPerforations; PerfIDMine++) {
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

  for (WellID = 0; WellID < MySim->NumberOfWells; WellID++) {
    for (PerfIDMine = 0; PerfIDMine
        < (MySim->Wells[WellID]).NumberOfPerforations; PerfIDMine++) {
      /* grab my index */
      MyI = (MySim->Wells[WellID]).Perforations[PerfIDMine].I;
      MyJ = (MySim->Wells[WellID]).Perforations[PerfIDMine].J;
      MyK = (MySim->Wells[WellID]).Perforations[PerfIDMine].K;
      /* set the row for the insert */
      row.i = MyI;
      row.j = MyJ;
      row.k = MyK;
      /* move terms to my RHS */
      qo = 0.0;
      qw = 0.0;
      /* move terms to my Ap */
      apo = 0.0; apw = 0.0;
      if ((MySim->Wells[WellID]).NumberOfPerforations == 1) {
        MyI = (MySim->Wells[WellID]).Perforations[PerfIDMine].I;
        MyJ = (MySim->Wells[WellID]).Perforations[PerfIDMine].J;
        MyK = (MySim->Wells[WellID]).Perforations[PerfIDMine].K;
        /* for all other perforations belonging to the same well do this */
        MyI = (MySim->Wells[WellID]).Perforations[PerfIDMine].I;
        MyJ = (MySim->Wells[WellID]).Perforations[PerfIDMine].J;
        MyK = (MySim->Wells[WellID]).Perforations[PerfIDMine].K;
        qo  = (MySim->Wells[WellID]).Perforations[PerfIDMine].WellIndex
            * LocalKro[MyK][MyJ][MyI] / LocalMuo[MyK][MyJ][MyI]
            * ((MySim->Wells[WellID]).Perforations[PerfIDMine].BHPo
            - MySim->GravAcc * LocalRhoo[MyK][MyJ][MyI]
            * ((MySim->Wells[WellID]).Perforations[PerfIDMine].zbh
            - Localx3[MyK][MyJ][MyI])) * Localh1[MyK][MyJ][MyI]
            * Localh2[MyK][MyJ][MyI] * Localh3[MyK][MyJ][MyI]
            / LocalBo[MyK][MyJ][MyI];
        qw  = (MySim->Wells[WellID]).Perforations[PerfIDMine].WellIndex
            * LocalKrw[MyK][MyJ][MyI] / LocalMuw[MyK][MyJ][MyI]
            * ((MySim->Wells[WellID]).Perforations[PerfIDMine].BHPw
            + LocalPcow[MyK][MyJ][MyI] - MySim->GravAcc
            * LocalRhow[MyK][MyJ][MyI]
            * ((MySim->Wells[WellID]).Perforations[PerfIDMine].zbh
            - Localx3[MyK][MyJ][MyI])) * Localh1[MyK][MyJ][MyI]
            * Localh2[MyK][MyJ][MyI] * Localh3[MyK][MyJ][MyI]
            / LocalBw[MyK][MyJ][MyI];
        qg  = (MySim->Wells[WellID]).Perforations[PerfIDMine].WellIndex
            * LocalKrg[MyK][MyJ][MyI] / LocalMug[MyK][MyJ][MyI]
            * ((MySim->Wells[WellID]).Perforations[PerfIDMine].BHPg
            - LocalPcog[MyK][MyJ][MyI] - MySim->GravAcc
            * LocalRhow[MyK][MyJ][MyI]
            * ((MySim->Wells[WellID]).Perforations[PerfIDMine].zbh
            - Localx3[MyK][MyJ][MyI])) * Localh1[MyK][MyJ][MyI]
            * Localh2[MyK][MyJ][MyI] * Localh3[MyK][MyJ][MyI]
            / LocalBg[MyK][MyJ][MyI];
        apo = (MySim->Wells[WellID]).Perforations[PerfIDMine].WellIndex
            * LocalKro[MyK][MyJ][MyI] / LocalMuo[MyK][MyJ][MyI]
            * Localh1[MyK][MyJ][MyI] * Localh2[MyK][MyJ][MyI] * Localh3[MyK][MyJ][MyI]
            / LocalBo[MyK][MyJ][MyI];
        apw = (MySim->Wells[WellID]).Perforations[PerfIDMine].WellIndex
            * LocalKrw[MyK][MyJ][MyI] / LocalMuw[MyK][MyJ][MyI]
            * Localh1[MyK][MyJ][MyI] * Localh2[MyK][MyJ][MyI] * Localh3[MyK][MyJ][MyI]
            / LocalBw[MyK][MyJ][MyI];
        apg = (MySim->Wells[WellID]).Perforations[PerfIDMine].WellIndex
            * LocalKrg[MyK][MyJ][MyI] / LocalMug[MyK][MyJ][MyI]
            * Localh1[MyK][MyJ][MyI] * Localh2[MyK][MyJ][MyI] * Localh3[MyK][MyJ][MyI]
            / LocalBg[MyK][MyJ][MyI];
      } else if ((MySim->Wells[WellID]).NumberOfPerforations > 1) {
        for (PerfIDOther = 0; PerfIDOther
            < (MySim->Wells[WellID]).NumberOfPerforations; PerfIDOther++) {
          OtherI = (MySim->Wells[WellID]).Perforations[PerfIDOther].I;
          OtherJ = (MySim->Wells[WellID]).Perforations[PerfIDOther].J;
          OtherK = (MySim->Wells[WellID]).Perforations[PerfIDOther].K;
          if (PerfIDOther != PerfIDMine) {
            /* for all other perforations belonging to the same well do this */
            OtherI = (MySim->Wells[WellID]).Perforations[PerfIDOther].I;
            OtherJ = (MySim->Wells[WellID]).Perforations[PerfIDOther].J;
            OtherK = (MySim->Wells[WellID]).Perforations[PerfIDOther].K;
            qo += (MySim->Wells[WellID]).Perforations[PerfIDOther].WellIndex
                * LocalKro[OtherK][OtherJ][OtherI]
                / LocalMuo[OtherK][OtherJ][OtherI]
                * ((MySim->Wells[WellID]).Perforations[PerfIDOther].BHPo
                - MySim->GravAcc * LocalRhoo[OtherK][OtherJ][OtherI]
                * ((MySim->Wells[WellID]).Perforations[PerfIDOther].zbh
                - Localx3[OtherK][OtherJ][OtherI]));
            qw += (MySim->Wells[WellID]).Perforations[PerfIDOther].WellIndex
                * LocalKrw[OtherK][OtherJ][OtherI]
                / LocalMuw[OtherK][OtherJ][OtherI]
                * ((MySim->Wells[WellID]).Perforations[PerfIDOther].BHPw
                + LocalPcow[OtherK][OtherJ][OtherI] - MySim->GravAcc
                * LocalRhow[OtherK][OtherJ][OtherI]
                * ((MySim->Wells[WellID]).Perforations[PerfIDOther].zbh
                - Localx3[OtherK][OtherJ][OtherI]));
            LocalRHS[MyK][MyJ][MyI] =   LocalRHS[MyK][MyJ][MyI] - qo - qw;
            apo = -(MySim->Wells[WellID]).Perforations[PerfIDOther].WellIndex
                * LocalKro[OtherK][OtherJ][OtherI] / LocalMuo[OtherK][OtherJ][OtherI]
                * Localh1[OtherK][OtherJ][OtherI] * Localh2[OtherK][OtherJ][OtherI] * Localh3[OtherK][OtherJ][OtherI];
            apw = - (MySim->Wells[WellID]).Perforations[PerfIDOther].WellIndex
                * LocalKrw[OtherK][OtherJ][OtherI] / LocalMuw[OtherK][OtherJ][OtherI]
                * Localh1[OtherK][OtherJ][OtherI] * Localh2[OtherK][OtherJ][OtherI] * Localh3[OtherK][OtherJ][OtherI];
          }
        }

        /* Store calculated values in perforation structure */
        MySim->Wells[WellID].Perforations[PerfIDMine].Apo = apo;
        MySim->Wells[WellID].Perforations[PerfIDMine].Apw = apw;
        MySim->Wells[WellID].Perforations[PerfIDMine].Apg = apg;

        LocalRHS[MyK][MyJ][MyI] +=  MySim->Wells[WellID].TotalQo/LocalRhoo[MyK][MyJ][MyI]
                                  + MySim->Wells[WellID].TotalQw/LocalRhow[MyK][MyJ][MyI];
      }
    }
  }

  /* Restore the RHS to it's place */
  ierr = DAVecRestoreArray(MySim->SimDA, MySim->RHS,&LocalRHS);CHKERRQ(ierr);
  /* Finalize by assembling the matrix */
  ierr = MatAssemblyBegin(MySim->A, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(MySim->A, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


