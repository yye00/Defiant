/*
 * Production.c
 *
 *  Created on: Oct 8, 2009
 *      Author: yye00
 */

#include "Defiant.h"

#undef __FUNCT__
#define __FUNCT__ "DefiantBlackOilSyncPerfOwners"
extern PetscErrorCode DefiantBlackOilSyncPerfOwners(BlackOilReservoirSimulation* MySim)
{
  PetscErrorCode ierr;
  PetscInt i, rank, TotalNumberOfPerfs;
  PetscInt mx, my, mz, xm, ym, zm, xs, ys, zs;
  /* Well handling variables */
  PetscInt PerfIDMine, WellID;
  PetscInt MyI, MyJ, MyK;
  /* temporary array */
  PetscInt *InBuffer, *OutBuffer;

  PetscFunctionBegin;
  /* Get dimensions and extents of the local vectors */
  ierr = DAGetInfo(MySim->SimDA, 0, &mx, &my, &mz, 0, 0, 0, 0, 0, 0, 0);CHKERRQ(ierr);
  ierr = DAGetCorners(MySim->SimDA, &xs, &ys, &zs, &xm, &ym, &zm);CHKERRQ(ierr);
  /* Get the current rank */
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);CHKMEMQ;

  /* Allocate memory to the in and out Buffers */
  /* Find the total number of perforations */
  TotalNumberOfPerfs = 0;
  for (WellID = 0; WellID < MySim->NumberOfWells; WellID++)
      for (PerfIDMine = 0; PerfIDMine < (MySim->Wells[WellID]).NumberOfPerforations; PerfIDMine++)
        TotalNumberOfPerfs++;
  /* Allocate the In and out Buffers */
  ierr = PetscMalloc(sizeof(PetscInt)*TotalNumberOfPerfs, &InBuffer);CHKERRQ(ierr);CHKMEMQ;
  ierr = PetscMalloc(sizeof(PetscInt)*TotalNumberOfPerfs, &OutBuffer);CHKERRQ(ierr);CHKMEMQ;

  /* Pack all the perf owner location int the in buffer */
  i = 0;
  for (WellID = 0; WellID < MySim->NumberOfWells; WellID++) {
    for (PerfIDMine = 0; PerfIDMine
        < (MySim->Wells[WellID]).NumberOfPerforations; PerfIDMine++) {

      MyI = (MySim->Wells[WellID]).Perforations[PerfIDMine].I;CHKMEMQ;
      MyJ = (MySim->Wells[WellID]).Perforations[PerfIDMine].J;CHKMEMQ;
      MyK = (MySim->Wells[WellID]).Perforations[PerfIDMine].K;CHKMEMQ;

      if (MyI >= xs && MyI < xs+xm && MyJ >= ys && MyJ < ys+ym && MyK >= zs && MyK < zs+zm)
        MySim->Wells[WellID].Perforations[PerfIDMine].OwnerRank = rank;
      else
        MySim->Wells[WellID].Perforations[PerfIDMine].OwnerRank = -1;

      InBuffer[i] = MySim->Wells[WellID].Perforations[PerfIDMine].OwnerRank;
      i++;
    }
  }

  /* Make the syncrhonization call */
  ierr = MPI_Allreduce(InBuffer, OutBuffer, TotalNumberOfPerfs, MPI_INT, MPI_MAX, PETSC_COMM_WORLD);CHKERRQ(ierr);CHKMEMQ;

  /* Now unpack all the perf owner location int the in buffer */
  i = 0;
  for (WellID = 0; WellID < MySim->NumberOfWells; WellID++) {
    for (PerfIDMine = 0; PerfIDMine
        < (MySim->Wells[WellID]).NumberOfPerforations; PerfIDMine++) {

       MySim->Wells[WellID].Perforations[PerfIDMine].OwnerRank = OutBuffer[i];
       i++;
    }
  }

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DefiantBlackOilComputePerfIndicesSyncPerfs"
extern PetscErrorCode DefiantBlackOilComputePerfIndicesSyncPerfs(BlackOilReservoirSimulation* MySim)
{
  PetscErrorCode ierr;
  PetscInt mx, my, mz, xm, ym, zm, xs, ys, zs;
  PetscInt i, rank;
  /* Well handling variables */
  PetscInt PerfIDMine, WellID;
  PetscInt MyI, MyJ, MyK;
  PetscScalar re;
  /* Local values for variables*/
  PetscScalar ***Localx1,   ***Localx2,   ***Localx3;
  PetscScalar ***Localh1,   ***Localh2,   ***Localh3;
  PetscScalar ***LocalSo,   ***LocalSw,   ***LocalSg;
  PetscScalar ***LocalPo,   ***LocalPw,   ***LocalPg;
  PetscScalar ***LocalMuo,  ***LocalMuw,  ***LocalMug;
  PetscScalar ***LocalRhoo, ***LocalRhow, ***LocalRhog;
  PetscScalar ***LocalK11,  ***LocalK22,  ***LocalK33;
  PetscScalar ***LocalKro,  ***LocalKrw,  ***LocalKrg;
  PetscScalar ***LocalBo,   ***LocalBw,    ***LocalBg;
  PetscScalar ***LocalPcow, ***LocalPcog;
  /* temporary array */
  PetscScalar *InBuffer;

  PetscFunctionBegin;
  /* Allocate the In and out Buffers, we have 27 doubles at each perforation */
  ierr = PetscMalloc(sizeof(PetscScalar)*27, &InBuffer);CHKERRQ(ierr);CHKMEMQ;

  /* Get dimensions and extents of the local vectors */
  ierr = DAGetInfo(MySim->SimDA, 0, &mx, &my, &mz, 0, 0, 0, 0, 0, 0, 0);CHKERRQ(ierr);
  ierr = DAGetCorners(MySim->SimDA, &xs, &ys, &zs, &xm, &ym, &zm);CHKERRQ(ierr);
  /* Grab the local coordiantes */
  ierr = DAVecGetArray(MySim->SimDA, MySim->x1, &Localx1);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAVecGetArray(MySim->SimDA, MySim->x2, &Localx2);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAVecGetArray(MySim->SimDA, MySim->x3, &Localx3);CHKERRQ(ierr);CHKMEMQ;
  /* Grab the local geometry */
  ierr = DAVecGetArray(MySim->SimDA, MySim->h1, &Localh1);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAVecGetArray(MySim->SimDA, MySim->h2, &Localh2);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAVecGetArray(MySim->SimDA, MySim->h3, &Localh3);CHKERRQ(ierr);CHKMEMQ;
  /* Grab the Saturations */
  ierr = DAVecGetArray(MySim->SimDA, MySim->So, &LocalSo);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAVecGetArray(MySim->SimDA, MySim->Sw, &LocalSw);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAVecGetArray(MySim->SimDA, MySim->Sg, &LocalSg);CHKERRQ(ierr);CHKMEMQ;
  /* Grab the Pressures */
  ierr = DAVecGetArray(MySim->SimDA, MySim->Po, &LocalPo);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAVecGetArray(MySim->SimDA, MySim->Pw, &LocalPw);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAVecGetArray(MySim->SimDA, MySim->Pg, &LocalPg);CHKERRQ(ierr);CHKMEMQ;
  /* Grab the Densities */
  ierr = DAVecGetArray(MySim->SimDA, MySim->Rhoo, &LocalRhoo);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAVecGetArray(MySim->SimDA, MySim->Rhow, &LocalRhow);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAVecGetArray(MySim->SimDA, MySim->Rhog, &LocalRhog);CHKERRQ(ierr);CHKMEMQ;
  /* Grab the local viscosities*/
  ierr = DAVecGetArray(MySim->SimDA, MySim->Muo, &LocalMuo);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAVecGetArray(MySim->SimDA, MySim->Muw, &LocalMuw);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAVecGetArray(MySim->SimDA, MySim->Mug, &LocalMug);CHKERRQ(ierr);CHKMEMQ;
  /* Grab the local permeabilities */
  ierr = DAVecGetArray(MySim->SimDA, MySim->K11, &LocalK11);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAVecGetArray(MySim->SimDA, MySim->K22, &LocalK22);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAVecGetArray(MySim->SimDA, MySim->K33, &LocalK33);CHKERRQ(ierr);CHKMEMQ;
  /* Grab the local permeabilities */
  ierr = DAVecGetArray(MySim->SimDA, MySim->Kro, &LocalKro);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAVecGetArray(MySim->SimDA, MySim->Krw, &LocalKrw);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAVecGetArray(MySim->SimDA, MySim->Krg, &LocalKrg);CHKERRQ(ierr);CHKMEMQ;
  /* Grab the local capillary pressure */
  ierr = DAVecGetArray(MySim->SimDA, MySim->Pcow, &LocalPcow);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAVecGetArray(MySim->SimDA, MySim->Pcog, &LocalPcog);CHKERRQ(ierr);CHKMEMQ;
  /* Grab the local data for volume factors at the cell centers */
  ierr = DAVecGetArray(MySim->SimDA, MySim->Bo, &LocalBo);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAVecGetArray(MySim->SimDA, MySim->Bw, &LocalBw);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAVecGetArray(MySim->SimDA, MySim->Bg, &LocalBg);CHKERRQ(ierr);CHKMEMQ;

  /* Get the current rank */
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);CHKMEMQ;
  /* Now we sync all the flow rates and BHP's across all perforations */
  /* At the end of this loop, all information across all wells is uniform and up-to-date */
  for (WellID = 0; WellID < MySim->NumberOfWells; WellID++) {
    for (PerfIDMine = 0; PerfIDMine
        < (MySim->Wells[WellID]).NumberOfPerforations; PerfIDMine++) {

      MyI = (MySim->Wells[WellID]).Perforations[PerfIDMine].I;CHKMEMQ;
      MyJ = (MySim->Wells[WellID]).Perforations[PerfIDMine].J;CHKMEMQ;
      MyK = (MySim->Wells[WellID]).Perforations[PerfIDMine].K;CHKMEMQ;

      if (MyI >= xs && MyI < xs+xm && MyJ >= ys && MyJ < ys+ym && MyK >= zs && MyK < zs+zm)
      {
        MySim->Wells[WellID].Perforations[PerfIDMine].So = LocalSo[MyK][MyJ][MyI];CHKMEMQ;
        MySim->Wells[WellID].Perforations[PerfIDMine].Sw = LocalSw[MyK][MyJ][MyI];CHKMEMQ;
        MySim->Wells[WellID].Perforations[PerfIDMine].Sg = LocalSg[MyK][MyJ][MyI];CHKMEMQ;

        MySim->Wells[WellID].Perforations[PerfIDMine].Po = LocalPo[MyK][MyJ][MyI];CHKMEMQ;
        MySim->Wells[WellID].Perforations[PerfIDMine].Pw = LocalPw[MyK][MyJ][MyI];CHKMEMQ;
        MySim->Wells[WellID].Perforations[PerfIDMine].Pg = LocalPg[MyK][MyJ][MyI];CHKMEMQ;

        MySim->Wells[WellID].Perforations[PerfIDMine].Kro = LocalKro[MyK][MyJ][MyI];CHKMEMQ;
        MySim->Wells[WellID].Perforations[PerfIDMine].Krw = LocalKrw[MyK][MyJ][MyI];CHKMEMQ;
        MySim->Wells[WellID].Perforations[PerfIDMine].Krg = LocalKrg[MyK][MyJ][MyI];CHKMEMQ;

        MySim->Wells[WellID].Perforations[PerfIDMine].Muo = LocalMuo[MyK][MyJ][MyI];CHKMEMQ;
        MySim->Wells[WellID].Perforations[PerfIDMine].Muw = LocalMuw[MyK][MyJ][MyI];CHKMEMQ;
        MySim->Wells[WellID].Perforations[PerfIDMine].Mug = LocalMug[MyK][MyJ][MyI];CHKMEMQ;

        MySim->Wells[WellID].Perforations[PerfIDMine].Rhoo = LocalRhoo[MyK][MyJ][MyI];CHKMEMQ;
        MySim->Wells[WellID].Perforations[PerfIDMine].Rhow = LocalRhow[MyK][MyJ][MyI];CHKMEMQ;
        MySim->Wells[WellID].Perforations[PerfIDMine].Rhog = LocalRhog[MyK][MyJ][MyI];CHKMEMQ;

        MySim->Wells[WellID].Perforations[PerfIDMine].Pcow = LocalPcow[MyK][MyJ][MyI];CHKMEMQ;
        MySim->Wells[WellID].Perforations[PerfIDMine].Pcog = LocalPcog[MyK][MyJ][MyI];CHKMEMQ;

        MySim->Wells[WellID].Perforations[PerfIDMine].Bo = LocalBo[MyK][MyJ][MyI];CHKMEMQ;
        MySim->Wells[WellID].Perforations[PerfIDMine].Bw = LocalBw[MyK][MyJ][MyI];CHKMEMQ;
        MySim->Wells[WellID].Perforations[PerfIDMine].Bg = LocalBg[MyK][MyJ][MyI];CHKMEMQ;

        MySim->Wells[WellID].Perforations[PerfIDMine].x1 = Localx1[MyK][MyJ][MyI];CHKMEMQ;
        MySim->Wells[WellID].Perforations[PerfIDMine].x2 = Localx2[MyK][MyJ][MyI];CHKMEMQ;
        MySim->Wells[WellID].Perforations[PerfIDMine].x3 = Localx3[MyK][MyJ][MyI];CHKMEMQ;

        MySim->Wells[WellID].Perforations[PerfIDMine].h1 = Localh1[MyK][MyJ][MyI];CHKMEMQ;
        MySim->Wells[WellID].Perforations[PerfIDMine].h2 = Localh2[MyK][MyJ][MyI];CHKMEMQ;
        MySim->Wells[WellID].Perforations[PerfIDMine].h3 = Localh3[MyK][MyJ][MyI];CHKMEMQ;

        /* we also compute the well indices */
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

      i = 0;
      InBuffer[i] = MySim->Wells[WellID].Perforations[PerfIDMine].Po;i++;
      InBuffer[i] = MySim->Wells[WellID].Perforations[PerfIDMine].Pw;i++;
      InBuffer[i] = MySim->Wells[WellID].Perforations[PerfIDMine].Pg;i++;

      InBuffer[i] = MySim->Wells[WellID].Perforations[PerfIDMine].So;i++;
      InBuffer[i] = MySim->Wells[WellID].Perforations[PerfIDMine].Sw;i++;
      InBuffer[i] = MySim->Wells[WellID].Perforations[PerfIDMine].Sg;i++;

      InBuffer[i] = MySim->Wells[WellID].Perforations[PerfIDMine].Rhoo;i++;
      InBuffer[i] = MySim->Wells[WellID].Perforations[PerfIDMine].Rhow;i++;
      InBuffer[i] = MySim->Wells[WellID].Perforations[PerfIDMine].Rhog;i++;

      InBuffer[i] = MySim->Wells[WellID].Perforations[PerfIDMine].Bo;i++;
      InBuffer[i] = MySim->Wells[WellID].Perforations[PerfIDMine].Bw;i++;
      InBuffer[i] = MySim->Wells[WellID].Perforations[PerfIDMine].Bg;i++;

      InBuffer[i] = MySim->Wells[WellID].Perforations[PerfIDMine].Kro;i++;
      InBuffer[i] = MySim->Wells[WellID].Perforations[PerfIDMine].Krw;i++;
      InBuffer[i] = MySim->Wells[WellID].Perforations[PerfIDMine].Krg;i++;

      InBuffer[i] = MySim->Wells[WellID].Perforations[PerfIDMine].Muo;i++;
      InBuffer[i] = MySim->Wells[WellID].Perforations[PerfIDMine].Muw;i++;
      InBuffer[i] = MySim->Wells[WellID].Perforations[PerfIDMine].Mug;i++;

      InBuffer[i] = MySim->Wells[WellID].Perforations[PerfIDMine].x1;i++;
      InBuffer[i] = MySim->Wells[WellID].Perforations[PerfIDMine].x2;i++;
      InBuffer[i] = MySim->Wells[WellID].Perforations[PerfIDMine].x3;i++;

      InBuffer[i] = MySim->Wells[WellID].Perforations[PerfIDMine].h1;i++;
      InBuffer[i] = MySim->Wells[WellID].Perforations[PerfIDMine].h2;i++;
      InBuffer[i] = MySim->Wells[WellID].Perforations[PerfIDMine].h3;i++;

      InBuffer[i] = MySim->Wells[WellID].Perforations[PerfIDMine].Pcow;i++;
      InBuffer[i] = MySim->Wells[WellID].Perforations[PerfIDMine].Pcog;i++;
      InBuffer[i] = MySim->Wells[WellID].Perforations[PerfIDMine].WellIndex;

      /* Sync properties at the perforation locations */
      ierr = MPI_Bcast(InBuffer, 27, MPI_DOUBLE, MySim->Wells[WellID].Perforations[PerfIDMine].OwnerRank, PETSC_COMM_WORLD);CHKERRQ(ierr);CHKMEMQ;

      /* and now unpack the data */
      i = 0;

      MySim->Wells[WellID].Perforations[PerfIDMine].Po = InBuffer[i];i++;
      MySim->Wells[WellID].Perforations[PerfIDMine].Pw = InBuffer[i];i++;
      MySim->Wells[WellID].Perforations[PerfIDMine].Pg = InBuffer[i];i++;

      MySim->Wells[WellID].Perforations[PerfIDMine].So = InBuffer[i];i++;
      MySim->Wells[WellID].Perforations[PerfIDMine].Sw = InBuffer[i];i++;
      MySim->Wells[WellID].Perforations[PerfIDMine].Sg = InBuffer[i];i++;

      MySim->Wells[WellID].Perforations[PerfIDMine].Rhoo = InBuffer[i];i++;
      MySim->Wells[WellID].Perforations[PerfIDMine].Rhow = InBuffer[i];i++;
      MySim->Wells[WellID].Perforations[PerfIDMine].Rhog = InBuffer[i];i++;

      MySim->Wells[WellID].Perforations[PerfIDMine].Bo = InBuffer[i];i++;
      MySim->Wells[WellID].Perforations[PerfIDMine].Bw = InBuffer[i];i++;
      MySim->Wells[WellID].Perforations[PerfIDMine].Bg = InBuffer[i];i++;

      MySim->Wells[WellID].Perforations[PerfIDMine].Kro = InBuffer[i];i++;
      MySim->Wells[WellID].Perforations[PerfIDMine].Krw = InBuffer[i];i++;
      MySim->Wells[WellID].Perforations[PerfIDMine].Krg = InBuffer[i];i++;

      MySim->Wells[WellID].Perforations[PerfIDMine].Muo = InBuffer[i];i++;
      MySim->Wells[WellID].Perforations[PerfIDMine].Muw = InBuffer[i];i++;
      MySim->Wells[WellID].Perforations[PerfIDMine].Mug = InBuffer[i];i++;

      MySim->Wells[WellID].Perforations[PerfIDMine].x1 = InBuffer[i];i++;
      MySim->Wells[WellID].Perforations[PerfIDMine].x2 = InBuffer[i];i++;
      MySim->Wells[WellID].Perforations[PerfIDMine].x3 = InBuffer[i];i++;

      MySim->Wells[WellID].Perforations[PerfIDMine].h1 = InBuffer[i];i++;
      MySim->Wells[WellID].Perforations[PerfIDMine].h2 = InBuffer[i];i++;
      MySim->Wells[WellID].Perforations[PerfIDMine].h3 = InBuffer[i];i++;

      MySim->Wells[WellID].Perforations[PerfIDMine].Pcow = InBuffer[i];i++;
      MySim->Wells[WellID].Perforations[PerfIDMine].Pcog = InBuffer[i];i++;
      MySim->Wells[WellID].Perforations[PerfIDMine].WellIndex = InBuffer[i];
       PetscPrintf(PETSC_COMM_WORLD, "\n At well: %d Value of WellIndex is:%g", WellID,MySim->Wells[WellID].Perforations[PerfIDMine].WellIndex);
      

    }
  }

  PetscFunctionReturn(0);
}

/* Compute production variables, update Observations vector */
extern PetscErrorCode DefiantBlackOil2PhProduction(BlackOilReservoirSimulation* MySim)
{
  PetscErrorCode ierr;
  PetscInt mx, my, mz, xm, ym, zm, xs, ys, zs;
  /* Well handling variables */
  PetscInt PerfIDMine, WellID;

  PetscFunctionBegin;
  /* Get dimensions and extents of the local vectors */
  ierr = DAGetInfo(MySim->SimDA, 0, &mx, &my, &mz, 0, 0, 0, 0, 0, 0, 0);CHKERRQ(ierr);
  ierr = DAGetCorners(MySim->SimDA, &xs, &ys, &zs, &xm, &ym, &zm);CHKERRQ(ierr);
  /* Sync everything before we start */
  ierr = DefiantBlackOilComputePerfIndicesSyncPerfs(MySim);CHKERRQ(ierr);

  for (WellID = 0; WellID < MySim->NumberOfWells; WellID++) {
    for (PerfIDMine = 0; PerfIDMine
        < (MySim->Wells[WellID]).NumberOfPerforations; PerfIDMine++) {
      if (MySim->Wells[WellID].Perforations[PerfIDMine].IsActive == PETSC_TRUE ){
        if (MySim->Wells[WellID].Perforations[PerfIDMine].Constraint == FLOW_RATE_CONSTRAINT) {
          MySim->Wells[WellID].Perforations[PerfIDMine].BHPo =
              MySim->Wells[WellID].Perforations[PerfIDMine].Bo
            / (MySim->Wells[WellID].Perforations[PerfIDMine].h1
            * MySim->Wells[WellID].Perforations[PerfIDMine].h2
            * MySim->Wells[WellID].Perforations[PerfIDMine].h3)
            * MySim->Wells[WellID].Perforations[PerfIDMine].Muo
            / MySim->Wells[WellID].Perforations[PerfIDMine].Kro
            / MySim->Wells[WellID].Perforations[PerfIDMine].WellIndex
            + MySim->Wells[WellID].Perforations[PerfIDMine].Po
            + MySim->GravAcc * MySim->Wells[WellID].Perforations[PerfIDMine].Rhoo
            * (MySim->Wells[WellID].Perforations[PerfIDMine].zbh
            - MySim->Wells[WellID].Perforations[PerfIDMine].x3);
          MySim->Wells[WellID].Perforations[PerfIDMine].BHPw =
              MySim->Wells[WellID].Perforations[PerfIDMine].Bw
            / (MySim->Wells[WellID].Perforations[PerfIDMine].h1
            * MySim->Wells[WellID].Perforations[PerfIDMine].h2
            * MySim->Wells[WellID].Perforations[PerfIDMine].h3)
            * MySim->Wells[WellID].Perforations[PerfIDMine].Muw
            / MySim->Wells[WellID].Perforations[PerfIDMine].Krw
            / MySim->Wells[WellID].Perforations[PerfIDMine].WellIndex
            + MySim->Wells[WellID].Perforations[PerfIDMine].Pw
            + MySim->GravAcc * MySim->Wells[WellID].Perforations[PerfIDMine].Rhow
            * (MySim->Wells[WellID].Perforations[PerfIDMine].zbh
            - MySim->Wells[WellID].Perforations[PerfIDMine].x3);
          PetscPrintf(PETSC_COMM_WORLD,"\nAt Well: %d Value of computed Qw is: %g",WellID,MySim->Wells[WellID].Perforations[PerfIDMine].Qw);
        } else if (MySim->Wells[WellID].Perforations[PerfIDMine].Constraint == BHP_CONSTRAINT) {
          MySim->Wells[WellID].Perforations[PerfIDMine].Qo =
                MySim->Wells[WellID].Perforations[PerfIDMine].WellIndex
              * MySim->Wells[WellID].Perforations[PerfIDMine].Kro
              / MySim->Wells[WellID].Perforations[PerfIDMine].Muo
              * (MySim->Wells[WellID].Perforations[PerfIDMine].BHPo
              - MySim->Wells[WellID].Perforations[PerfIDMine].Po
              - MySim->GravAcc * MySim->Wells[WellID].Perforations[PerfIDMine].Rhoo
              * (MySim->Wells[WellID].Perforations[PerfIDMine].zbh
              - MySim->Wells[WellID].Perforations[PerfIDMine].x3))
              * MySim->Wells[WellID].Perforations[PerfIDMine].h1
              * MySim->Wells[WellID].Perforations[PerfIDMine].h2
              * MySim->Wells[WellID].Perforations[PerfIDMine].h3
              / MySim->Wells[WellID].Perforations[PerfIDMine].Bo;
          MySim->Wells[WellID].Perforations[PerfIDMine].Qw =
                MySim->Wells[WellID].Perforations[PerfIDMine].WellIndex
              * MySim->Wells[WellID].Perforations[PerfIDMine].Krw
              / MySim->Wells[WellID].Perforations[PerfIDMine].Muw
              * (MySim->Wells[WellID].Perforations[PerfIDMine].BHPw
              - MySim->Wells[WellID].Perforations[PerfIDMine].Pw
              - MySim->GravAcc * MySim->Wells[WellID].Perforations[PerfIDMine].Rhow
              * (MySim->Wells[WellID].Perforations[PerfIDMine].zbh
              - MySim->Wells[WellID].Perforations[PerfIDMine].x3))
              * MySim->Wells[WellID].Perforations[PerfIDMine].h1
              * MySim->Wells[WellID].Perforations[PerfIDMine].h2
              * MySim->Wells[WellID].Perforations[PerfIDMine].h3
              / MySim->Wells[WellID].Perforations[PerfIDMine].Bw;

              PetscPrintf(PETSC_COMM_WORLD,"\nAt Well:%d Value of computed Qo is: %g",WellID, MySim->Wells[WellID].Perforations[PerfIDMine].Qo);
        }
      }
    }
  }


  PetscFunctionReturn(0);
}

/* Compute production variables, update Observations vector */
extern PetscErrorCode DefiantBlackOil3PhProduction(BlackOilReservoirSimulation* MySim)
{
  PetscErrorCode ierr;
  PetscInt mx, my, mz, xm, ym, zm, xs, ys, zs;
  /* Well handling variables */
  PetscInt PerfIDMine, WellID;

  PetscFunctionBegin;
  /* Get dimensions and extents of the local vectors */
  ierr = DAGetInfo(MySim->SimDA, 0, &mx, &my, &mz, 0, 0, 0, 0, 0, 0, 0);CHKERRQ(ierr);
  ierr = DAGetCorners(MySim->SimDA, &xs, &ys, &zs, &xm, &ym, &zm);CHKERRQ(ierr);
  /* Sync everything before we start */
  ierr = DefiantBlackOilComputePerfIndicesSyncPerfs(MySim);CHKERRQ(ierr);

  for (WellID = 0; WellID < MySim->NumberOfWells; WellID++) {
    for (PerfIDMine = 0; PerfIDMine
        < (MySim->Wells[WellID]).NumberOfPerforations; PerfIDMine++) {
      if (MySim->Wells[WellID].Perforations[PerfIDMine].IsActive == PETSC_TRUE ){
        if (MySim->Wells[WellID].Perforations[PerfIDMine].Constraint == FLOW_RATE_CONSTRAINT) {
          MySim->Wells[WellID].Perforations[PerfIDMine].BHPo =
              MySim->Wells[WellID].Perforations[PerfIDMine].Bo
            / (MySim->Wells[WellID].Perforations[PerfIDMine].h1
            * MySim->Wells[WellID].Perforations[PerfIDMine].h2
            * MySim->Wells[WellID].Perforations[PerfIDMine].h3)
            * MySim->Wells[WellID].Perforations[PerfIDMine].Muo
            / MySim->Wells[WellID].Perforations[PerfIDMine].Kro
            / MySim->Wells[WellID].Perforations[PerfIDMine].WellIndex
            + MySim->Wells[WellID].Perforations[PerfIDMine].Po
            + MySim->GravAcc * MySim->Wells[WellID].Perforations[PerfIDMine].Rhoo
            * (MySim->Wells[WellID].Perforations[PerfIDMine].zbh
            - MySim->Wells[WellID].Perforations[PerfIDMine].x3);
          MySim->Wells[WellID].Perforations[PerfIDMine].BHPw =
              MySim->Wells[WellID].Perforations[PerfIDMine].Bw
            / (MySim->Wells[WellID].Perforations[PerfIDMine].h1
            * MySim->Wells[WellID].Perforations[PerfIDMine].h2
            * MySim->Wells[WellID].Perforations[PerfIDMine].h3)
            * MySim->Wells[WellID].Perforations[PerfIDMine].Muw
            / MySim->Wells[WellID].Perforations[PerfIDMine].Krw
            / MySim->Wells[WellID].Perforations[PerfIDMine].WellIndex
            + MySim->Wells[WellID].Perforations[PerfIDMine].Pw
            + MySim->GravAcc * MySim->Wells[WellID].Perforations[PerfIDMine].Rhow
            * (MySim->Wells[WellID].Perforations[PerfIDMine].zbh
            - MySim->Wells[WellID].Perforations[PerfIDMine].x3);
          MySim->Wells[WellID].Perforations[PerfIDMine].BHPg =
              MySim->Wells[WellID].Perforations[PerfIDMine].Bg
            / (MySim->Wells[WellID].Perforations[PerfIDMine].h1
            * MySim->Wells[WellID].Perforations[PerfIDMine].h2
            * MySim->Wells[WellID].Perforations[PerfIDMine].h3)
            * MySim->Wells[WellID].Perforations[PerfIDMine].Mug
            / MySim->Wells[WellID].Perforations[PerfIDMine].Krg
            / MySim->Wells[WellID].Perforations[PerfIDMine].WellIndex
            + MySim->Wells[WellID].Perforations[PerfIDMine].Pg
            + MySim->GravAcc * MySim->Wells[WellID].Perforations[PerfIDMine].Rhog
            * (MySim->Wells[WellID].Perforations[PerfIDMine].zbh
            - MySim->Wells[WellID].Perforations[PerfIDMine].x3);
        } else if (MySim->Wells[WellID].Perforations[PerfIDMine].Constraint == BHP_CONSTRAINT) {
          MySim->Wells[WellID].Perforations[PerfIDMine].Qo =
                MySim->Wells[WellID].Perforations[PerfIDMine].WellIndex
              * MySim->Wells[WellID].Perforations[PerfIDMine].Kro
              / MySim->Wells[WellID].Perforations[PerfIDMine].Muo
              * (MySim->Wells[WellID].Perforations[PerfIDMine].BHPo
              - MySim->Wells[WellID].Perforations[PerfIDMine].Po
              - MySim->GravAcc * MySim->Wells[WellID].Perforations[PerfIDMine].Rhoo
              * (MySim->Wells[WellID].Perforations[PerfIDMine].zbh
              - MySim->Wells[WellID].Perforations[PerfIDMine].x3))
              * MySim->Wells[WellID].Perforations[PerfIDMine].h1
              * MySim->Wells[WellID].Perforations[PerfIDMine].h2
              * MySim->Wells[WellID].Perforations[PerfIDMine].h3
              / MySim->Wells[WellID].Perforations[PerfIDMine].Bo;
          MySim->Wells[WellID].Perforations[PerfIDMine].Qw =
                MySim->Wells[WellID].Perforations[PerfIDMine].WellIndex
              * MySim->Wells[WellID].Perforations[PerfIDMine].Krw
              / MySim->Wells[WellID].Perforations[PerfIDMine].Muw
              * (MySim->Wells[WellID].Perforations[PerfIDMine].BHPw
              - MySim->Wells[WellID].Perforations[PerfIDMine].Pw
              - MySim->GravAcc * MySim->Wells[WellID].Perforations[PerfIDMine].Rhow
              * (MySim->Wells[WellID].Perforations[PerfIDMine].zbh
              - MySim->Wells[WellID].Perforations[PerfIDMine].x3))
              * MySim->Wells[WellID].Perforations[PerfIDMine].h1
              * MySim->Wells[WellID].Perforations[PerfIDMine].h2
              * MySim->Wells[WellID].Perforations[PerfIDMine].h3
              / MySim->Wells[WellID].Perforations[PerfIDMine].Bw;
          MySim->Wells[WellID].Perforations[PerfIDMine].Qg =
                MySim->Wells[WellID].Perforations[PerfIDMine].WellIndex
              * MySim->Wells[WellID].Perforations[PerfIDMine].Krg
              / MySim->Wells[WellID].Perforations[PerfIDMine].Mug
              * (MySim->Wells[WellID].Perforations[PerfIDMine].BHPg
              - MySim->Wells[WellID].Perforations[PerfIDMine].Pg
              - MySim->GravAcc * MySim->Wells[WellID].Perforations[PerfIDMine].Rhog
              * (MySim->Wells[WellID].Perforations[PerfIDMine].zbh
              - MySim->Wells[WellID].Perforations[PerfIDMine].x3))
              * MySim->Wells[WellID].Perforations[PerfIDMine].h1
              * MySim->Wells[WellID].Perforations[PerfIDMine].h2
              * MySim->Wells[WellID].Perforations[PerfIDMine].h3
              / MySim->Wells[WellID].Perforations[PerfIDMine].Bg;
        }
      }
    }
  }


  PetscFunctionReturn(0);
}
