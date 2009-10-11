/*
 * ValiantInterface.c
 *
 *  Created on: Oct 8, 2009
 *      Author: yye00
 */

#include "Defiant.h"

#undef __FUNCT__
#define __FUNCT__ "DefiantBlackOil2PhValiantLoadVecs"
extern PetscErrorCode DefiantBlackOil2PhValiantLoadVecs(BlackOilReservoirSimulation* MySim)
{
  PetscErrorCode ierr;
  PetscViewer viewer;

  PetscFunctionBegin;

  /* Load the porosity */
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"Phi.Valiant",FILE_MODE_READ,&viewer);CHKERRQ(ierr);
  ierr = VecLoad(viewer,PETSC_NULL,&MySim->Phi);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);

  /* Load the LnK 11 22 33 */
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"LnK11.Valiant",FILE_MODE_READ,&viewer);CHKERRQ(ierr);
  ierr = VecLoad(viewer,PETSC_NULL,&MySim->K11);CHKERRQ(ierr);
  ierr = VecExp(MySim->K11);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"LnK22.Valiant",FILE_MODE_READ,&viewer);CHKERRQ(ierr);
  ierr = VecLoad(viewer,PETSC_NULL,&MySim->K22);CHKERRQ(ierr);
  ierr = VecExp(MySim->K22);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"LnK33.Valiant",FILE_MODE_READ,&viewer);CHKERRQ(ierr);
  ierr = VecLoad(viewer,PETSC_NULL,&MySim->K33);CHKERRQ(ierr);
  ierr = VecExp(MySim->K33);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);

  /* Load the pressure */
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"Pw.Valiant",FILE_MODE_READ,&viewer);CHKERRQ(ierr);
  ierr = VecLoad(viewer,PETSC_NULL,&MySim->Pw);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);

  /* Load the Saturation */
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"Sw.Valiant",FILE_MODE_READ,&viewer);CHKERRQ(ierr);
  ierr = VecLoad(viewer,PETSC_NULL,&MySim->Sw);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "DefiantBlackOil3PhValiantLoadVecs"
extern PetscErrorCode DefiantBlackOil3PhValiantLoadVecs(BlackOilReservoirSimulation* MySim)
{
  PetscErrorCode ierr;
  PetscViewer viewer;

  PetscFunctionBegin;

  /* Load the porosity */
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"Phi.Valiant",FILE_MODE_READ,&viewer);CHKERRQ(ierr);
  ierr = VecLoad(viewer,PETSC_NULL,&MySim->Phi);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);

  /* Load the LnK 11 22 33 */
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"LnK11.Valiant",FILE_MODE_READ,&viewer);CHKERRQ(ierr);
  ierr = VecLoad(viewer,PETSC_NULL,&MySim->K11);CHKERRQ(ierr);
  ierr = VecExp(MySim->K11);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"LnK22.Valiant",FILE_MODE_READ,&viewer);CHKERRQ(ierr);
  ierr = VecLoad(viewer,PETSC_NULL,&MySim->K22);CHKERRQ(ierr);
  ierr = VecExp(MySim->K22);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"LnK33.Valiant",FILE_MODE_READ,&viewer);CHKERRQ(ierr);
  ierr = VecLoad(viewer,PETSC_NULL,&MySim->K33);CHKERRQ(ierr);
  ierr = VecExp(MySim->K33);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);

  /* Load the pressure */
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"Pw.Valiant",FILE_MODE_READ,&viewer);CHKERRQ(ierr);
  ierr = VecLoad(viewer,PETSC_NULL,&MySim->Pw);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"Po.Valiant",FILE_MODE_READ,&viewer);CHKERRQ(ierr);
  ierr = VecLoad(viewer,PETSC_NULL,&MySim->Po);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"Pg.Valiant",FILE_MODE_READ,&viewer);CHKERRQ(ierr);
  ierr = VecLoad(viewer,PETSC_NULL,&MySim->Pg);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);

  /* Load the Saturation */
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"Sw.Valiant",FILE_MODE_READ,&viewer);CHKERRQ(ierr);
  ierr = VecLoad(viewer,PETSC_NULL,&MySim->Sw);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"So.Valiant",FILE_MODE_READ,&viewer);CHKERRQ(ierr);
  ierr = VecLoad(viewer,PETSC_NULL,&MySim->Sw);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DefiantBlackOil2PhValiantWriteVecs"
extern PetscErrorCode DefiantBlackOil2PhValiantWriteVecs(BlackOilReservoirSimulation* MySim)
{
  PetscErrorCode ierr;
  PetscInt i, TempInt, ObsSize, WellID, PerfID;
  Vec ObsVec;
  PetscViewer viewer;
  Vec TempVec;
  PetscFunctionBegin;
  /* create the tempvec */
  ierr = VecDuplicate(MySim->Po,&TempVec);CHKERRQ(ierr);CHKMEMQ;

  /* Write the porosity */
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"Phi.Valiant",FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
  ierr = VecView(MySim->Phi,viewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);

   /* Write the LnK 11 22 33 */
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"LnK11.Valiant",FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
  ierr = VecCopy(MySim->K11, TempVec);CHKERRQ(ierr);
  ierr = VecLog(TempVec);CHKERRQ(ierr);
  ierr = VecView(TempVec,viewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"LnK22.Valiant",FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
  ierr = VecCopy(MySim->K22, TempVec);CHKERRQ(ierr);
  ierr = VecLog(TempVec);CHKERRQ(ierr);
  ierr = VecView(TempVec,viewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"LnK33.Valiant",FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
  ierr = VecCopy(MySim->K33, TempVec);CHKERRQ(ierr);
  ierr = VecLog(TempVec);CHKERRQ(ierr);
  ierr = VecView(TempVec,viewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);

  ierr = VecDestroy(TempVec);CHKERRQ(ierr);

   /* Write the pressure */
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"Pw.Valiant",FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
  ierr = VecView(MySim->Pw,viewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);

  /* Write the Saturation */
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"Sw.Valiant",FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
  ierr = VecView(MySim->Sw,viewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);

  /* Now for the actual observations */
  /* Make sure we are sync'ed before we write the observations */
  ierr = DefiantBlackOilComputePerfIndicesSyncPerfs(MySim);CHKERRQ(ierr);CHKMEMQ;
  /* Find out how many perforations are there */
  i=0;
  for (WellID=0;WellID<MySim->NumberOfWells;WellID++)
    for (PerfID=0;PerfID<MySim->Wells[WellID].NumberOfPerforations;PerfID++)
      i++;
  ObsSize = i*2*2; /* 2 phases, flow rate and Bhp per phase */

  /* Now we have everything on proc zero, it outputs the observations vector */
  ierr = VecCreate(PETSC_COMM_WORLD, &ObsVec);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSetSizes(ObsVec,PETSC_DECIDE,ObsSize);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSetFromOptions(ObsVec);CHKERRQ(ierr);CHKMEMQ;
  i=0;
  for (WellID=0;WellID<MySim->NumberOfWells;WellID++)
    for (PerfID=0;PerfID<MySim->Wells[WellID].NumberOfPerforations;PerfID++)
    {
      TempInt = i;
      ierr = VecSetValues(ObsVec, 1, &TempInt, &MySim->Wells[WellID].Perforations[PerfID].Qo, INSERT_VALUES);CHKERRQ(ierr);CHKMEMQ;
      TempInt = i + 1;
      ierr = VecSetValues(ObsVec, 1, &TempInt, &MySim->Wells[WellID].Perforations[PerfID].Qw, INSERT_VALUES);CHKERRQ(ierr);CHKMEMQ;
      TempInt = i + 2;
      ierr = VecSetValues(ObsVec, 1, &TempInt, &MySim->Wells[WellID].Perforations[PerfID].BHPo, INSERT_VALUES);CHKERRQ(ierr);CHKMEMQ;
      TempInt = i + 3;
      ierr = VecSetValues(ObsVec, 1, &TempInt, &MySim->Wells[WellID].Perforations[PerfID].BHPw, INSERT_VALUES);CHKERRQ(ierr);CHKMEMQ;
      i = i + 4;
    }

  ierr = VecAssemblyBegin(ObsVec);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecAssemblyEnd(ObsVec);CHKERRQ(ierr);CHKMEMQ;

  /* Write the Observations */
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"Observations.Valiant",FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecView(ObsVec,viewer);CHKERRQ(ierr);CHKMEMQ;
  /* Destroy the viewer */
  ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);CHKMEMQ;

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DefiantBlackOil3PhValiantWriteVecs"
extern PetscErrorCode DefiantBlackOil3PhValiantWriteVecs(BlackOilReservoirSimulation* MySim)
{
  PetscErrorCode ierr;
  PetscInt i, TempInt, ObsSize, WellID, PerfID;
  Vec ObsVec;
  PetscViewer viewer;
  Vec TempVec;
  PetscFunctionBegin;

  /* create the tempvec */
  ierr = VecDuplicate(MySim->Po,&TempVec);CHKERRQ(ierr);CHKMEMQ;

  /* Write the porosity */
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"Phi.Valiant",FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
  ierr = VecView(MySim->Phi,viewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);

   /* Write the LnK 11 22 33 */
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"LnK11.Valiant",FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
  ierr = VecCopy(MySim->K11, TempVec);CHKERRQ(ierr);
  ierr = VecLog(TempVec);CHKERRQ(ierr);
  ierr = VecView(TempVec,viewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"LnK22.Valiant",FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
  ierr = VecCopy(MySim->K22, TempVec);CHKERRQ(ierr);
  ierr = VecLog(TempVec);CHKERRQ(ierr);
  ierr = VecView(TempVec,viewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"LnK33.Valiant",FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
  ierr = VecCopy(MySim->K33, TempVec);CHKERRQ(ierr);
  ierr = VecLog(TempVec);CHKERRQ(ierr);
  ierr = VecView(TempVec,viewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);

  ierr = VecDestroy(TempVec);CHKERRQ(ierr);

   /* Write the pressures */
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"Pw.Valiant",FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
  ierr = VecView(MySim->Pw,viewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"Po.Valiant",FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
  ierr = VecView(MySim->Po,viewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"Pg.Valiant",FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
  ierr = VecView(MySim->Pg,viewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);

  /* Write the Saturations */
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"Sw.Valiant",FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
  ierr = VecView(MySim->Sw,viewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"Sg.Valiant",FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
  ierr = VecView(MySim->Sg,viewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);

  /* Now for the actual observations */
  /* Make sure we are sync'ed before we write the observations */
  ierr = DefiantBlackOilComputePerfIndicesSyncPerfs(MySim);CHKERRQ(ierr);CHKMEMQ;
  /* Find out how many perforations are there */
  i=0;
  for (WellID=0;WellID<MySim->NumberOfWells;WellID++)
    for (PerfID=0;PerfID<MySim->Wells[WellID].NumberOfPerforations;PerfID++)
      i++;
  ObsSize = i*3*2; /* 3 phases, flow rate and Bhp per phase */

  /* Now we have everything on proc zero, it outputs the observations vector */
  ierr = VecCreate(PETSC_COMM_WORLD, &ObsVec);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSetSizes(ObsVec,PETSC_DECIDE,ObsSize);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSetFromOptions(ObsVec);CHKERRQ(ierr);CHKMEMQ;
  i=0;
  for (WellID=0;WellID<MySim->NumberOfWells;WellID++)
    for (PerfID=0;PerfID<MySim->Wells[WellID].NumberOfPerforations;PerfID++)
    {
      TempInt = i;
      ierr = VecSetValues(ObsVec, 1, &TempInt, &MySim->Wells[WellID].Perforations[PerfID].Qo, INSERT_VALUES);CHKERRQ(ierr);CHKMEMQ;
      TempInt = i + 1;
      ierr = VecSetValues(ObsVec, 1, &TempInt, &MySim->Wells[WellID].Perforations[PerfID].Qw, INSERT_VALUES);CHKERRQ(ierr);CHKMEMQ;
      TempInt = i + 2;
      ierr = VecSetValues(ObsVec, 1, &TempInt, &MySim->Wells[WellID].Perforations[PerfID].Qg, INSERT_VALUES);CHKERRQ(ierr);CHKMEMQ;
      TempInt = i + 3;
      ierr = VecSetValues(ObsVec, 1, &TempInt, &MySim->Wells[WellID].Perforations[PerfID].BHPo, INSERT_VALUES);CHKERRQ(ierr);CHKMEMQ;
      TempInt = i + 4;
      ierr = VecSetValues(ObsVec, 1, &TempInt, &MySim->Wells[WellID].Perforations[PerfID].BHPw, INSERT_VALUES);CHKERRQ(ierr);CHKMEMQ;
      TempInt = i + 5;
      ierr = VecSetValues(ObsVec, 1, &TempInt, &MySim->Wells[WellID].Perforations[PerfID].BHPg, INSERT_VALUES);CHKERRQ(ierr);CHKMEMQ;
      i = i + 6;
    }

  ierr = VecAssemblyBegin(ObsVec);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecAssemblyEnd(ObsVec);CHKERRQ(ierr);CHKMEMQ;

  /* Write the Observations */
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"Observations.Valiant",FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecView(ObsVec,viewer);CHKERRQ(ierr);CHKMEMQ;
  /* Destroy the viewer */
  ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);CHKMEMQ;

  PetscFunctionReturn(0);
}
