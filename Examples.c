/*
 * Examples.c
 *
 *  Created on: Sep 13, 2009
 *      Author: yye00
 */

#include "Defiant.h"

#undef __FUNCT__
#define __FUNCT__ "DefiantIMPES2PhBuckleyLeverett"
extern PetscErrorCode DefiantIMPES2PhBuckleyLeverett()
{
  PetscErrorCode ierr;

  BlackOilReservoirSimulation BuckleyLeverett;

  PetscFunctionBegin;
  /* Initialize the simulation */
  ierr = DefiantCreateSimulationVecs(&BuckleyLeverett);CHKERRQ(ierr);CHKMEMQ;
  ierr = DefiantCreate2PhIMPESSystem(&BuckleyLeverett);CHKERRQ(ierr);CHKMEMQ;

  /* Set FlowMask */
  ierr = VecSet(BuckleyLeverett.FlowMask, FLUID_FLOW);CHKERRQ(ierr);CHKMEMQ;

  /* Set the coordinates */
  ierr = DASetUniformCoordinates(BuckleyLeverett.SimDA, 0.0, 90.0, 0.0, 2.0, 0.0, 2.0);CHKERRQ(ierr);CHKMEMQ;
  ierr = DefiantGetDACoords(&BuckleyLeverett);CHKERRQ(ierr);CHKMEMQ;
  /* set some geometry */
  ierr = VecSet(BuckleyLeverett.A1, 1.0);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSet(BuckleyLeverett.A2, 10.0);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSet(BuckleyLeverett.A3, 10.0);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSet(BuckleyLeverett.h1, 10.0);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSet(BuckleyLeverett.h2, 1.0);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSet(BuckleyLeverett.h3, 1.0);CHKERRQ(ierr);CHKMEMQ;
  /* set initial values to vectors */
  ierr = VecSet(BuckleyLeverett.Po, 4000.0);CHKERRQ(ierr);CHKMEMQ;
  /* set residual oil and connate water */
  BuckleyLeverett.Sor = 0.2;
  BuckleyLeverett.Swc = 0.2;
  /* Set saturations */
  ierr = VecSet(BuckleyLeverett.So, 1.0-0.2 );CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSet(BuckleyLeverett.Sw, 0.2 );CHKERRQ(ierr);CHKMEMQ;
  /* Set densities */
  BuckleyLeverett.Rhoos = 1.0;
  BuckleyLeverett.Rhows = 1.0;
  ierr = VecSet(BuckleyLeverett.Rhoo, BuckleyLeverett.Rhoos );CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSet(BuckleyLeverett.Rhow, BuckleyLeverett.Rhows );CHKERRQ(ierr);CHKMEMQ;
  /* Set viscosities */
  ierr = VecSet(BuckleyLeverett.Muo, 1);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSet(BuckleyLeverett.Muw, 1);CHKERRQ(ierr);CHKMEMQ;
  /* Set permeabilities */
  ierr = VecSet(BuckleyLeverett.K11, 300.0);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSet(BuckleyLeverett.K22, 300.0);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSet(BuckleyLeverett.K33, 300.0);CHKERRQ(ierr);CHKMEMQ;
  /* set the relative permeabilities */
  ierr = VecSet(BuckleyLeverett.Kro, 1.0);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSet(BuckleyLeverett.Krw, 1.0);CHKERRQ(ierr);CHKMEMQ;
  /* Set the table for PcowSw */
  BuckleyLeverett.PcowSw.NumberOfEntries = 2;
  ierr = PetscMalloc(BuckleyLeverett.PcowSw.NumberOfEntries*sizeof(PetscScalar),&BuckleyLeverett.PcowSw.X);CHKERRQ(ierr);CHKMEMQ;
  ierr = PetscMalloc(BuckleyLeverett.PcowSw.NumberOfEntries*sizeof(PetscScalar),&BuckleyLeverett.PcowSw.Y);CHKERRQ(ierr);CHKMEMQ;
  BuckleyLeverett.PcowSw.X[0] = 0;BuckleyLeverett.PcowSw.X[1] = 0;
  BuckleyLeverett.PcowSw.Y[0] = 0;BuckleyLeverett.PcowSw.Y[1] = 1;
  /* Set the graviational acceleration */
  BuckleyLeverett.GravAcc = 32.2;
  /* Set the rock properties */
  ierr = VecSet(BuckleyLeverett.Phi, 0.2);CHKERRQ(ierr);CHKMEMQ;
  BuckleyLeverett.RockCompressibility = 1.0e-15;
  BuckleyLeverett.RockCompRefPressure = 14.7;
  BuckleyLeverett.RockCompRefPorosity = 0.3;
  /* Set the fluid compressibilities */
  BuckleyLeverett.OilCompressibility = 1.0e-5;
  BuckleyLeverett.OilCompPressure = 14.7;
  BuckleyLeverett.WaterCompressibility = 1.0e-5;
  BuckleyLeverett.WaterCompPressure = 14.7;
  /* FIXME Set the volume factors */
  ierr = VecSet(BuckleyLeverett.Bo, 1.0);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSet(BuckleyLeverett.Bw, 1.0);CHKERRQ(ierr);CHKMEMQ;
  /* Set the well information */
  /* allocate wells */
  BuckleyLeverett.NumberOfWells = 2;
  ierr = PetscMalloc(BuckleyLeverett.NumberOfWells*sizeof(Well),&BuckleyLeverett.Wells);CHKERRQ(ierr);CHKMEMQ;
  /* allocate perforations for wells */
  BuckleyLeverett.Wells[0].NumberOfPerforations = 1;
  ierr = PetscMalloc(BuckleyLeverett.Wells[0].NumberOfPerforations*sizeof(Perforation),&BuckleyLeverett.Wells[0].Perforations);CHKERRQ(ierr);CHKMEMQ;
  BuckleyLeverett.Wells[1].NumberOfPerforations = 1;
  ierr = PetscMalloc(BuckleyLeverett.Wells[1].NumberOfPerforations*sizeof(Perforation),&BuckleyLeverett.Wells[1].Perforations);CHKERRQ(ierr);CHKMEMQ;
  /* Set information for well zero,the injector */
  BuckleyLeverett.Wells[0].Perforations[0].I = 1;
  BuckleyLeverett.Wells[0].Perforations[0].J = 1;
  BuckleyLeverett.Wells[0].Perforations[0].K = 1;
  BuckleyLeverett.Wells[0].Perforations[0].Qw = 1.0;
  BuckleyLeverett.Wells[0].Perforations[0].Constraint = FLOW_RATE_CONSTRAINT;
  BuckleyLeverett.Wells[0].Perforations[0].WellType = WATER_INJECTOR;
  BuckleyLeverett.Wells[0].Perforations[0].IsActive = PETSC_TRUE;
  BuckleyLeverett.Wells[0].Perforations[0].Orientation = PERF_ORIENTATION_X3X3;
  BuckleyLeverett.Wells[0].Perforations[0].Rw = 0.3;
  BuckleyLeverett.Wells[0].Perforations[0].S = 0.25;
  BuckleyLeverett.Wells[0].Perforations[0].zbh = 0.0;

  /* Set information for well zero,the producer */
  BuckleyLeverett.Wells[1].Perforations[0].I = 8;
  BuckleyLeverett.Wells[1].Perforations[0].J = 1;
  BuckleyLeverett.Wells[1].Perforations[0].K = 1;
  BuckleyLeverett.Wells[1].Perforations[0].BHPw = 3999;
  BuckleyLeverett.Wells[1].Perforations[0].BHPo = 3999;
  BuckleyLeverett.Wells[1].Perforations[0].Constraint = BHP_CONSTRAINT;
  BuckleyLeverett.Wells[1].Perforations[0].IsActive = PETSC_TRUE;
  BuckleyLeverett.Wells[1].Perforations[0].Orientation = PERF_ORIENTATION_X3X3;
  BuckleyLeverett.Wells[1].Perforations[0].Rw = 0.3;
  BuckleyLeverett.Wells[1].Perforations[0].S = 0.25;
  BuckleyLeverett.Wells[1].Perforations[0].zbh = 0.0;

  /* Not do adaptive time-step */
  BuckleyLeverett.AdaptiveTimeStep = PETSC_FALSE;

  /* Set the Times */
  BuckleyLeverett.DeltaTP = 0.5;
  BuckleyLeverett.DeltaTS = 0.1;
  BuckleyLeverett.DSmax = 0.1;

  /* Set the start and end time */
  BuckleyLeverett.StartTime = 0.0;
  BuckleyLeverett.EndTime = 7.0;

  /* Now we solve */
  ierr = DefiantIMPES2PhIterate(&BuckleyLeverett);CHKERRQ(ierr);CHKMEMQ;
  ierr = DefiantBlackOil2PhValiantWriteVecs(&BuckleyLeverett);CHKERRQ(ierr);CHKMEMQ;

  ierr = PetscPrintf(PETSC_COMM_WORLD,"Sw is: \n");
  ierr = VecView(BuckleyLeverett.Sw,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
  ierr = PetscPrintf(PETSC_COMM_WORLD,"x1 is: \n");
  ierr = VecView(BuckleyLeverett.x1,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
  ierr = PetscPrintf(PETSC_COMM_WORLD,"h1 is: \n");
  ierr = VecView(BuckleyLeverett.h1,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;

  /* Now we clean up */
  ierr = DefiantDestroySimulationVecs(&BuckleyLeverett);CHKERRQ(ierr);CHKMEMQ;
  ierr = DefiantDestroy2PhIMPESSystem(&BuckleyLeverett);CHKERRQ(ierr);CHKMEMQ;


  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DefiantIMPES3PhBuckleyLeverett"
extern PetscErrorCode DefiantIMPES3PhBuckleyLeverett();

#undef __FUNCT__
#define __FUNCT__ "DefiantImprovedIMPES2PhBuckleyLeverett"
extern PetscErrorCode DefiantImprovedIMPES2PhBuckleyLeverett();

#undef __FUNCT__
#define __FUNCT__ "DefiantImprovedIMPES3PhBuckleyLeverett"
extern PetscErrorCode DefiantImprovedIMPES3PhBuckleyLeverett();

#undef __FUNCT__
#define __FUNCT__ "DefiantNewton2PhBuckleyLeverett"
extern PetscErrorCode DefiantNewton2PhBuckleyLeverett();

#undef __FUNCT__
#define __FUNCT__ "DefiantNewton3PhBuckleyLeverett"
extern PetscErrorCode DefiantNewton3PhBuckleyLeverett();

#undef __FUNCT__
#define __FUNCT__ "DefiantIMPES2PhFivePoint"
extern PetscErrorCode DefiantIMPES2PhFivePoint()
{
  PetscErrorCode ierr;
  PetscViewer viewer;

  BlackOilReservoirSimulation FivePoint;

  PetscFunctionBegin;
  /* Initialize the simulation */
  ierr = DefiantCreateSimulationVecs(&FivePoint);CHKERRQ(ierr);CHKMEMQ;
  ierr = DefiantCreate2PhIMPESSystem(&FivePoint);CHKERRQ(ierr);CHKMEMQ;

  /* Set FlowMask */
  ierr = VecSet(FivePoint.FlowMask, FLUID_FLOW);CHKERRQ(ierr);CHKMEMQ;

  /* Set the coordinates */
  ierr = DASetUniformCoordinates(FivePoint.SimDA, 0.0, 180.0, 0.0, 180.0, 0.0, 2.0);CHKERRQ(ierr);CHKMEMQ;
  ierr = DefiantGetDACoords(&FivePoint);CHKERRQ(ierr);CHKMEMQ;
  /* set some geometry */
  ierr = VecSet(FivePoint.A1, 10.0);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSet(FivePoint.A2, 10.0);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSet(FivePoint.A3, 1.0);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSet(FivePoint.h1, 10.0);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSet(FivePoint.h2, 10.0);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSet(FivePoint.h3, 1.0);CHKERRQ(ierr);CHKMEMQ;
  /* set initial values to vectors */
  ierr = VecSet(FivePoint.Po, 4000.0);CHKERRQ(ierr);CHKMEMQ;
  /* set residual oil and connate water */
  FivePoint.Sor = 0.2;
  FivePoint.Swc = 0.2;
  /* Set saturations */
  ierr = VecSet(FivePoint.So, 1.0-0.2 );CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSet(FivePoint.Sw, 0.2 );CHKERRQ(ierr);CHKMEMQ;
  /* Set densities */
  FivePoint.Rhoos = 1.0;
  FivePoint.Rhows = 1.0;
  ierr = VecSet(FivePoint.Rhoo, FivePoint.Rhoos );CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSet(FivePoint.Rhow, FivePoint.Rhows );CHKERRQ(ierr);CHKMEMQ;
  /* Set viscosities */
  ierr = VecSet(FivePoint.Muo, 1);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSet(FivePoint.Muw, 1);CHKERRQ(ierr);CHKMEMQ;
  /* Set permeabilities */
  ierr = VecSet(FivePoint.K11, 300.0);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSet(FivePoint.K22, 300.0);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSet(FivePoint.K33, 300.0);CHKERRQ(ierr);CHKMEMQ;
  /* set the relative permeabilities */
  ierr = VecSet(FivePoint.Kro, 1.0);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSet(FivePoint.Krw, 1.0);CHKERRQ(ierr);CHKMEMQ;
  /* Set the table for PcowSw */
  FivePoint.PcowSw.NumberOfEntries = 2;
  ierr = PetscMalloc(FivePoint.PcowSw.NumberOfEntries*sizeof(PetscScalar),&FivePoint.PcowSw.X);CHKERRQ(ierr);CHKMEMQ;
  ierr = PetscMalloc(FivePoint.PcowSw.NumberOfEntries*sizeof(PetscScalar),&FivePoint.PcowSw.Y);CHKERRQ(ierr);CHKMEMQ;
  FivePoint.PcowSw.X[0] = 0;FivePoint.PcowSw.X[1] = 0;
  FivePoint.PcowSw.Y[0] = 0;FivePoint.PcowSw.Y[1] = 1;
  /* Set the graviational acceleration */
  FivePoint.GravAcc = 32.2;
  /* Set the rock properties */
  ierr = VecSet(FivePoint.Phi, 0.2);CHKERRQ(ierr);CHKMEMQ;
  FivePoint.RockCompressibility = 1.0e-15;
  FivePoint.RockCompRefPressure = 14.7;
  FivePoint.RockCompRefPorosity = 0.3;
  /* Set the fluid compressibilities */
  FivePoint.OilCompressibility = 1.0e-5;
  FivePoint.OilCompPressure = 14.7;
  FivePoint.WaterCompressibility = 1.0e-5;
  FivePoint.WaterCompPressure = 14.7;
  /* FIXME Set the volume factors */
  ierr = VecSet(FivePoint.Bo, 1.0);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSet(FivePoint.Bw, 1.0);CHKERRQ(ierr);CHKMEMQ;
  /* Set the well information */
  /* allocate wells */
  FivePoint.NumberOfWells = 2;
  ierr = PetscMalloc(FivePoint.NumberOfWells*sizeof(Well),&FivePoint.Wells);CHKERRQ(ierr);CHKMEMQ;
  /* allocate perforations for wells */
  FivePoint.Wells[0].NumberOfPerforations = 1;
  ierr = PetscMalloc(FivePoint.Wells[0].NumberOfPerforations*sizeof(Perforation),&FivePoint.Wells[0].Perforations);CHKERRQ(ierr);CHKMEMQ;
  FivePoint.Wells[1].NumberOfPerforations = 1;
  ierr = PetscMalloc(FivePoint.Wells[1].NumberOfPerforations*sizeof(Perforation),&FivePoint.Wells[1].Perforations);CHKERRQ(ierr);CHKMEMQ;
  /* Set information for well zero,the injector */
  FivePoint.Wells[0].Perforations[0].I = 1;
  FivePoint.Wells[0].Perforations[0].J = 1;
  FivePoint.Wells[0].Perforations[0].K = 1;
  FivePoint.Wells[0].Perforations[0].Qw = 2.0;
  FivePoint.Wells[0].Perforations[0].Constraint = FLOW_RATE_CONSTRAINT;
  FivePoint.Wells[0].Perforations[0].WellType = WATER_INJECTOR;
  FivePoint.Wells[0].Perforations[0].IsActive = PETSC_TRUE;
  FivePoint.Wells[0].Perforations[0].Orientation = PERF_ORIENTATION_X3X3;
  FivePoint.Wells[0].Perforations[0].Rw = 0.3;
  FivePoint.Wells[0].Perforations[0].S = 0.25;
  FivePoint.Wells[0].Perforations[0].zbh = 0.0;

  /* Set information for well zero,the producer */
  FivePoint.Wells[1].Perforations[0].I = 18;
  FivePoint.Wells[1].Perforations[0].J = 18;
  FivePoint.Wells[1].Perforations[0].K = 1;
  FivePoint.Wells[1].Perforations[0].BHPw = 3999;
  FivePoint.Wells[1].Perforations[0].BHPo = 3999;
  FivePoint.Wells[1].Perforations[0].Constraint = BHP_CONSTRAINT;
  FivePoint.Wells[1].Perforations[0].IsActive = PETSC_TRUE;
  FivePoint.Wells[1].Perforations[0].Orientation = PERF_ORIENTATION_X3X3;
  FivePoint.Wells[1].Perforations[0].Rw = 0.3;
  FivePoint.Wells[1].Perforations[0].S = 0.25;
  FivePoint.Wells[1].Perforations[0].zbh = 0.0;

  /* Not do adaptive time-step */
  FivePoint.AdaptiveTimeStep = PETSC_FALSE;

  /* Set the Times */
  FivePoint.DeltaTP = 2.0;
  FivePoint.DeltaTS = 0.5;
  FivePoint.DSmax = 0.1;

  /* Set the start and end time */
  FivePoint.StartTime = 0.0;
  FivePoint.EndTime = 30.0;

  /* Now we solve */
  ierr = DefiantIMPES2PhIterate(&FivePoint);CHKERRQ(ierr);CHKMEMQ;
  ierr = DefiantBlackOil2PhValiantWriteVecs(&FivePoint);CHKERRQ(ierr);CHKMEMQ;

#if DEFIANT_DEBUG
  ierr = PetscViewerCreate(PETSC_COMM_WORLD, &viewer);CHKERRQ(ierr);CHKMEMQ;
  ierr = PetscViewerSetType(viewer, PETSC_VIEWER_ASCII);CHKERRQ(ierr);CHKMEMQ;
  ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_VTK);CHKERRQ(ierr);CHKMEMQ;
  ierr = PetscViewerFileSetName(viewer, "FivePoint.vtk");CHKERRQ(ierr);CHKMEMQ;
  ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);CHKMEMQ;
#endif

  ierr = PetscPrintf(PETSC_COMM_WORLD,"Sw is: \n");
  ierr = VecView(FivePoint.Sw,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
  ierr = PetscPrintf(PETSC_COMM_WORLD,"x1 is: \n");
  ierr = VecView(FivePoint.x1,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;

  /* Now we clean up */
  ierr = DefiantDestroySimulationVecs(&FivePoint);CHKERRQ(ierr);CHKMEMQ;
  ierr = DefiantDestroy2PhIMPESSystem(&FivePoint);CHKERRQ(ierr);CHKMEMQ;


  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DefiantIMPES3PhFivePoint"
extern PetscErrorCode DefiantIMPES3PhFivePoint(BlackOilReservoirSimulation* MySim);

#undef __FUNCT__
#define __FUNCT__ "DefiantImprovedIMPES2PhFivePoint"
extern PetscErrorCode DefiantImprovedIMPES2PhFivePoint(BlackOilReservoirSimulation* MySim);

#undef __FUNCT__
#define __FUNCT__ "DefiantImprovedIMPES3PhFivePoint"
extern PetscErrorCode DefiantImprovedIMPES3PhFivePoint(BlackOilReservoirSimulation* MySim);

#undef __FUNCT__
#define __FUNCT__ "DefiantNewton2PhFivePoint"
extern PetscErrorCode DefiantNewton2PhFivePoint(BlackOilReservoirSimulation* MySim);

#undef __FUNCT__
#define __FUNCT__ "DefiantNewton3PhFivePoint"
extern PetscErrorCode DefiantNewton3PhFivePoint(BlackOilReservoirSimulation* MySim);
