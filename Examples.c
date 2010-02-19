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
  PetscInt mx, my, mz, xm, ym, zm, xs, ys, zs;
  PetscScalar ***Localx1;
  PetscScalar TempScalar;

  BlackOilReservoirSimulation BuckleyLeverett;

  PetscFunctionBegin;
  /* Initialize the simulation */
  ierr = DefiantCreateSimulationVecs(&BuckleyLeverett);CHKERRQ(ierr);CHKMEMQ;
  ierr = DefiantCreate2PhIMPESSystem(&BuckleyLeverett);CHKERRQ(ierr);CHKMEMQ;

  /* Set FlowMask */
  ierr = VecSet(BuckleyLeverett.FlowMask, FLUID_FLOW);CHKERRQ(ierr);CHKMEMQ;

  /* Set the coordinates */
  ierr = DAGetInfo(BuckleyLeverett.SimDA, 0, &mx, &my, &mz, 0, 0, 0, 0, 0, 0, 0);CHKERRQ(ierr);
  ierr = DAGetCorners(BuckleyLeverett.SimDA, &xs, &ys, &zs, &xm, &ym, &zm);CHKERRQ(ierr);
  ierr = DASetUniformCoordinates(BuckleyLeverett.SimDA, -1000.0/(mx-3), 1000.0+1000.0/(mx-3), 0.0, 2.0, 0.0, 2.0);CHKERRQ(ierr);CHKMEMQ;
  ierr = DefiantGetDACoords(&BuckleyLeverett);CHKERRQ(ierr);CHKMEMQ;
  /* set some geometry */
  ierr = VecSet(BuckleyLeverett.A1, 50000.0);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSet(BuckleyLeverett.A2, 1000.0/(mx-3)*100);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSet(BuckleyLeverett.A3, 1000.0/(mx-3)*500);CHKERRQ(ierr);CHKMEMQ;

  ierr = DAVecGetArray(BuckleyLeverett.SimDA,BuckleyLeverett.x1, &Localx1);CHKERRQ(ierr);CHKMEMQ;
  TempScalar = Localx1[1][1][1]-Localx1[1][1][0];
  ierr = PetscPrintf(PETSC_COMM_WORLD,"We Set h1 to: %g\n", TempScalar);


  ierr = VecSet(BuckleyLeverett.h1, TempScalar);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSet(BuckleyLeverett.h2, 500.0);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSet(BuckleyLeverett.h3, 100.0);CHKERRQ(ierr);CHKMEMQ;
  /* set initial values to vectors */
  ierr = VecSet(BuckleyLeverett.Po, 4000.0);CHKERRQ(ierr);CHKMEMQ;
  /* set residual oil and connate water */
  BuckleyLeverett.Sor = 0.1;
  BuckleyLeverett.Swc = 0.2;
  /* Set saturations */
  ierr = VecSet(BuckleyLeverett.So, 1.0-0.2);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSet(BuckleyLeverett.Sw, 0.2 );CHKERRQ(ierr);CHKMEMQ;
  /* Set densities */
  BuckleyLeverett.Rhoos = 1.0;
  BuckleyLeverett.Rhows = 1.0;
  ierr = VecSet(BuckleyLeverett.Rhoo, BuckleyLeverett.Rhoos );CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSet(BuckleyLeverett.Rhow, BuckleyLeverett.Rhows );CHKERRQ(ierr);CHKMEMQ;
  /* Set viscosities */
  ierr = VecSet(BuckleyLeverett.Muo, 6);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSet(BuckleyLeverett.Muw, 0.8);CHKERRQ(ierr);CHKMEMQ;
  /* Set permeabilities */
  ierr = VecSet(BuckleyLeverett.K11, 100.0);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSet(BuckleyLeverett.K22, 100.0);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSet(BuckleyLeverett.K33, 100.0);CHKERRQ(ierr);CHKMEMQ;
  /* set the relative permeabilities */
  ierr = VecSet(BuckleyLeverett.Kro, 1.0);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSet(BuckleyLeverett.Krw, 0.0);CHKERRQ(ierr);CHKMEMQ;
  /* Set the table for PcowSw */
  BuckleyLeverett.PcowSw.NumberOfEntries = 2;
  ierr = PetscMalloc(BuckleyLeverett.PcowSw.NumberOfEntries*sizeof(PetscScalar),&BuckleyLeverett.PcowSw.X);CHKERRQ(ierr);CHKMEMQ;
  ierr = PetscMalloc(BuckleyLeverett.PcowSw.NumberOfEntries*sizeof(PetscScalar),&BuckleyLeverett.PcowSw.Y);CHKERRQ(ierr);CHKMEMQ;
  BuckleyLeverett.PcowSw.X[0] = 0;BuckleyLeverett.PcowSw.X[1] = 0;
  BuckleyLeverett.PcowSw.Y[0] = 0;BuckleyLeverett.PcowSw.Y[1] = 1;
  /* Set the graviational acceleration */
  BuckleyLeverett.GravAcc = 0.0;
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
  BuckleyLeverett.Wells[0].Perforations[0].Qw = 100/0.001127;
  BuckleyLeverett.Wells[0].Perforations[0].Constraint = FLOW_RATE_CONSTRAINT;
  BuckleyLeverett.Wells[0].Perforations[0].WellType = WATER_INJECTOR;
  BuckleyLeverett.Wells[0].Perforations[0].IsActive = PETSC_TRUE;
  BuckleyLeverett.Wells[0].Perforations[0].Orientation = PERF_ORIENTATION_X3X3;
  BuckleyLeverett.Wells[0].Perforations[0].Rw = 0.3;
  BuckleyLeverett.Wells[0].Perforations[0].S = 0.25;
  BuckleyLeverett.Wells[0].Perforations[0].zbh = 0.0;

  /* Set information for well zero,the producer */
  BuckleyLeverett.Wells[1].Perforations[0].I = mx - 2;
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
  BuckleyLeverett.DeltaTS = 0.5;
  BuckleyLeverett.DSmax = 0.1;

  /* Set the start and end time */
  BuckleyLeverett.StartTime = 0.0;
  BuckleyLeverett.EndTime = 15.0;

  ierr = PetscPrintf(PETSC_COMM_WORLD,"x1 is: \n");
  ierr = VecView(BuckleyLeverett.x1,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;


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
  BlackOilReservoirSimulation FivePoint;
  PetscViewer viewer;

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
  ierr = VecSet(FivePoint.Muw, 0.5);CHKERRQ(ierr);CHKMEMQ;
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
  FivePoint.Wells[0].Perforations[0].Qw = 40.0;
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
  FivePoint.Wells[1].Perforations[0].BHPw = 3500;
  FivePoint.Wells[1].Perforations[0].BHPo = 3500;
  FivePoint.Wells[1].Perforations[0].Constraint = BHP_CONSTRAINT;
  FivePoint.Wells[1].Perforations[0].IsActive = PETSC_TRUE;
  FivePoint.Wells[1].Perforations[0].Orientation = PERF_ORIENTATION_X3X3;
  FivePoint.Wells[1].Perforations[0].Rw = 0.3;
  FivePoint.Wells[1].Perforations[0].S = 0.25;
  FivePoint.Wells[1].Perforations[0].zbh = 0.0;

  /* Not do adaptive time-step */
  FivePoint.AdaptiveTimeStep = PETSC_FALSE;

  /* Set the Times */
  FivePoint.DeltaTP = 1.0;
  FivePoint.DeltaTS = 0.1;
  FivePoint.DSmax = 0.1;

  /* Set the start and end time */
  FivePoint.StartTime = 0.0;
  FivePoint.EndTime = 40.0;

  /* Now we solve */
  ierr = DefiantIMPES2PhIterate(&FivePoint);CHKERRQ(ierr);CHKMEMQ;
  ierr = DefiantBlackOil2PhValiantWriteVecs(&FivePoint);CHKERRQ(ierr);CHKMEMQ;

#if 1
  ierr = PetscViewerCreate(PETSC_COMM_WORLD, &viewer);CHKERRQ(ierr);CHKMEMQ;
  ierr = PetscViewerSetType(viewer, PETSC_VIEWER_ASCII);CHKERRQ(ierr);CHKMEMQ;
  ierr = PetscViewerSetFormat(viewer, 	PETSC_VIEWER_ASCII_MATLAB);CHKERRQ(ierr);CHKMEMQ;
  ierr = PetscViewerFileSetName(viewer, "FivePoint.m");CHKERRQ(ierr);CHKMEMQ;
  ierr = VecView(FivePoint.Sw,viewer);CHKERRQ(ierr);CHKMEMQ;
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

#undef __FUNCT__
#define __FUNCT__ "DefiantIMPES2PhBenchmark"
extern PetscErrorCode DefiantIMPES2PhBenchmark()
{
  PetscErrorCode ierr;
  PetscInt mx, my, mz, xm, ym, zm, xs, ys, zs;
  BlackOilReservoirSimulation Benchmark;

  PetscFunctionBegin;
  /* Initialize the simulation */
  ierr = DefiantCreateSimulationVecs(&Benchmark);CHKERRQ(ierr);CHKMEMQ;
  ierr = DefiantCreate2PhIMPESSystem(&Benchmark);CHKERRQ(ierr);CHKMEMQ;

  /* Set FlowMask */
  ierr = VecSet(Benchmark.FlowMask, FLUID_FLOW);CHKERRQ(ierr);CHKMEMQ;

  /* Set the coordinates */
  ierr = DASetUniformCoordinates(Benchmark.SimDA, 0.0, 90.0, 0.0, 2.0, 0.0, 2.0);CHKERRQ(ierr);CHKMEMQ;
  ierr = DefiantGetDACoords(&Benchmark);CHKERRQ(ierr);CHKMEMQ;
  /* Get dimensions and extents of the local vectors */
  ierr = DAGetInfo(Benchmark.SimDA, 0, &mx, &my, &mz, 0, 0, 0, 0, 0, 0, 0);CHKERRQ(ierr);
  ierr = DAGetCorners(Benchmark.SimDA, &xs, &ys, &zs, &xm, &ym, &zm);CHKERRQ(ierr);
  /* set some geometry */
  ierr = VecSet(Benchmark.A1, 1.0);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSet(Benchmark.A2, 10.0);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSet(Benchmark.A3, 10.0);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSet(Benchmark.h1, 10.0);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSet(Benchmark.h2, 1.0);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSet(Benchmark.h3, 1.0);CHKERRQ(ierr);CHKMEMQ;
  /* set initial values to vectors */
  ierr = VecSet(Benchmark.Po, 4000.0);CHKERRQ(ierr);CHKMEMQ;
  /* set residual oil and connate water */
  Benchmark.Sor = 0.2;
  Benchmark.Swc = 0.2;
  /* Set saturations */
  ierr = VecSet(Benchmark.So, 1.0-0.2 );CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSet(Benchmark.Sw, 0.2 );CHKERRQ(ierr);CHKMEMQ;
  /* Set densities */
  Benchmark.Rhoos = 1.0;
  Benchmark.Rhows = 1.0;
  ierr = VecSet(Benchmark.Rhoo, Benchmark.Rhoos );CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSet(Benchmark.Rhow, Benchmark.Rhows );CHKERRQ(ierr);CHKMEMQ;
  /* Set viscosities */
  ierr = VecSet(Benchmark.Muo, 1);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSet(Benchmark.Muw, 0.3);CHKERRQ(ierr);CHKMEMQ;
  /* Set permeabilities */
  ierr = VecSet(Benchmark.K11, 300.0);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSet(Benchmark.K22, 300.0);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSet(Benchmark.K33, 300.0);CHKERRQ(ierr);CHKMEMQ;
  /* set the relative permeabilities */
  ierr = VecSet(Benchmark.Kro, 1.0);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSet(Benchmark.Krw, 1.0);CHKERRQ(ierr);CHKMEMQ;
  /* Set the table for PcowSw */
  Benchmark.PcowSw.NumberOfEntries = 2;
  ierr = PetscMalloc(Benchmark.PcowSw.NumberOfEntries*sizeof(PetscScalar),&Benchmark.PcowSw.X);CHKERRQ(ierr);CHKMEMQ;
  ierr = PetscMalloc(Benchmark.PcowSw.NumberOfEntries*sizeof(PetscScalar),&Benchmark.PcowSw.Y);CHKERRQ(ierr);CHKMEMQ;
  Benchmark.PcowSw.X[0] = 0;Benchmark.PcowSw.X[1] = 0;
  Benchmark.PcowSw.Y[0] = 0;Benchmark.PcowSw.Y[1] = 1;
  /* Set the graviational acceleration */
  Benchmark.GravAcc = 32.2;
  /* Set the rock properties */
  ierr = VecSet(Benchmark.Phi, 0.2);CHKERRQ(ierr);CHKMEMQ;
  Benchmark.RockCompressibility = 1.0e-15;
  Benchmark.RockCompRefPressure = 14.7;
  Benchmark.RockCompRefPorosity = 0.3;
  /* Set the fluid compressibilities */
  Benchmark.OilCompressibility = 1.0e-5;
  Benchmark.OilCompPressure = 14.7;
  Benchmark.WaterCompressibility = 1.0e-5;
  Benchmark.WaterCompPressure = 14.7;
  /* FIXME Set the volume factors */
  ierr = VecSet(Benchmark.Bo, 1.0);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSet(Benchmark.Bw, 1.0);CHKERRQ(ierr);CHKMEMQ;
  /* Set the well information */
  /* allocate wells */
  Benchmark.NumberOfWells = 2;
  ierr = PetscMalloc(Benchmark.NumberOfWells*sizeof(Well),&Benchmark.Wells);CHKERRQ(ierr);CHKMEMQ;
  /* allocate perforations for wells */
  Benchmark.Wells[0].NumberOfPerforations = 1;
  ierr = PetscMalloc(Benchmark.Wells[0].NumberOfPerforations*sizeof(Perforation),&Benchmark.Wells[0].Perforations);CHKERRQ(ierr);CHKMEMQ;
  Benchmark.Wells[1].NumberOfPerforations = 1;
  ierr = PetscMalloc(Benchmark.Wells[1].NumberOfPerforations*sizeof(Perforation),&Benchmark.Wells[1].Perforations);CHKERRQ(ierr);CHKMEMQ;
  /* Set information for well zero,the injector */
  Benchmark.Wells[0].Perforations[0].I = 1;
  Benchmark.Wells[0].Perforations[0].J = 1;
  Benchmark.Wells[0].Perforations[0].K = 1;
  Benchmark.Wells[0].Perforations[0].Qw = 1.0;
  Benchmark.Wells[0].Perforations[0].Constraint = FLOW_RATE_CONSTRAINT;
  Benchmark.Wells[0].Perforations[0].WellType = WATER_INJECTOR;
  Benchmark.Wells[0].Perforations[0].IsActive = PETSC_TRUE;
  Benchmark.Wells[0].Perforations[0].Orientation = PERF_ORIENTATION_X3X3;
  Benchmark.Wells[0].Perforations[0].Rw = 0.3;
  Benchmark.Wells[0].Perforations[0].S = 0.25;
  Benchmark.Wells[0].Perforations[0].zbh = 0.0;

  /* Set information for well zero,the producer */
  Benchmark.Wells[1].Perforations[0].I = mx - 2;
  Benchmark.Wells[1].Perforations[0].J = my - 2;
  Benchmark.Wells[1].Perforations[0].K = mz - 2;
  Benchmark.Wells[1].Perforations[0].BHPw = 3999;
  Benchmark.Wells[1].Perforations[0].BHPo = 3999;
  Benchmark.Wells[1].Perforations[0].Constraint = BHP_CONSTRAINT;
  Benchmark.Wells[1].Perforations[0].IsActive = PETSC_TRUE;
  Benchmark.Wells[1].Perforations[0].Orientation = PERF_ORIENTATION_X3X3;
  Benchmark.Wells[1].Perforations[0].Rw = 0.3;
  Benchmark.Wells[1].Perforations[0].S = 0.25;
  Benchmark.Wells[1].Perforations[0].zbh = 0.0;

  /* Not do adaptive time-step */
  Benchmark.AdaptiveTimeStep = PETSC_FALSE;

  /* Set the Times */
  Benchmark.DeltaTP = 0.1;
  Benchmark.DeltaTS = 0.1;
  Benchmark.DSmax = 0.1;

  /* Set the start and end time */
  Benchmark.StartTime = 0.0;
  Benchmark.EndTime = 0.5;

  /* Now we solve */
  ierr = DefiantIMPES2PhIterate(&Benchmark);CHKERRQ(ierr);CHKMEMQ;
  //ierr = DefiantBlackOil2PhValiantWriteVecs(&Benchmark);CHKERRQ(ierr);CHKMEMQ;

  //ierr = PetscPrintf(PETSC_COMM_WORLD,"Sw is: \n");
  //ierr = VecView(Benchmark.Sw,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
  //ierr = PetscPrintf(PETSC_COMM_WORLD,"x1 is: \n");
  //ierr = VecView(Benchmark.x1,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
  //ierr = PetscPrintf(PETSC_COMM_WORLD,"h1 is: \n");
  //ierr = VecView(Benchmark.h1,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;

  /* Now we clean up */
  ierr = DefiantDestroySimulationVecs(&Benchmark);CHKERRQ(ierr);CHKMEMQ;
  ierr = DefiantDestroy2PhIMPESSystem(&Benchmark);CHKERRQ(ierr);CHKMEMQ;


  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DefiantIMPES2PhTestValiant"
extern PetscErrorCode DefiantIMPES2PhTestValiant()
{
  PetscErrorCode ierr;
  PetscInt       MySeed1 = 1;
  PetscInt       MySeed2 = 2;
  PetscRandom    rctx;
  PetscTruth     BoolPerturb = PETSC_FALSE;
  PetscTruth     BoolRunPerturbed = PETSC_FALSE;
  PetscTruth     BoolRunNormal = PETSC_FALSE;
  PetscInt mx, my, mz, xm, ym, zm, xs, ys, zs;
  Vec TempVec, RandomVec;  /* temporary vector */
  PetscScalar Percentage1, Percentage2, TempScalar;
  BlackOilReservoirSimulation TestValiant;

  PetscFunctionBegin;
  /* Initialize the simulation */
  ierr = DefiantCreateSimulationVecs(&TestValiant);CHKERRQ(ierr);CHKMEMQ;
  ierr = DefiantCreate2PhIMPESSystem(&TestValiant);CHKERRQ(ierr);CHKMEMQ;

  /* Set FlowMask */
  ierr = VecSet(TestValiant.FlowMask, FLUID_FLOW);CHKERRQ(ierr);CHKMEMQ;

  /* Set the coordinates */
  ierr = DASetUniformCoordinates(TestValiant.SimDA, 0.0, 90.0, 0.0, 2.0, 0.0, 2.0);CHKERRQ(ierr);CHKMEMQ;
  ierr = DefiantGetDACoords(&TestValiant);CHKERRQ(ierr);CHKMEMQ;
  /* Get dimensions and extents of the local vectors */
  ierr = DAGetInfo(TestValiant.SimDA, 0, &mx, &my, &mz, 0, 0, 0, 0, 0, 0, 0);CHKERRQ(ierr);
  ierr = DAGetCorners(TestValiant.SimDA, &xs, &ys, &zs, &xm, &ym, &zm);CHKERRQ(ierr);
  /* set some geometry */
  ierr = VecSet(TestValiant.A1, 1.0);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSet(TestValiant.A2, 10.0);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSet(TestValiant.A3, 10.0);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSet(TestValiant.h1, 10.0);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSet(TestValiant.h2, 1.0);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSet(TestValiant.h3, 1.0);CHKERRQ(ierr);CHKMEMQ;
  /* set initial values to vectors */
  ierr = VecSet(TestValiant.Po, 10.0);CHKERRQ(ierr);CHKMEMQ;
  /* set residual oil and connate water */
  TestValiant.Sor = 0.2;
  TestValiant.Swc = 0.2;
  /* Set saturations */
  ierr = VecSet(TestValiant.So, 1.0-0.4 );CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSet(TestValiant.Sw, 0.4 );CHKERRQ(ierr);CHKMEMQ;
  /* Set densities */
  TestValiant.Rhoos = 1.0;
  TestValiant.Rhows = 1.0;
  ierr = VecSet(TestValiant.Rhoo, TestValiant.Rhoos );CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSet(TestValiant.Rhow, TestValiant.Rhows );CHKERRQ(ierr);CHKMEMQ;
  /* Set viscosities */
  ierr = VecSet(TestValiant.Muo, 1);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSet(TestValiant.Muw, 1);CHKERRQ(ierr);CHKMEMQ;
  /* Set permeabilities */
  ierr = VecSet(TestValiant.K11, 300.0);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSet(TestValiant.K22, 300.0);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSet(TestValiant.K33, 300.0);CHKERRQ(ierr);CHKMEMQ;
  /* set the relative permeabilities */
  ierr = VecSet(TestValiant.Kro, 1.0);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSet(TestValiant.Krw, 1.0);CHKERRQ(ierr);CHKMEMQ;
  /* Set the table for PcowSw */
  TestValiant.PcowSw.NumberOfEntries = 2;
  ierr = PetscMalloc(TestValiant.PcowSw.NumberOfEntries*sizeof(PetscScalar),&TestValiant.PcowSw.X);CHKERRQ(ierr);CHKMEMQ;
  ierr = PetscMalloc(TestValiant.PcowSw.NumberOfEntries*sizeof(PetscScalar),&TestValiant.PcowSw.Y);CHKERRQ(ierr);CHKMEMQ;
  TestValiant.PcowSw.X[0] = 0;TestValiant.PcowSw.X[1] = 0;
  TestValiant.PcowSw.Y[0] = 0;TestValiant.PcowSw.Y[1] = 1;
  /* Set the graviational acceleration */
  TestValiant.GravAcc = 32.2;
  /* Set the rock properties */
  ierr = VecSet(TestValiant.Phi, 0.2);CHKERRQ(ierr);CHKMEMQ;
  TestValiant.RockCompressibility = 0.0;
  TestValiant.RockCompRefPressure = 14.7;
  TestValiant.RockCompRefPorosity = 0.3;
  /* Set the fluid compressibilities */
  TestValiant.OilCompressibility = 1.0e-5;
  TestValiant.OilCompPressure = 14.7;
  TestValiant.WaterCompressibility = 1.0e-5;
  TestValiant.WaterCompPressure = 14.7;
  /* FIXME Set the volume factors */
  ierr = VecSet(TestValiant.Bo, 1.0);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSet(TestValiant.Bw, 1.0);CHKERRQ(ierr);CHKMEMQ;
  /* Set the well information */
  /* allocate wells */
  TestValiant.NumberOfWells = 3;
  ierr = PetscMalloc(TestValiant.NumberOfWells*sizeof(Well),&TestValiant.Wells);CHKERRQ(ierr);CHKMEMQ;
  /* allocate perforations for wells */
  TestValiant.Wells[0].NumberOfPerforations = 1;
  ierr = PetscMalloc(TestValiant.Wells[0].NumberOfPerforations*sizeof(Perforation),&TestValiant.Wells[0].Perforations);CHKERRQ(ierr);CHKMEMQ;
  TestValiant.Wells[1].NumberOfPerforations = 1;
  ierr = PetscMalloc(TestValiant.Wells[1].NumberOfPerforations*sizeof(Perforation),&TestValiant.Wells[1].Perforations);CHKERRQ(ierr);CHKMEMQ;
  TestValiant.Wells[2].NumberOfPerforations = 1;
  ierr = PetscMalloc(TestValiant.Wells[2].NumberOfPerforations*sizeof(Perforation),&TestValiant.Wells[2].Perforations);CHKERRQ(ierr);CHKMEMQ;
  /* Set information for well zero,the injector */
  TestValiant.Wells[0].Perforations[0].I = 1;
  TestValiant.Wells[0].Perforations[0].J = 1;
  TestValiant.Wells[0].Perforations[0].K = 1;
  TestValiant.Wells[0].Perforations[0].Qw = 1.0;
  TestValiant.Wells[0].Perforations[0].Constraint = FLOW_RATE_CONSTRAINT;
  TestValiant.Wells[0].Perforations[0].WellType = WATER_INJECTOR;
  TestValiant.Wells[0].Perforations[0].IsActive = PETSC_TRUE;
  TestValiant.Wells[0].Perforations[0].Orientation = PERF_ORIENTATION_X3X3;
  TestValiant.Wells[0].Perforations[0].Rw = 0.3;
  TestValiant.Wells[0].Perforations[0].S = 0.25;
  TestValiant.Wells[0].Perforations[0].zbh = 0.0;

  /* Set information for well one,the producer */
  TestValiant.Wells[1].Perforations[0].I = mx - 2;
  TestValiant.Wells[1].Perforations[0].J = my - 2;
  TestValiant.Wells[1].Perforations[0].K = mz - 2;
  TestValiant.Wells[1].Perforations[0].BHPw = 5;
  TestValiant.Wells[1].Perforations[0].BHPo = 5;
  TestValiant.Wells[1].Perforations[0].Constraint = BHP_CONSTRAINT;
  TestValiant.Wells[1].Perforations[0].IsActive = PETSC_TRUE;
  TestValiant.Wells[1].Perforations[0].Orientation = PERF_ORIENTATION_X3X3;
  TestValiant.Wells[1].Perforations[0].Rw = 0.3;
  TestValiant.Wells[1].Perforations[0].S = 0.25;
  TestValiant.Wells[1].Perforations[0].zbh = 0.0;


  /* Set information for well two,the monitoring well */
  TestValiant.Wells[2].Perforations[0].I = (mx - 2)/2;
  TestValiant.Wells[2].Perforations[0].J = (my - 2)/2;
  TestValiant.Wells[2].Perforations[0].K = (mz - 2)/2;
  TestValiant.Wells[2].Perforations[0].Constraint = MONITORING;
  TestValiant.Wells[2].Perforations[0].IsActive = PETSC_TRUE;
  TestValiant.Wells[2].Perforations[0].Orientation = PERF_ORIENTATION_X3X3;

  /* Not do adaptive time-step */
  TestValiant.AdaptiveTimeStep = PETSC_FALSE;

  /* Check if we perturb */
  ierr = PetscOptionsGetTruth(PETSC_NULL, "-perturb", &BoolPerturb, PETSC_NULL);CHKERRQ(ierr);
  /* Check if we run the perturbed case to completion */
  ierr = PetscOptionsGetTruth(PETSC_NULL, "-runperturbed", &BoolRunPerturbed, PETSC_NULL);CHKERRQ(ierr);
  /* Check if we run normally */
  ierr = PetscOptionsGetTruth(PETSC_NULL, "-runnormal", &BoolRunNormal, PETSC_NULL);CHKERRQ(ierr);

  if(BoolPerturb)
  {
    ierr = PetscOptionsGetInt(PETSC_NULL,"-seed_phi",&MySeed1,PETSC_NULL);CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(PETSC_NULL,"-seed_k11",&MySeed2,PETSC_NULL);CHKERRQ(ierr);

    ierr = PetscOptionsGetReal(PETSC_NULL,"-percentage_phi",&Percentage1,PETSC_NULL);CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(PETSC_NULL,"-percentage_k11",&Percentage2,PETSC_NULL);CHKERRQ(ierr);

    /* create the random vector */
    ierr = VecDuplicate(TestValiant.Phi,&RandomVec);CHKERRQ(ierr);

    /* perturb the porosity */
    ierr = VecDuplicate(TestValiant.Phi,&TempVec);CHKERRQ(ierr);
    ierr = VecCopy(TestValiant.Phi, TempVec);CHKERRQ(ierr);
    ierr = VecScale(TempVec, Percentage1);CHKERRQ(ierr);
    ierr = PetscRandomCreate(PETSC_COMM_WORLD,&rctx);CHKERRQ(ierr);
    ierr = PetscRandomSetSeed(rctx,MySeed1);CHKERRQ(ierr);
    ierr = PetscRandomSetFromOptions(rctx);CHKERRQ(ierr);
    ierr = VecSetRandom(RandomVec,rctx);CHKERRQ(ierr);
    ierr = PetscRandomDestroy(rctx);CHKERRQ(ierr);
    ierr = VecPointwiseMult(TempVec,TempVec,RandomVec);CHKERRQ(ierr);
    TempScalar = 1.0;
    ierr = VecAYPX(TestValiant.Phi,TempScalar,TempVec);CHKERRQ(ierr);
    ierr = VecAssemblyBegin(TestValiant.Phi);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(TestValiant.Phi);CHKERRQ(ierr);

    /* perturb the permeability along the x11 direction */
    ierr = VecLog(TestValiant.K11);CHKERRQ(ierr);
    ierr = VecCopy(TestValiant.K11, TempVec);CHKERRQ(ierr);
    ierr = VecScale(TempVec, Percentage2);CHKERRQ(ierr);
    ierr = PetscRandomCreate(PETSC_COMM_WORLD,&rctx);CHKERRQ(ierr);
    ierr = PetscRandomSetSeed(rctx,MySeed2);CHKERRQ(ierr);
    ierr = PetscRandomSetFromOptions(rctx);CHKERRQ(ierr);
    ierr = VecSetRandom(RandomVec,rctx);CHKERRQ(ierr);
    ierr = PetscRandomDestroy(rctx);CHKERRQ(ierr);
    ierr = VecPointwiseMult(TempVec,TempVec,RandomVec);CHKERRQ(ierr);
    TempScalar = 1.0;
    ierr = VecAYPX(TestValiant.K11,TempScalar,TempVec);CHKERRQ(ierr);
    ierr = VecExp(TestValiant.K11);CHKERRQ(ierr);
    ierr = VecAssemblyBegin(TestValiant.K11);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(TestValiant.K11);CHKERRQ(ierr);

    ierr = VecDestroy(TempVec);CHKERRQ(ierr);
    ierr = VecDestroy(RandomVec);CHKERRQ(ierr);

    ierr = DefiantBlackOil2PhValiantWriteVecs(&TestValiant);CHKERRQ(ierr);CHKMEMQ;
  }

  if (BoolRunPerturbed)
  {

    //ierr = PetscPrintf(PETSC_COMM_WORLD,"Post Perturbed Phi is: with seed: %d \n", MySeed1);
    //ierr = VecView(TestValiant.Phi,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
    //ierr = PetscPrintf(PETSC_COMM_WORLD,"Post Perturbed K11 is: with seed: %d \n", MySeed2);
    //ierr = VecView(TestValiant.K11,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;

    /* Set the Times */
    TestValiant.DeltaTP = 0.5;
    TestValiant.DeltaTS = 0.1;
    TestValiant.DSmax = 0.1;

    /* Set the start and end time */
    TestValiant.StartTime = 0.0;
    TestValiant.EndTime = 20.0;

    /* Now we solve */
    ierr = DefiantIMPES2PhIterate(&TestValiant);CHKERRQ(ierr);CHKMEMQ;
    /* Now we write the fields */
    ierr = DefiantBlackOil2PhValiantWriteVecs(&TestValiant);CHKERRQ(ierr);CHKMEMQ;

    ierr = PetscPrintf(PETSC_COMM_WORLD,"After perturbed run Sw is: \n");
    ierr = VecView(TestValiant.Sw,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"After perturbed run x1 is: \n");
    ierr = VecView(TestValiant.x1,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
  } else  if  (BoolRunNormal) {
    /* Set the Times */
    TestValiant.DeltaTP = 0.5;
    TestValiant.DeltaTS = 0.1;
    TestValiant.DSmax = 0.1;

    /* Set the start and end time */
    TestValiant.StartTime = 0.0;
    TestValiant.EndTime = 1.0;

    /* we first load the modified fields */
    ierr = DefiantBlackOil2PhValiantLoadVecs(&TestValiant);CHKERRQ(ierr);CHKMEMQ;

    //ierr = PetscPrintf(PETSC_COMM_WORLD,"Post Perturbed Phi is: with seed: %d \n", MySeed1);
    //ierr = VecView(TestValiant.Phi,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
    //ierr = PetscPrintf(PETSC_COMM_WORLD,"Post Perturbed K11 is: with seed: %d \n", MySeed2);
    //ierr = VecView(TestValiant.K11,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;

    /* Now we solve */
    ierr = DefiantIMPES2PhIterate(&TestValiant);CHKERRQ(ierr);CHKMEMQ;
    /* Now we write the fields */
    ierr = DefiantBlackOil2PhValiantWriteVecs(&TestValiant);CHKERRQ(ierr);CHKMEMQ;

    //ierr = PetscPrintf(PETSC_COMM_WORLD,"Sw is: \n");
    //ierr = VecView(TestValiant.Sw,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
    //ierr = PetscPrintf(PETSC_COMM_WORLD,"x1 is: \n");
    //ierr = VecView(TestValiant.x1,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
  }

  /* Now we clean up */
  ierr = DefiantDestroySimulationVecs(&TestValiant);CHKERRQ(ierr);CHKMEMQ;
  ierr = DefiantDestroy2PhIMPESSystem(&TestValiant);CHKERRQ(ierr);CHKMEMQ;


  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "DefiantIMPES2PhAziz"
extern PetscErrorCode DefiantIMPES2PhAziz()
{
  PetscErrorCode ierr;
  PetscInt mx, my, mz, xm, ym, zm, xs, ys, zs;
  PetscScalar ***Localx1;
  PetscScalar TempScalar;

  BlackOilReservoirSimulation BuckleyLeverett;

  PetscFunctionBegin;
  /* Initialize the simulation */
  ierr = DefiantCreateSimulationVecs(&BuckleyLeverett);CHKERRQ(ierr);CHKMEMQ;
  ierr = DefiantCreate2PhIMPESSystem(&BuckleyLeverett);CHKERRQ(ierr);CHKMEMQ;

  /* Set FlowMask */
  ierr = VecSet(BuckleyLeverett.FlowMask, FLUID_FLOW);CHKERRQ(ierr);CHKMEMQ;

  /* Set the coordinates */
  ierr = DASetUniformCoordinates(BuckleyLeverett.SimDA, -100.0, 1100.0, 0.0, 2.0, 0.0, 2.0);CHKERRQ(ierr);CHKMEMQ;
  ierr = DAGetInfo(BuckleyLeverett.SimDA, 0, &mx, &my, &mz, 0, 0, 0, 0, 0, 0, 0);CHKERRQ(ierr);
  ierr = DAGetCorners(BuckleyLeverett.SimDA, &xs, &ys, &zs, &xm, &ym, &zm);CHKERRQ(ierr);
  ierr = DefiantGetDACoords(&BuckleyLeverett);CHKERRQ(ierr);CHKMEMQ;
  /* set some geometry */
  ierr = VecSet(BuckleyLeverett.A1, 10000.0);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSet(BuckleyLeverett.A2, 10.0);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSet(BuckleyLeverett.A3, 10.0);CHKERRQ(ierr);CHKMEMQ;

  ierr = DAVecGetArray(BuckleyLeverett.SimDA,BuckleyLeverett.x1, &Localx1);CHKERRQ(ierr);CHKMEMQ;
  TempScalar = Localx1[1][1][1]-Localx1[1][1][0];
  ierr = PetscPrintf(PETSC_COMM_WORLD,"We Set h1 to: %g\n", TempScalar);

  ierr = VecSet(BuckleyLeverett.h1, TempScalar);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSet(BuckleyLeverett.h2, 100.0);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSet(BuckleyLeverett.h3, 100.0);CHKERRQ(ierr);CHKMEMQ;
  /* set initial values to vectors */
  ierr = VecSet(BuckleyLeverett.Po, 4000.0);CHKERRQ(ierr);CHKMEMQ;
  /* set initial values to vectors */
  ierr = VecSet(BuckleyLeverett.Pw, 4000.0);CHKERRQ(ierr);CHKMEMQ;
  /* set residual oil and connate water */
  BuckleyLeverett.Sor = 0.2;
  BuckleyLeverett.Swc = 0.16;
  /* Set saturations */
  ierr = VecSet(BuckleyLeverett.So, 1.0-0.16 );CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSet(BuckleyLeverett.Sw, 0.16 );CHKERRQ(ierr);CHKMEMQ;
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
  ierr = VecSet(BuckleyLeverett.Krw, 0.0);CHKERRQ(ierr);CHKMEMQ;
  /* Set the table for PcowSw */
  BuckleyLeverett.PcowSw.NumberOfEntries = 2;
  ierr = PetscMalloc(BuckleyLeverett.PcowSw.NumberOfEntries*sizeof(PetscScalar),&BuckleyLeverett.PcowSw.X);CHKERRQ(ierr);CHKMEMQ;
  ierr = PetscMalloc(BuckleyLeverett.PcowSw.NumberOfEntries*sizeof(PetscScalar),&BuckleyLeverett.PcowSw.Y);CHKERRQ(ierr);CHKMEMQ;
  BuckleyLeverett.PcowSw.X[0] = 0.1;BuckleyLeverett.PcowSw.X[1] = 0.16;
  BuckleyLeverett.PcowSw.Y[0] = 0.0;BuckleyLeverett.PcowSw.Y[1] = 0.8;
  /* Set the graviational acceleration */
  BuckleyLeverett.GravAcc = 0.0;
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
  BuckleyLeverett.Wells[0].Perforations[0].Qw = 426.5;
  BuckleyLeverett.Wells[0].Perforations[0].Constraint = FLOW_RATE_CONSTRAINT;
  BuckleyLeverett.Wells[0].Perforations[0].WellType = WATER_INJECTOR;
  BuckleyLeverett.Wells[0].Perforations[0].IsActive = PETSC_TRUE;
  BuckleyLeverett.Wells[0].Perforations[0].Orientation = PERF_ORIENTATION_X3X3;
  BuckleyLeverett.Wells[0].Perforations[0].Rw = 0.3;
  BuckleyLeverett.Wells[0].Perforations[0].S = 0.0;
  BuckleyLeverett.Wells[0].Perforations[0].zbh = 0.0;

  /* Set information for well zero,the producer */
  BuckleyLeverett.Wells[1].Perforations[0].I = mx-2;
  BuckleyLeverett.Wells[1].Perforations[0].J = 1;
  BuckleyLeverett.Wells[1].Perforations[0].K = 1;
  BuckleyLeverett.Wells[1].Perforations[0].BHPw = 3300;
  BuckleyLeverett.Wells[1].Perforations[0].BHPo = 3300;
  BuckleyLeverett.Wells[1].Perforations[0].Constraint = BHP_CONSTRAINT;
  BuckleyLeverett.Wells[1].Perforations[0].WellType = OIL_PRODUCER;  
  BuckleyLeverett.Wells[1].Perforations[0].IsActive = PETSC_TRUE;
  BuckleyLeverett.Wells[1].Perforations[0].Orientation = PERF_ORIENTATION_X3X3;
  BuckleyLeverett.Wells[1].Perforations[0].Rw = 0.3;
  BuckleyLeverett.Wells[1].Perforations[0].S = 0.0;
  BuckleyLeverett.Wells[1].Perforations[0].zbh = 0.0;

  /* Not do adaptive time-step */
  BuckleyLeverett.AdaptiveTimeStep = PETSC_FALSE;

  /* Set the Times */
  BuckleyLeverett.DeltaTP = 10.0;
  BuckleyLeverett.DeltaTS = 1.0;
  BuckleyLeverett.DSmax = 0.1;

  /* Set the start and end time */
  BuckleyLeverett.StartTime = 0.0;
  BuckleyLeverett.EndTime = 1500.0;

  /* Now we solve */
  ierr = DefiantIMPES2PhIterate(&BuckleyLeverett);CHKERRQ(ierr);CHKMEMQ;
  ierr = DefiantBlackOil2PhValiantWriteVecs(&BuckleyLeverett);CHKERRQ(ierr);CHKMEMQ;

  ierr = PetscPrintf(PETSC_COMM_WORLD,"Sw is: \n");
  ierr = VecView(BuckleyLeverett.Sw,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Pw is: \n");
  ierr = VecView(BuckleyLeverett.Pw,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
  ierr = PetscPrintf(PETSC_COMM_WORLD,"x1 is: \n");
  ierr = VecView(BuckleyLeverett.x1,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;

  /* Now we clean up */
  ierr = DefiantDestroySimulationVecs(&BuckleyLeverett);CHKERRQ(ierr);CHKMEMQ;
  ierr = DefiantDestroy2PhIMPESSystem(&BuckleyLeverett);CHKERRQ(ierr);CHKMEMQ;


  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "DefiantIMPES2PhPUNQ3"
extern PetscErrorCode DefiantIMPES2PhPUNQ3()
{
  PetscErrorCode ierr;
  PetscInt       MySeed1 = 1;
  PetscInt       MySeed2 = 2;
  PetscRandom    rctx;
  PetscTruth     BoolPerturb = PETSC_FALSE;
  PetscTruth     BoolRunPerturbed = PETSC_FALSE;
  PetscTruth     BoolRunNormal = PETSC_FALSE;
  PetscInt mx, my, mz, xm, ym, zm, xs, ys, zs;
  Vec TempVec, RandomVec;  /* temporary vector */
  PetscScalar Percentage1, Percentage2, TempScalar;
  BlackOilReservoirSimulation PUNQ3;

  PetscFunctionBegin;
  /* Initialize the simulation */
  ierr = DefiantCreateSimulationVecs(&PUNQ3);CHKERRQ(ierr);CHKMEMQ;
  ierr = DefiantCreate2PhIMPESSystem(&PUNQ3);CHKERRQ(ierr);CHKMEMQ;

  /* Set FlowMask */
  ierr = VecSet(PUNQ3.FlowMask, FLUID_FLOW);CHKERRQ(ierr);CHKMEMQ;

  /* Set the coordinates */
  ierr = DASetUniformCoordinates(PUNQ3.SimDA, 0.0, 90.0, 0.0, 2.0, 0.0, 2.0);CHKERRQ(ierr);CHKMEMQ;
  ierr = DefiantGetDACoords(&PUNQ3);CHKERRQ(ierr);CHKMEMQ;
  /* Get dimensions and extents of the local vectors */
  ierr = DAGetInfo(PUNQ3.SimDA, 0, &mx, &my, &mz, 0, 0, 0, 0, 0, 0, 0);CHKERRQ(ierr);
  ierr = DAGetCorners(PUNQ3.SimDA, &xs, &ys, &zs, &xm, &ym, &zm);CHKERRQ(ierr);
  /* set some geometry */
  ierr = VecSet(PUNQ3.A1, 1.0);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSet(PUNQ3.A2, 10.0);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSet(PUNQ3.A3, 10.0);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSet(PUNQ3.h1, 10.0);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSet(PUNQ3.h2, 1.0);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSet(PUNQ3.h3, 1.0);CHKERRQ(ierr);CHKMEMQ;
  /* set initial values to vectors */
  ierr = VecSet(PUNQ3.Po, 10.0);CHKERRQ(ierr);CHKMEMQ;
  /* set residual oil and connate water */
  PUNQ3.Sor = 0.2;
  PUNQ3.Swc = 0.2;
  /* Set saturations */
  ierr = VecSet(PUNQ3.So, 1.0-0.2 );CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSet(PUNQ3.Sw, 0.2 );CHKERRQ(ierr);CHKMEMQ;
  /* Set densities */
  PUNQ3.Rhoos = 1.0;
  PUNQ3.Rhows = 1.0;
  ierr = VecSet(PUNQ3.Rhoo, PUNQ3.Rhoos );CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSet(PUNQ3.Rhow, PUNQ3.Rhows );CHKERRQ(ierr);CHKMEMQ;
  /* Set viscosities */
  ierr = VecSet(PUNQ3.Muo, 1);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSet(PUNQ3.Muw, 1);CHKERRQ(ierr);CHKMEMQ;
  /* Set permeabilities */
  ierr = VecSet(PUNQ3.K11, 300.0);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSet(PUNQ3.K22, 300.0);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSet(PUNQ3.K33, 300.0);CHKERRQ(ierr);CHKMEMQ;
  /* set the relative permeabilities */
  ierr = VecSet(PUNQ3.Kro, 1.0);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSet(PUNQ3.Krw, 1.0);CHKERRQ(ierr);CHKMEMQ;
  /* Set the table for PcowSw */
  PUNQ3.PcowSw.NumberOfEntries = 2;
  ierr = PetscMalloc(PUNQ3.PcowSw.NumberOfEntries*sizeof(PetscScalar),&PUNQ3.PcowSw.X);CHKERRQ(ierr);CHKMEMQ;
  ierr = PetscMalloc(PUNQ3.PcowSw.NumberOfEntries*sizeof(PetscScalar),&PUNQ3.PcowSw.Y);CHKERRQ(ierr);CHKMEMQ;
  PUNQ3.PcowSw.X[0] = 0;PUNQ3.PcowSw.X[1] = 0;
  PUNQ3.PcowSw.Y[0] = 0;PUNQ3.PcowSw.Y[1] = 1;
  /* Set the graviational acceleration */
  PUNQ3.GravAcc = 32.2;
  /* Set the rock properties */
  ierr = VecSet(PUNQ3.Phi, 0.2);CHKERRQ(ierr);CHKMEMQ;
  PUNQ3.RockCompressibility = 0.0;
  PUNQ3.RockCompRefPressure = 14.7;
  PUNQ3.RockCompRefPorosity = 0.3;
  /* Set the fluid compressibilities */
  PUNQ3.OilCompressibility = 1.0e-5;
  PUNQ3.OilCompPressure = 14.7;
  PUNQ3.WaterCompressibility = 1.0e-5;
  PUNQ3.WaterCompPressure = 14.7;
  /* FIXME Set the volume factors */
  ierr = VecSet(PUNQ3.Bo, 1.0);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSet(PUNQ3.Bw, 1.0);CHKERRQ(ierr);CHKMEMQ;
  /* Set the well information */
  /* allocate wells */
  PUNQ3.NumberOfWells = 3;
  ierr = PetscMalloc(PUNQ3.NumberOfWells*sizeof(Well),&PUNQ3.Wells);CHKERRQ(ierr);CHKMEMQ;
  /* allocate perforations for wells */
  PUNQ3.Wells[0].NumberOfPerforations = 1;
  ierr = PetscMalloc(PUNQ3.Wells[0].NumberOfPerforations*sizeof(Perforation),&PUNQ3.Wells[0].Perforations);CHKERRQ(ierr);CHKMEMQ;
  PUNQ3.Wells[1].NumberOfPerforations = 1;
  ierr = PetscMalloc(PUNQ3.Wells[1].NumberOfPerforations*sizeof(Perforation),&PUNQ3.Wells[1].Perforations);CHKERRQ(ierr);CHKMEMQ;
  PUNQ3.Wells[2].NumberOfPerforations = 1;
  ierr = PetscMalloc(PUNQ3.Wells[2].NumberOfPerforations*sizeof(Perforation),&PUNQ3.Wells[2].Perforations);CHKERRQ(ierr);CHKMEMQ;
  /* Set information for well zero,the injector */
  PUNQ3.Wells[0].Perforations[0].I = 1;
  PUNQ3.Wells[0].Perforations[0].J = 1;
  PUNQ3.Wells[0].Perforations[0].K = 1;
  PUNQ3.Wells[0].Perforations[0].Qw = 1.0;
  PUNQ3.Wells[0].Perforations[0].Constraint = FLOW_RATE_CONSTRAINT;
  PUNQ3.Wells[0].Perforations[0].WellType = WATER_INJECTOR;
  PUNQ3.Wells[0].Perforations[0].IsActive = PETSC_TRUE;
  PUNQ3.Wells[0].Perforations[0].Orientation = PERF_ORIENTATION_X3X3;
  PUNQ3.Wells[0].Perforations[0].Rw = 0.3;
  PUNQ3.Wells[0].Perforations[0].S = 0.25;
  PUNQ3.Wells[0].Perforations[0].zbh = 0.0;

  /* Set information for well one,the producer */
  PUNQ3.Wells[1].Perforations[0].I = mx - 2;
  PUNQ3.Wells[1].Perforations[0].J = my - 2;
  PUNQ3.Wells[1].Perforations[0].K = mz - 2;
  PUNQ3.Wells[1].Perforations[0].BHPw = 5;
  PUNQ3.Wells[1].Perforations[0].BHPo = 5;
  PUNQ3.Wells[1].Perforations[0].Constraint = BHP_CONSTRAINT;
  PUNQ3.Wells[1].Perforations[0].IsActive = PETSC_TRUE;
  PUNQ3.Wells[1].Perforations[0].Orientation = PERF_ORIENTATION_X3X3;
  PUNQ3.Wells[1].Perforations[0].Rw = 0.3;
  PUNQ3.Wells[1].Perforations[0].S = 0.25;
  PUNQ3.Wells[1].Perforations[0].zbh = 0.0;


  /* Set information for well two,the monitoring well */
  PUNQ3.Wells[2].Perforations[0].I = (mx - 2)/2;
  PUNQ3.Wells[2].Perforations[0].J = (my - 2)/2;
  PUNQ3.Wells[2].Perforations[0].K = (mz - 2)/2;
  PUNQ3.Wells[2].Perforations[0].Constraint = MONITORING;
  PUNQ3.Wells[2].Perforations[0].IsActive = PETSC_TRUE;
  PUNQ3.Wells[2].Perforations[0].Orientation = PERF_ORIENTATION_X3X3;

  /* Not do adaptive time-step */
  PUNQ3.AdaptiveTimeStep = PETSC_FALSE;

  /* Check if we perturb */
  ierr = PetscOptionsGetTruth(PETSC_NULL, "-perturb", &BoolPerturb, PETSC_NULL);CHKERRQ(ierr);
  /* Check if we run the perturbed case to completion */
  ierr = PetscOptionsGetTruth(PETSC_NULL, "-runperturbed", &BoolRunPerturbed, PETSC_NULL);CHKERRQ(ierr);
  /* Check if we run normally */
  ierr = PetscOptionsGetTruth(PETSC_NULL, "-runnormal", &BoolRunNormal, PETSC_NULL);CHKERRQ(ierr);

  if(BoolPerturb)
  {
    ierr = PetscOptionsGetInt(PETSC_NULL,"-seed_phi",&MySeed1,PETSC_NULL);CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(PETSC_NULL,"-seed_k11",&MySeed2,PETSC_NULL);CHKERRQ(ierr);

    ierr = PetscOptionsGetReal(PETSC_NULL,"-percentage_phi",&Percentage1,PETSC_NULL);CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(PETSC_NULL,"-percentage_k11",&Percentage2,PETSC_NULL);CHKERRQ(ierr);

    /* create the random vector */
    ierr = VecDuplicate(PUNQ3.Phi,&RandomVec);CHKERRQ(ierr);

    /* perturb the porosity */
    ierr = VecDuplicate(PUNQ3.Phi,&TempVec);CHKERRQ(ierr);
    ierr = VecCopy(PUNQ3.Phi, TempVec);CHKERRQ(ierr);
    ierr = VecScale(TempVec, Percentage1);CHKERRQ(ierr);
    ierr = PetscRandomCreate(PETSC_COMM_WORLD,&rctx);CHKERRQ(ierr);
    ierr = PetscRandomSetSeed(rctx,MySeed1);CHKERRQ(ierr);
    ierr = PetscRandomSetFromOptions(rctx);CHKERRQ(ierr);
    ierr = VecSetRandom(RandomVec,rctx);CHKERRQ(ierr);
    ierr = PetscRandomDestroy(rctx);CHKERRQ(ierr);
    ierr = VecPointwiseMult(TempVec,TempVec,RandomVec);CHKERRQ(ierr);
    TempScalar = 1.0;
    ierr = VecAYPX(PUNQ3.Phi,TempScalar,TempVec);CHKERRQ(ierr);
    ierr = VecAssemblyBegin(PUNQ3.Phi);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(PUNQ3.Phi);CHKERRQ(ierr);

    /* perturb the permeability along the x11 direction */
    ierr = VecLog(PUNQ3.K11);CHKERRQ(ierr);
    ierr = VecCopy(PUNQ3.K11, TempVec);CHKERRQ(ierr);
    ierr = VecScale(TempVec, Percentage2);CHKERRQ(ierr);
    ierr = PetscRandomCreate(PETSC_COMM_WORLD,&rctx);CHKERRQ(ierr);
    ierr = PetscRandomSetSeed(rctx,MySeed2);CHKERRQ(ierr);
    ierr = PetscRandomSetFromOptions(rctx);CHKERRQ(ierr);
    ierr = VecSetRandom(RandomVec,rctx);CHKERRQ(ierr);
    ierr = PetscRandomDestroy(rctx);CHKERRQ(ierr);
    ierr = VecPointwiseMult(TempVec,TempVec,RandomVec);CHKERRQ(ierr);
    TempScalar = 1.0;
    ierr = VecAYPX(PUNQ3.K11,TempScalar,TempVec);CHKERRQ(ierr);
    ierr = VecExp(PUNQ3.K11);CHKERRQ(ierr);
    ierr = VecAssemblyBegin(PUNQ3.K11);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(PUNQ3.K11);CHKERRQ(ierr);

    ierr = VecDestroy(TempVec);CHKERRQ(ierr);
    ierr = VecDestroy(RandomVec);CHKERRQ(ierr);

    ierr = DefiantBlackOil2PhValiantWriteVecs(&PUNQ3);CHKERRQ(ierr);CHKMEMQ;
  }

  if (BoolRunPerturbed)
  {

    ierr = PetscPrintf(PETSC_COMM_WORLD,"Post Perturbed Phi is: with seed: %d \n", MySeed1);
    ierr = VecView(PUNQ3.Phi,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Post Perturbed K11 is: with seed: %d \n", MySeed2);
    ierr = VecView(PUNQ3.K11,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;

    /* Set the Times */
    PUNQ3.DeltaTP = 0.5;
    PUNQ3.DeltaTS = 0.1;
    PUNQ3.DSmax = 0.1;

    /* Set the start and end time */
    PUNQ3.StartTime = 0.0;
    PUNQ3.EndTime = 20.0;

    /* Now we solve */
    ierr = DefiantIMPES2PhIterate(&PUNQ3);CHKERRQ(ierr);CHKMEMQ;
    /* Now we write the fields */
    ierr = DefiantBlackOil2PhValiantWriteVecs(&PUNQ3);CHKERRQ(ierr);CHKMEMQ;

    ierr = PetscPrintf(PETSC_COMM_WORLD,"After perturbed run Sw is: \n");
    ierr = VecView(PUNQ3.Sw,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"After perturbed run x1 is: \n");
    ierr = VecView(PUNQ3.x1,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
  } else  if  (BoolRunNormal) {
    /* Set the Times */
    PUNQ3.DeltaTP = 0.5;
    PUNQ3.DeltaTS = 0.1;
    PUNQ3.DSmax = 0.1;

    /* Set the start and end time */
    PUNQ3.StartTime = 0.0;
    PUNQ3.EndTime = 1.0;

    /* we first load the modified fields */
    ierr = DefiantBlackOil2PhValiantLoadVecs(&PUNQ3);CHKERRQ(ierr);CHKMEMQ;

    ierr = PetscPrintf(PETSC_COMM_WORLD,"Post Perturbed Phi is: with seed: %d \n", MySeed1);
    ierr = VecView(PUNQ3.Phi,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Post Perturbed K11 is: with seed: %d \n", MySeed2);
    ierr = VecView(PUNQ3.K11,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;

    /* Now we solve */
    ierr = DefiantIMPES2PhIterate(&PUNQ3);CHKERRQ(ierr);CHKMEMQ;
    /* Now we write the fields */
    ierr = DefiantBlackOil2PhValiantWriteVecs(&PUNQ3);CHKERRQ(ierr);CHKMEMQ;

    ierr = PetscPrintf(PETSC_COMM_WORLD,"Sw is: \n");
    ierr = VecView(PUNQ3.Sw,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"x1 is: \n");
    ierr = VecView(PUNQ3.x1,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
  }

  /* Now we clean up */
  ierr = DefiantDestroySimulationVecs(&PUNQ3);CHKERRQ(ierr);CHKMEMQ;
  ierr = DefiantDestroy2PhIMPESSystem(&PUNQ3);CHKERRQ(ierr);CHKMEMQ;


  PetscFunctionReturn(0);
}


//#undef __FUNCT__
//#define __FUNCT__ "DefiantIMPES2PhVictor"
//extern PetscErrorCode DefiantIMPES2PhVictor()
//{
  //PetscErrorCode ierr;
  //PetscInt       MySeed1 = 1;
  //PetscInt       MySeed2 = 2;
  //PetscInt       MySeed3 = 2;
  //PetscInt       MySeed4 = 2;
  //PetscInt       MySeed5 = 2;
  //PetscRandom    rctx;
  //PetscTruth     BoolPerturb = PETSC_FALSE;
  //PetscTruth     BoolRunPerturbed = PETSC_FALSE;
  //PetscTruth     BoolRunNormal = PETSC_FALSE;
  //PetscInt mx, my, mz, xm, ym, zm, xs, ys, zs;
  //Vec TempVec, RandomVec;  /* temporary vector */
  //PetscScalar Percentage1, Percentage2, Percentage3;
  //PetscScalar Percentage4, Percentage5, TempScalar;
  //PetscRandom    rctx;
  //PetscScalar  value;
  //PetscScalar ***LocalFlowMask;
  //BlackOilReservoirSimulation PUNQ3;

  //PetscFunctionBegin;
  ///* Initialize the simulation */
  //ierr = DefiantCreateSimulationVecs(&PUNQ3);CHKERRQ(ierr);CHKMEMQ;
  //ierr = DefiantCreate2PhIMPESSystem(&PUNQ3);CHKERRQ(ierr);CHKMEMQ;

  ///* Set FlowMask */
  //ierr = VecSet(PUNQ3.FlowMask, FLUID_FLOW);CHKERRQ(ierr);CHKMEMQ;

  ///* Set the coordinates */
  //ierr = DASetUniformCoordinates(PUNQ3.SimDA, 0.0, 90.0, 0.0, 2.0, 0.0, 2.0);CHKERRQ(ierr);CHKMEMQ;
  //ierr = DefiantGetDACoords(&PUNQ3);CHKERRQ(ierr);CHKMEMQ;
  ///* Get dimensions and extents of the local vectors */
  //ierr = DAGetInfo(PUNQ3.SimDA, 0, &mx, &my, &mz, 0, 0, 0, 0, 0, 0, 0);CHKERRQ(ierr);
  //ierr = DAGetCorners(PUNQ3.SimDA, &xs, &ys, &zs, &xm, &ym, &zm);CHKERRQ(ierr);

  //ierr = DAVecGetArray(PUNQ3.SimDA, PUNQ3.FlowMask, &LocalFlowMask);CHKERRQ(ierr);

  ///* set some geometry */
  //ierr = VecSet(PUNQ3.A1, 1.0);CHKERRQ(ierr);CHKMEMQ;
  //ierr = VecSet(PUNQ3.A2, 10.0);CHKERRQ(ierr);CHKMEMQ;
  //ierr = VecSet(PUNQ3.A3, 10.0);CHKERRQ(ierr);CHKMEMQ;
  //ierr = VecSet(PUNQ3.h1, 10.0);CHKERRQ(ierr);CHKMEMQ;
  //ierr = VecSet(PUNQ3.h2, 1.0);CHKERRQ(ierr);CHKMEMQ;
  //ierr = VecSet(PUNQ3.h3, 1.0);CHKERRQ(ierr);CHKMEMQ;
  ///* set initial values to vectors */
  //ierr = VecSet(PUNQ3.Po, 10.0);CHKERRQ(ierr);CHKMEMQ;
  ///* set residual oil and connate water */
  //PUNQ3.Sor = 0.2;
  //PUNQ3.Swc = 0.2;
  ///* Set saturations */
  //ierr = VecSet(PUNQ3.So, 1.0-0.2 );CHKERRQ(ierr);CHKMEMQ;
  //ierr = VecSet(PUNQ3.Sw, 0.2 );CHKERRQ(ierr);CHKMEMQ;
  ///* Set densities */
  //PUNQ3.Rhoos = 1.0;
  //PUNQ3.Rhows = 1.0;
  //ierr = VecSet(PUNQ3.Rhoo, PUNQ3.Rhoos );CHKERRQ(ierr);CHKMEMQ;
  //ierr = VecSet(PUNQ3.Rhow, PUNQ3.Rhows );CHKERRQ(ierr);CHKMEMQ;
  ///* Set viscosities */
  //ierr = VecSet(PUNQ3.Muo, 1);CHKERRQ(ierr);CHKMEMQ;
  //ierr = VecSet(PUNQ3.Muw, 1);CHKERRQ(ierr);CHKMEMQ;
  ///* Set permeabilities */
  //ierr = VecSet(PUNQ3.K11, 300.0);CHKERRQ(ierr);CHKMEMQ;
  //ierr = VecSet(PUNQ3.K22, 300.0);CHKERRQ(ierr);CHKMEMQ;
  //ierr = VecSet(PUNQ3.K33, 300.0);CHKERRQ(ierr);CHKMEMQ;
  ///* set the relative permeabilities */
  //ierr = VecSet(PUNQ3.Kro, 1.0);CHKERRQ(ierr);CHKMEMQ;
  //ierr = VecSet(PUNQ3.Krw, 1.0);CHKERRQ(ierr);CHKMEMQ;
  ///* Set the table for PcowSw */
  //PUNQ3.PcowSw.NumberOfEntries = 2;
  //ierr = PetscMalloc(PUNQ3.PcowSw.NumberOfEntries*sizeof(PetscScalar),&PUNQ3.PcowSw.X);CHKERRQ(ierr);CHKMEMQ;
  //ierr = PetscMalloc(PUNQ3.PcowSw.NumberOfEntries*sizeof(PetscScalar),&PUNQ3.PcowSw.Y);CHKERRQ(ierr);CHKMEMQ;
  //PUNQ3.PcowSw.X[0] = 0;PUNQ3.PcowSw.X[1] = 0;
  //PUNQ3.PcowSw.Y[0] = 0;PUNQ3.PcowSw.Y[1] = 1;
  ///* Set the graviational acceleration */
  //PUNQ3.GravAcc = 32.2;
  ///* Set the rock properties */
  //ierr = VecSet(PUNQ3.Phi, 0.2);CHKERRQ(ierr);CHKMEMQ;
  //PUNQ3.RockCompressibility = 0.0;
  //PUNQ3.RockCompRefPressure = 14.7;
  //PUNQ3.RockCompRefPorosity = 0.3;
  ///* Set the fluid compressibilities */
  //PUNQ3.OilCompressibility = 1.0e-5;
  //PUNQ3.OilCompPressure = 14.7;
  //PUNQ3.WaterCompressibility = 1.0e-5;
  //PUNQ3.WaterCompPressure = 14.7;
  ///* FIXME Set the volume factors */
  //ierr = VecSet(PUNQ3.Bo, 1.0);CHKERRQ(ierr);CHKMEMQ;
  //ierr = VecSet(PUNQ3.Bw, 1.0);CHKERRQ(ierr);CHKMEMQ;
  ///* Set the well information */
  ///* allocate wells */
  //PUNQ3.NumberOfWells = 3;
  //ierr = PetscMalloc(PUNQ3.NumberOfWells*sizeof(Well),&PUNQ3.Wells);CHKERRQ(ierr);CHKMEMQ;
  ///* allocate perforations for wells */
  //PUNQ3.Wells[0].NumberOfPerforations = 1;
  //ierr = PetscMalloc(PUNQ3.Wells[0].NumberOfPerforations*sizeof(Perforation),&PUNQ3.Wells[0].Perforations);CHKERRQ(ierr);CHKMEMQ;
  //PUNQ3.Wells[1].NumberOfPerforations = 1;
  //ierr = PetscMalloc(PUNQ3.Wells[1].NumberOfPerforations*sizeof(Perforation),&PUNQ3.Wells[1].Perforations);CHKERRQ(ierr);CHKMEMQ;
  //PUNQ3.Wells[2].NumberOfPerforations = 1;
  //ierr = PetscMalloc(PUNQ3.Wells[2].NumberOfPerforations*sizeof(Perforation),&PUNQ3.Wells[2].Perforations);CHKERRQ(ierr);CHKMEMQ;
  ///* Set information for well zero,the injector */
  //PUNQ3.Wells[0].Perforations[0].I = 1;
  //PUNQ3.Wells[0].Perforations[0].J = 1;
  //PUNQ3.Wells[0].Perforations[0].K = 1;
  //PUNQ3.Wells[0].Perforations[0].Qw = 1.0;
  //PUNQ3.Wells[0].Perforations[0].Constraint = FLOW_RATE_CONSTRAINT;
  //PUNQ3.Wells[0].Perforations[0].WellType = WATER_INJECTOR;
  //PUNQ3.Wells[0].Perforations[0].IsActive = PETSC_TRUE;
  //PUNQ3.Wells[0].Perforations[0].Orientation = PERF_ORIENTATION_X3X3;
  //PUNQ3.Wells[0].Perforations[0].Rw = 0.3;
  //PUNQ3.Wells[0].Perforations[0].S = 0.25;
  //PUNQ3.Wells[0].Perforations[0].zbh = 0.0;

  ///* Set information for well one,the producer */
  //PUNQ3.Wells[1].Perforations[0].I = mx - 2;
  //PUNQ3.Wells[1].Perforations[0].J = my - 2;
  //PUNQ3.Wells[1].Perforations[0].K = mz - 2;
  //PUNQ3.Wells[1].Perforations[0].BHPw = 5;
  //PUNQ3.Wells[1].Perforations[0].BHPo = 5;
  //PUNQ3.Wells[1].Perforations[0].Constraint = BHP_CONSTRAINT;
  //PUNQ3.Wells[1].Perforations[0].IsActive = PETSC_TRUE;
  //PUNQ3.Wells[1].Perforations[0].Orientation = PERF_ORIENTATION_X3X3;
  //PUNQ3.Wells[1].Perforations[0].Rw = 0.3;
  //PUNQ3.Wells[1].Perforations[0].S = 0.25;
  //PUNQ3.Wells[1].Perforations[0].zbh = 0.0;


  ///* Set information for well two,the monitoring well */
  //PUNQ3.Wells[2].Perforations[0].I = (mx - 2)/2;
  //PUNQ3.Wells[2].Perforations[0].J = (my - 2)/2;
  //PUNQ3.Wells[2].Perforations[0].K = (mz - 2)/2;
  //PUNQ3.Wells[2].Perforations[0].Constraint = MONITORING;
  //PUNQ3.Wells[2].Perforations[0].IsActive = PETSC_TRUE;
  //PUNQ3.Wells[2].Perforations[0].Orientation = PERF_ORIENTATION_X3X3;

  ///* Not do adaptive time-step */
  //PUNQ3.AdaptiveTimeStep = PETSC_FALSE;

  ///* Check if we perturb */
  //ierr = PetscOptionsGetTruth(PETSC_NULL, "-perturb", &BoolPerturb, PETSC_NULL);CHKERRQ(ierr);
  ///* Check if we run the perturbed case to completion */
  //ierr = PetscOptionsGetTruth(PETSC_NULL, "-runperturbed", &BoolRunPerturbed, PETSC_NULL);CHKERRQ(ierr);
  ///* Check if we run normally */
  //ierr = PetscOptionsGetTruth(PETSC_NULL, "-runnormal", &BoolRunNormal, PETSC_NULL);CHKERRQ(ierr);

  //if(BoolPerturb)
  //{
    //ierr = PetscOptionsGetInt(PETSC_NULL,"-seed_phi",&MySeed1,PETSC_NULL);CHKERRQ(ierr);
    //ierr = PetscOptionsGetInt(PETSC_NULL,"-seed_k11",&MySeed2,PETSC_NULL);CHKERRQ(ierr);
    //ierr = PetscOptionsGetInt(PETSC_NULL,"-seed_k22",&MySeed3,PETSC_NULL);CHKERRQ(ierr);
    //ierr = PetscOptionsGetInt(PETSC_NULL,"-seed_k33",&MySeed4,PETSC_NULL);CHKERRQ(ierr);
    //ierr = PetscOptionsGetInt(PETSC_NULL,"-seed_flowmask",&MySeed5,PETSC_NULL);CHKERRQ(ierr);

    //ierr = PetscOptionsGetReal(PETSC_NULL,"-percentage_phi",&Percentage1,PETSC_NULL);CHKERRQ(ierr);
    //ierr = PetscOptionsGetReal(PETSC_NULL,"-percentage_k11",&Percentage2,PETSC_NULL);CHKERRQ(ierr);
    //ierr = PetscOptionsGetReal(PETSC_NULL,"-percentage_k22",&Percentage3,PETSC_NULL);CHKERRQ(ierr);
    //ierr = PetscOptionsGetReal(PETSC_NULL,"-percentage_k33",&Percentage4,PETSC_NULL);CHKERRQ(ierr);
    //ierr = PetscOptionsGetReal(PETSC_NULL,"-percentage_flowmask",&Percentage5,PETSC_NULL);CHKERRQ(ierr);

    ///* create the random vector */
    //ierr = VecDuplicate(PUNQ3.Phi,&RandomVec);CHKERRQ(ierr);

    ///* perturb the porosity */
    //ierr = VecDuplicate(PUNQ3.Phi,&TempVec);CHKERRQ(ierr);
    //ierr = VecCopy(PUNQ3.Phi, TempVec);CHKERRQ(ierr);
    //ierr = VecScale(TempVec, Percentage1);CHKERRQ(ierr);
    //ierr = PetscRandomCreate(PETSC_COMM_WORLD,&rctx);CHKERRQ(ierr);
    //ierr = PetscRandomSetSeed(rctx,MySeed1);CHKERRQ(ierr);
    //ierr = PetscRandomSetFromOptions(rctx);CHKERRQ(ierr);
    //ierr = VecSetRandom(RandomVec,rctx);CHKERRQ(ierr);
    //ierr = PetscRandomDestroy(rctx);CHKERRQ(ierr);
    //ierr = VecPointwiseMult(TempVec,TempVec,RandomVec);CHKERRQ(ierr);
    //TempScalar = 1.0;
    //ierr = VecAYPX(PUNQ3.Phi,TempScalar,TempVec);CHKERRQ(ierr);
    //ierr = VecAssemblyBegin(PUNQ3.Phi);CHKERRQ(ierr);
    //ierr = VecAssemblyEnd(PUNQ3.Phi);CHKERRQ(ierr);

    ///* perturb the permeability along the x11 direction */
    //ierr = VecLog(PUNQ3.K11);CHKERRQ(ierr);
    //ierr = VecCopy(PUNQ3.K11, TempVec);CHKERRQ(ierr);
    //ierr = VecScale(TempVec, Percentage2);CHKERRQ(ierr);
    //ierr = PetscRandomCreate(PETSC_COMM_WORLD,&rctx);CHKERRQ(ierr);
    //ierr = PetscRandomSetSeed(rctx,MySeed2);CHKERRQ(ierr);
    //ierr = PetscRandomSetFromOptions(rctx);CHKERRQ(ierr);
    //ierr = VecSetRandom(RandomVec,rctx);CHKERRQ(ierr);
    //ierr = PetscRandomDestroy(rctx);CHKERRQ(ierr);
    //ierr = VecPointwiseMult(TempVec,TempVec,RandomVec);CHKERRQ(ierr);
    //TempScalar = 1.0;
    //ierr = VecAYPX(PUNQ3.K11,TempScalar,TempVec);CHKERRQ(ierr);
    //ierr = VecExp(PUNQ3.K11);CHKERRQ(ierr);
    //ierr = VecAssemblyBegin(PUNQ3.K11);CHKERRQ(ierr);
    //ierr = VecAssemblyEnd(PUNQ3.K11);CHKERRQ(ierr);

    ///* perturb the permeability along the x22 direction */
    //ierr = VecLog(PUNQ3.K22);CHKERRQ(ierr);
    //ierr = VecCopy(PUNQ3.K22, TempVec);CHKERRQ(ierr);
    //ierr = VecScale(TempVec, Percentage3);CHKERRQ(ierr);
    //ierr = PetscRandomCreate(PETSC_COMM_WORLD,&rctx);CHKERRQ(ierr);
    //ierr = PetscRandomSetSeed(rctx,MySeed3);CHKERRQ(ierr);
    //ierr = PetscRandomSetFromOptions(rctx);CHKERRQ(ierr);
    //ierr = VecSetRandom(RandomVec,rctx);CHKERRQ(ierr);
    //ierr = PetscRandomDestroy(rctx);CHKERRQ(ierr);
    //ierr = VecPointwiseMult(TempVec,TempVec,RandomVec);CHKERRQ(ierr);
    //TempScalar = 1.0;
    //ierr = VecAYPX(PUNQ3.K22,TempScalar,TempVec);CHKERRQ(ierr);
    //ierr = VecExp(PUNQ3.K22);CHKERRQ(ierr);
    //ierr = VecAssemblyBegin(PUNQ3.K22);CHKERRQ(ierr);
    //ierr = VecAssemblyEnd(PUNQ3.K22);CHKERRQ(ierr);

    ///* perturb the permeability along the x33 direction */
    //ierr = VecLog(PUNQ3.K33);CHKERRQ(ierr);
    //ierr = VecCopy(PUNQ3.K33, TempVec);CHKERRQ(ierr);
    //ierr = VecScale(TempVec, Percentage4);CHKERRQ(ierr);
    //ierr = PetscRandomCreate(PETSC_COMM_WORLD,&rctx);CHKERRQ(ierr);
    //ierr = PetscRandomSetSeed(rctx,MySeed4);CHKERRQ(ierr);
    //ierr = PetscRandomSetFromOptions(rctx);CHKERRQ(ierr);
    //ierr = VecSetRandom(RandomVec,rctx);CHKERRQ(ierr);
    //ierr = PetscRandomDestroy(rctx);CHKERRQ(ierr);
    //ierr = VecPointwiseMult(TempVec,TempVec,RandomVec);CHKERRQ(ierr);
    //TempScalar = 1.0;
    //ierr = VecAYPX(PUNQ3.K33,TempScalar,TempVec);CHKERRQ(ierr);
    //ierr = VecExp(PUNQ3.K33);CHKERRQ(ierr);
    //ierr = VecAssemblyBegin(PUNQ3.K33);CHKERRQ(ierr);
    //ierr = VecAssemblyEnd(PUNQ3.K33);CHKERRQ(ierr);

    //ierr = VecDestroy(TempVec);CHKERRQ(ierr);
    //ierr = VecDestroy(RandomVec);CHKERRQ(ierr);

    ///* perturb the flow mask */
    //ierr = VecSet(PUNQ3.FlowMask,FLUID_FLOW);CHKERRQ(ierr);
    //ierr = VecLog(PUNQ3.FlowMask);CHKERRQ(ierr);
    //ierr = VecCopy(PUNQ3.FlowMask, TempVec);CHKERRQ(ierr);
    //ierr = VecScale(TempVec, Percentage4);CHKERRQ(ierr);
    //ierr = PetscRandomCreate(PETSC_COMM_WORLD,&rctx);CHKERRQ(ierr);
    //ierr = PetscRandomSetSeed(rctx,MySeed4);CHKERRQ(ierr);
    //ierr = PetscRandomSetFromOptions(rctx);CHKERRQ(ierr);
    //ierr = VecSetRandom(RandomVec,rctx);CHKERRQ(ierr);
    //ierr = PetscRandomDestroy(rctx);CHKERRQ(ierr);
    //ierr = VecPointwiseMult(TempVec,TempVec,RandomVec);CHKERRQ(ierr);
    //TempScalar = 1.0;
    //ierr = VecAYPX(PUNQ3.K33,TempScalar,TempVec);CHKERRQ(ierr);
    //ierr = VecExp(PUNQ3.K33);CHKERRQ(ierr);
    //ierr = VecAssemblyBegin(PUNQ3.K33);CHKERRQ(ierr);
    //ierr = VecAssemblyEnd(PUNQ3.K33);CHKERRQ(ierr);

    //ierr = VecDestroy(TempVec);CHKERRQ(ierr);
    //ierr = VecDestroy(RandomVec);CHKERRQ(ierr);

    //ierr = PetscRandomCreate(PETSC_COMM_WORLD,&rctx);
    //for (k = zs; k < zs + zm; k++) {
      //for (j = ys; j < ys + ym; j++) {
        //for (i = xs; i < xs + xm; i++) {
          //if (i == 0 || j == 0 || k == 0 || i == mx - 1 || j == my - 1 || k == mz - 1) {
          //} else {
    
            //ierr = PetscRandomGetValue(rctx,&value);
            //if (value <= Percentage5)
              //LocalFlowMask[k][j][i] = NO_FLUID_FLOW;
          //}
        //}
      //}
    //}
    //ierr = DAVecRestoreArray(PUNQ3.SimDA, PUNQ3.FlowMask, &LocalFlowMask);CHKERRQ(ierr);
    //ierr = VecAssemblyBegin(PUNQ3.FlowMask);CHKERRQ(ierr);
    //ierr = VecAssemblyEnd(PUNQ3.FlowMask);CHKERRQ(ierr);

    //ierr = PetscRandomDestroy(rctx);
      
    //ierr = DefiantBlackOil2PhValiantWriteVecs(&PUNQ3);CHKERRQ(ierr);CHKMEMQ;
  //}

  //if (BoolRunPerturbed)
  //{

    //ierr = PetscPrintf(PETSC_COMM_WORLD,"Post Perturbed Phi is: with seed: %d \n", MySeed1);
    //ierr = VecView(PUNQ3.Phi,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
    //ierr = PetscPrintf(PETSC_COMM_WORLD,"Post Perturbed K11 is: with seed: %d \n", MySeed2);
    //ierr = VecView(PUNQ3.K11,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;

    ///* Set the Times */
    //PUNQ3.DeltaTP = 0.5;
    //PUNQ3.DeltaTS = 0.1;
    //PUNQ3.DSmax = 0.1;

    ///* Set the start and end time */
    //PUNQ3.StartTime = 0.0;
    //PUNQ3.EndTime = 20.0;

    ///* Now we solve */
    //ierr = DefiantIMPES2PhIterate(&PUNQ3);CHKERRQ(ierr);CHKMEMQ;
    ///* Now we write the fields */
    //ierr = DefiantBlackOil2PhValiantWriteVecs(&PUNQ3);CHKERRQ(ierr);CHKMEMQ;

    //ierr = PetscPrintf(PETSC_COMM_WORLD,"After perturbed run Sw is: \n");
    //ierr = VecView(PUNQ3.Sw,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
    //ierr = PetscPrintf(PETSC_COMM_WORLD,"After perturbed run x1 is: \n");
    //ierr = VecView(PUNQ3.x1,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
  //} else  if  (BoolRunNormal) {
    ///* Set the Times */
    //PUNQ3.DeltaTP = 0.5;
    //PUNQ3.DeltaTS = 0.1;
    //PUNQ3.DSmax = 0.1;

    ///* Set the start and end time */
    //PUNQ3.StartTime = 0.0;
    //PUNQ3.EndTime = 1.0;

    ///* we first load the modified fields */
    //ierr = DefiantBlackOil2PhValiantLoadVecs(&PUNQ3);CHKERRQ(ierr);CHKMEMQ;

    //ierr = PetscPrintf(PETSC_COMM_WORLD,"Post Perturbed Phi is: with seed: %d \n", MySeed1);
    //ierr = VecView(PUNQ3.Phi,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
    //ierr = PetscPrintf(PETSC_COMM_WORLD,"Post Perturbed K11 is: with seed: %d \n", MySeed2);
    //ierr = VecView(PUNQ3.K11,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;

    ///* Now we solve */
    //ierr = DefiantIMPES2PhIterate(&PUNQ3);CHKERRQ(ierr);CHKMEMQ;
    ///* Now we write the fields */
    //ierr = DefiantBlackOil2PhValiantWriteVecs(&PUNQ3);CHKERRQ(ierr);CHKMEMQ;

    //ierr = PetscPrintf(PETSC_COMM_WORLD,"Sw is: \n");
    //ierr = VecView(PUNQ3.Sw,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
    //ierr = PetscPrintf(PETSC_COMM_WORLD,"x1 is: \n");
    //ierr = VecView(PUNQ3.x1,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
  //}

  ///* Now we clean up */
  //ierr = DefiantDestroySimulationVecs(&PUNQ3);CHKERRQ(ierr);CHKMEMQ;
  //ierr = DefiantDestroy2PhIMPESSystem(&PUNQ3);CHKERRQ(ierr);CHKMEMQ;


  //PetscFunctionReturn(0);
//}
