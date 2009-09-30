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
  ierr = DASetUniformCoordinates(BuckleyLeverett.SimDA, 0.0, 90.0, 0.0, 3.0, 0.0, 3.0);CHKERRQ(ierr);CHKMEMQ;
  ierr = DefiantGetDACoords(&BuckleyLeverett);CHKERRQ(ierr);CHKMEMQ;
  /* set some geometry */
  ierr = VecSet(BuckleyLeverett.A1, 1.0);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSet(BuckleyLeverett.A2, 1.0);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSet(BuckleyLeverett.A3, 1.0);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSet(BuckleyLeverett.h1, 10.0);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSet(BuckleyLeverett.h2, 1.0);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSet(BuckleyLeverett.h3, 1.0);CHKERRQ(ierr);CHKMEMQ;
  /* set initial values to vectors */
  ierr = VecSet(BuckleyLeverett.Po, 4000.0);CHKERRQ(ierr);CHKMEMQ;
  /* set residual oil and connate water */
  BuckleyLeverett.Sor = 0.2;
  BuckleyLeverett.Swc = 0.2;
  /* Set saturations */
  ierr = VecSet(BuckleyLeverett.So, 1.0-BuckleyLeverett.Swc );CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSet(BuckleyLeverett.Sw,     BuckleyLeverett.Swc );CHKERRQ(ierr);CHKMEMQ;
  /* Set densities */
  BuckleyLeverett.Rhoos = 1.0;
  BuckleyLeverett.Rhows = 1.0;
  ierr = VecSet(BuckleyLeverett.Rhoo, BuckleyLeverett.Rhoos );CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSet(BuckleyLeverett.Rhow, BuckleyLeverett.Rhows );CHKERRQ(ierr);CHKMEMQ;
  /* Set viscosities */
  ierr = VecSet(BuckleyLeverett.Muo, 0.00001);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSet(BuckleyLeverett.Muw, 0.00001);CHKERRQ(ierr);CHKMEMQ;
  /* Set permeabilities */
  ierr = VecSet(BuckleyLeverett.K11, 300.0);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSet(BuckleyLeverett.K22, 300.0);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSet(BuckleyLeverett.K33, 300.0);CHKERRQ(ierr);CHKMEMQ;
  /* set the relative permeabilities */
  ierr = VecSet(BuckleyLeverett.Kro, 1.0);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSet(BuckleyLeverett.Krw, 1.0);CHKERRQ(ierr);CHKMEMQ;
  /* Set the table for PcowSw */
  BuckleyLeverett.PcowSw.NumberOfEntries = 2;
  ierr = PetscMalloc(BuckleyLeverett.PcowSw.NumberOfEntries*sizeof(PetscScalar),&BuckleyLeverett.PcowSw.X);
  ierr = PetscMalloc(BuckleyLeverett.PcowSw.NumberOfEntries*sizeof(PetscScalar),&BuckleyLeverett.PcowSw.Y);
  BuckleyLeverett.PcowSw.X[0] = 0;BuckleyLeverett.PcowSw.X[1] = 0;
  BuckleyLeverett.PcowSw.Y[0] = 0;BuckleyLeverett.PcowSw.Y[1] = 1;
  /* Set the graviational acceleration */
  BuckleyLeverett.GravAcc = 32.2;
  /* Set the rock properties */
  ierr = VecSet(BuckleyLeverett.Phi, 0.2);CHKERRQ(ierr);CHKMEMQ;
  BuckleyLeverett.RockCompressibility = 1.0e-5;
  BuckleyLeverett.RockCompRefPressure = 14.7;
  BuckleyLeverett.RockCompRefPorosity = 0.3;
  /* Set the Times */
  BuckleyLeverett.StartTime = 0.0;
  BuckleyLeverett.EndTime = 0.3;
  BuckleyLeverett.DeltaTP = 0.1;
  BuckleyLeverett.DeltaTS = 0.01;
  BuckleyLeverett.DSmax = 0.01;
  /* Set the fluid compressibilities */
  BuckleyLeverett.OilCompressibility = 1.0e-5;
  BuckleyLeverett.OilCompPressure = 14.7;
  BuckleyLeverett.WaterCompressibility = 1.0e-6;
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
  BuckleyLeverett.Wells[0].Perforations[0].Qw = 3.1e11;
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
  BuckleyLeverett.Wells[1].Perforations[0].BHPw = 3000;
  BuckleyLeverett.Wells[1].Perforations[0].BHPo = 3000;
  BuckleyLeverett.Wells[1].Perforations[0].Constraint = BHP_CONSTRAINT;
  BuckleyLeverett.Wells[1].Perforations[0].IsActive = PETSC_TRUE;
  BuckleyLeverett.Wells[1].Perforations[0].Orientation = PERF_ORIENTATION_X3X3;
  BuckleyLeverett.Wells[1].Perforations[0].Rw = 0.3;
  BuckleyLeverett.Wells[1].Perforations[0].S = 0.25;
  BuckleyLeverett.Wells[1].Perforations[0].zbh = 0.0;

  /* Now we solve */
  ierr = DefiantIMPES2PhIterate(&BuckleyLeverett);CHKERRQ(ierr);CHKMEMQ;

  /* Now we clean up */
  ierr = DefiantDestroySimulationVecs(&BuckleyLeverett);CHKERRQ(ierr);CHKMEMQ;
  ierr = DefiantDestroy2PhIMPESSystem(&BuckleyLeverett);CHKERRQ(ierr);CHKMEMQ;


  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DefiantIMPES3PhBuckleyLeverett"
extern PetscErrorCode DefiantIMPES3PhBuckleyLeverett(BlackOilReservoirSimulation* MySim);

#undef __FUNCT__
#define __FUNCT__ "DefiantImprovedIMPES2PhBuckleyLeverett"
extern PetscErrorCode DefiantImprovedIMPES2PhBuckleyLeverett(BlackOilReservoirSimulation* MySim);

#undef __FUNCT__
#define __FUNCT__ "DefiantImprovedIMPES3PhBuckleyLeverett"
extern PetscErrorCode DefiantImprovedIMPES3PhBuckleyLeverett(BlackOilReservoirSimulation* MySim);

#undef __FUNCT__
#define __FUNCT__ "DefiantNewton2PhBuckleyLeverett"
extern PetscErrorCode DefiantNewton2PhBuckleyLeverett(BlackOilReservoirSimulation* MySim);

#undef __FUNCT__
#define __FUNCT__ "DefiantNewton3PhBuckleyLeverett"
extern PetscErrorCode DefiantNewton3PhBuckleyLeverett(BlackOilReservoirSimulation* MySim);

#undef __FUNCT__
#define __FUNCT__ "DefiantIMPES2PhFivePoint"
extern PetscErrorCode DefiantIMPES2PhFivePoint(BlackOilReservoirSimulation* MySim);

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
