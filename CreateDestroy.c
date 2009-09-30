/*
 * CreateDestroy.c
 *
 *  Created on: Sep 6, 2009
 *      Author: yye00
 */

#include "Defiant.h"

#undef __FUNCT__
#define __FUNCT__ "DefiantIMPES2PhCreateSimulation"
PetscErrorCode DefiantCreateSimulationVecs(BlackOilReservoirSimulation* MySim) {

  PetscErrorCode ierr;
  PetscFunctionBegin;

  ierr = DMMGCreate(PETSC_COMM_WORLD, 3, PETSC_NULL, &(MySim->SimDMMG));CHKERRQ(ierr);
  ierr = DACreate3d(PETSC_COMM_WORLD, DA_NONPERIODIC, DA_STENCIL_STAR, -3, -3,
      -3, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, 1, 1, 0, 0, 0,
      &(MySim->SimDA));CHKERRQ(ierr);
  ierr = DMMGSetDM(MySim->SimDMMG, (DM) MySim->SimDA);CHKERRQ(ierr);CHKMEMQ;CHKMEMQ;

  /* Pressure Vectors */
  ierr = DACreateGlobalVector(MySim->SimDA, &(MySim->Po));CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Pw);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Pg);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Pb);CHKERRQ(ierr);CHKMEMQ;
  /* Capillary Pressure Vectors */
  ierr = VecDuplicate(MySim->Po, &MySim->Pcow);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Pcog);CHKERRQ(ierr);CHKMEMQ;
  /* Saturation Vectors */
  ierr = VecDuplicate(MySim->Po, &MySim->So);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Sw);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Sg);CHKERRQ(ierr);CHKMEMQ;
  /* Density Vectors */
  ierr = VecDuplicate(MySim->Po, &MySim->Rhoo);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Rhow);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Rhog);CHKERRQ(ierr);CHKMEMQ;
  /* Viscosity Vectors */
  ierr = VecDuplicate(MySim->Po, &MySim->Muo);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Muw);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Mug);CHKERRQ(ierr);CHKMEMQ;
  /* Permeabilities Vectors */
  ierr = VecDuplicate(MySim->Po, &MySim->K11);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->K22);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->K33);CHKERRQ(ierr);CHKMEMQ;
  /* Relative permeabilities Vectors */
  ierr = VecDuplicate(MySim->Po, &MySim->Kro);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Krw);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Krg);CHKERRQ(ierr);CHKMEMQ;
  /* Transmissibility Vectors */
  ierr = VecDuplicate(MySim->Po, &MySim->Tox1p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Tox2p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Tox3p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Tox1m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Tox2m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Tox3m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Twx1p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Twx2p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Twx3p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Twx1m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Twx2m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Twx3m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Tgx1p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Tgx2p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Tgx3p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Tgx1m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Tgx2m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Tgx3m);CHKERRQ(ierr);CHKMEMQ;
  /* Flow rates in standard conditions volumes at cell centers */
  ierr = VecDuplicate(MySim->Po, &MySim->Qo);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Qw);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Qg);CHKERRQ(ierr);CHKMEMQ;
  /* Volume factors at cell centers */
  ierr = VecDuplicate(MySim->Po, &MySim->Bo);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Bw);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Bg);CHKERRQ(ierr);CHKMEMQ;
  /* Volume factors at cell faces */
  ierr = VecDuplicate(MySim->Po, &MySim->Box1p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Box2p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Box3p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Box1m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Box2m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Box3m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Bwx1p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Bwx2p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Bwx3p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Bwx1m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Bwx2m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Bwx3m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Bgx1p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Bgx2p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Bgx3p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Bgx1m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Bgx2m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Bgx3m);CHKERRQ(ierr);CHKMEMQ;
  /* The gas solubility */
  ierr = VecDuplicate(MySim->Po, &MySim->Rso);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Rsox1p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Rsox1m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Rsox2p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Rsox2m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Rsox3p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Rsox3m);CHKERRQ(ierr);CHKMEMQ;
  /* Area*Perm divided by height, transmissibility components */
  ierr = VecDuplicate(MySim->Po, &MySim->AKHx1m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->AKHx2m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->AKHx3m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->AKHx1p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->AKHx2p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->AKHx3p);CHKERRQ(ierr);CHKMEMQ;
  /* Density At the block interfaces */
  ierr = VecDuplicate(MySim->Po, &MySim->Rhoox1m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Rhoox1p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Rhowx1m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Rhowx1p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Rhogx1m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Rhogx1p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Rhoox2m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Rhoox2p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Rhowx2m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Rhowx2p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Rhogx2m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Rhogx2p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Rhoox3m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Rhoox3p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Rhowx3m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Rhowx3p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Rhogx3m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Rhogx3p);CHKERRQ(ierr);CHKMEMQ;
  /* Viscosity At the block interfaces */
  ierr = VecDuplicate(MySim->Po, &MySim->Muox1m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Muox1p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Muwx1m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Muwx1p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Mugx1m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Mugx1p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Muox2m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Muox2p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Muwx2m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Muwx2p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Mugx2m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Mugx2p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Muox3m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Muox3p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Muwx3m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Muwx3p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Mugx3m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->Mugx3p);CHKERRQ(ierr);CHKMEMQ;
  /* Density by Viscosity at the faces */
  ierr = VecDuplicate(MySim->Po, &MySim->RhoByMuox1m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->RhoByMuox1p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->RhoByMuwx1m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->RhoByMuwx1p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->RhoByMugx1m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->RhoByMugx1p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->RhoByMuox2m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->RhoByMuox2p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->RhoByMuwx2m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->RhoByMuwx2p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->RhoByMugx2m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->RhoByMugx2p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->RhoByMuox3m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->RhoByMuox3p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->RhoByMuwx3m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->RhoByMuwx3p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->RhoByMugx3m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->RhoByMugx3p);CHKERRQ(ierr);CHKMEMQ;
  /* Relative Permeability across the faces */
  ierr = VecDuplicate(MySim->Po, &MySim->RelPermox1m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->RelPermox1p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->RelPermox2m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->RelPermox2p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->RelPermox3m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->RelPermox3p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->RelPermwx1m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->RelPermwx1p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->RelPermwx2m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->RelPermwx2p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->RelPermwx3m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->RelPermwx3p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->RelPermgx1m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->RelPermgx1p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->RelPermgx2m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->RelPermgx2p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->RelPermgx3m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->RelPermgx3p);CHKERRQ(ierr);CHKMEMQ;

  /* Geometry vectors */
  ierr = VecDuplicate(MySim->Po, &MySim->x1);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->x2);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->x3);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->A1);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->A2);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->A3);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->h1);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->h2);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MySim->Po, &MySim->h3);CHKERRQ(ierr);CHKMEMQ;
  /* Porosity Vector */
  ierr = VecDuplicate(MySim->Po, &MySim->Phi);CHKERRQ(ierr);CHKMEMQ;
  /* Flow-mask Vector */
  ierr = VecDuplicate(MySim->Po, &MySim->FlowMask);CHKERRQ(ierr);CHKMEMQ;

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DefiantCreate2PhIMPESSystem"
extern PetscErrorCode DefiantCreate2PhIMPESSystem(BlackOilReservoirSimulation* MySim)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;

  /* Gravity collected terms */
  ierr = VecDuplicate(MySim->Po, &MySim->Gravity);CHKERRQ(ierr);CHKMEMQ;
  /* Capillary Pressure collected term */
  ierr = VecDuplicate(MySim->Po, &MySim->CapillaryPressure);CHKERRQ(ierr);CHKMEMQ;
  /* Matrix and RHS */
  ierr = DAGetMatrix(MySim->SimDA, MATMPIAIJ, &MySim->A);CHKERRQ(ierr);CHKMEMQ;
  ierr = DACreateGlobalVector(MySim->SimDA, &(MySim->RHS));CHKERRQ(ierr);CHKMEMQ;
  ierr = KSPCreate(PETSC_COMM_WORLD, &MySim->ksp);CHKERRQ(ierr);CHKMEMQ;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DefiantCreate3PhIMPESSystem"
extern PetscErrorCode DefiantCreate3PhIMPESSystem(BlackOilReservoirSimulation* MySim)
{

  PetscFunctionBegin;

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DefiantCreate2PhNewtonRaphsonSystem"
extern PetscErrorCode DefiantCreate2PhNewtonRaphsonSystem(BlackOilReservoirSimulation* MySim)
{
  PetscErrorCode ierr;
  PetscInt       M,N,P,m,n,p;
  PetscFunctionBegin;
  /* retrieve info from normal da */
  ierr = DAGetInfo(MySim->SimDA, 0, &M, &N, &P, &m, &n, &p, 0, 0, 0, 0);CHKERRQ(ierr);
  /* setup the new 3DOF da */
  ierr = DACreate3d(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR,M,N,P,m,n,p,2,1,0,0,0, &MySim->NewtonRaphsonDA);CHKERRQ(ierr);
  /* Matrix and RHS */
  ierr = DAGetMatrix(MySim->NewtonRaphsonDA, MATMPIAIJ, &MySim->NewtonRaphsonF);CHKERRQ(ierr);CHKMEMQ;
  ierr = DACreateGlobalVector(MySim->NewtonRaphsonDA, &(MySim->NewtonRaphsonSolution));CHKERRQ(ierr);CHKMEMQ;
  ierr = DACreateGlobalVector(MySim->NewtonRaphsonDA, &(MySim->NewtonRaphsonRHS));CHKERRQ(ierr);CHKMEMQ;
  /* SNES */
  ierr = SNESCreate(PETSC_COMM_WORLD,&MySim->snes);CHKERRQ(ierr);CHKMEMQ;
  /* temporary vectors */
  ierr = DACreateGlobalVector(MySim->SimDA, &(MySim->RHS0));CHKERRQ(ierr);CHKMEMQ;
  ierr = DACreateGlobalVector(MySim->SimDA, &(MySim->RHS1));CHKERRQ(ierr);CHKMEMQ;
  ierr = DACreateGlobalVector(MySim->SimDA, &(MySim->NRCapPressure));CHKERRQ(ierr);CHKMEMQ;
  ierr = DACreateGlobalVector(MySim->SimDA, &(MySim->Grav0));CHKERRQ(ierr);CHKMEMQ;
  ierr = DACreateGlobalVector(MySim->SimDA, &(MySim->Grav1));CHKERRQ(ierr);CHKMEMQ;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DefiantCreate3PhNewtonRaphsonSystem"
extern PetscErrorCode DefiantCreate3PhNewtonRaphsonSystem(BlackOilReservoirSimulation* MySim)
{
  PetscErrorCode ierr;
  PetscInt       M,N,P,m,n,p;
  PetscFunctionBegin;
  /* retrieve info from normal da */
  ierr = DAGetInfo(MySim->SimDA, 0, &M, &N, &P, &m, &n, &p, 0, 0, 0, 0);CHKERRQ(ierr);
  /* setup the new 3DOF da */
  ierr = DACreate3d(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR,M,N,P,m,n,p,3,1,0,0,0, &MySim->NewtonRaphsonDA);CHKERRQ(ierr);
  /* Matrix and RHS */
  ierr = DAGetMatrix(MySim->NewtonRaphsonDA, MATMPIAIJ, &MySim->NewtonRaphsonF);CHKERRQ(ierr);CHKMEMQ;
  ierr = DACreateGlobalVector(MySim->NewtonRaphsonDA, &(MySim->NewtonRaphsonSolution));CHKERRQ(ierr);CHKMEMQ;
  ierr = DACreateGlobalVector(MySim->NewtonRaphsonDA, &(MySim->NewtonRaphsonRHS));CHKERRQ(ierr);CHKMEMQ;
  /* SNES */
  ierr = SNESCreate(PETSC_COMM_WORLD,&MySim->snes);CHKERRQ(ierr);CHKMEMQ;
  /* temporary vectors */
  ierr = DACreateGlobalVector(MySim->SimDA, &(MySim->RHS0));CHKERRQ(ierr);CHKMEMQ;
  ierr = DACreateGlobalVector(MySim->SimDA, &(MySim->RHS1));CHKERRQ(ierr);CHKMEMQ;
  ierr = DACreateGlobalVector(MySim->SimDA, &(MySim->RHS2));CHKERRQ(ierr);CHKMEMQ;
  ierr = DACreateGlobalVector(MySim->SimDA, &(MySim->NRCapPressure));CHKERRQ(ierr);CHKMEMQ;
  ierr = DACreateGlobalVector(MySim->SimDA, &(MySim->NRGasCapPressure));CHKERRQ(ierr);CHKMEMQ;
  ierr = DACreateGlobalVector(MySim->SimDA, &(MySim->Grav0));CHKERRQ(ierr);CHKMEMQ;
  ierr = DACreateGlobalVector(MySim->SimDA, &(MySim->Grav1));CHKERRQ(ierr);CHKMEMQ;
  ierr = DACreateGlobalVector(MySim->SimDA, &(MySim->Grav2));CHKERRQ(ierr);CHKMEMQ;
  /* variable temporary vectors */
  ierr = DACreateGlobalVector(MySim->SimDA, &(MySim->PhiSoByBoOld));CHKERRQ(ierr);CHKMEMQ;
  ierr = DACreateGlobalVector(MySim->SimDA, &(MySim->PhiSwByBwOld));CHKERRQ(ierr);CHKMEMQ;
  ierr = DACreateGlobalVector(MySim->SimDA, &(MySim->PhiSgByBgPlusRsoSoByBoOld));CHKERRQ(ierr);CHKMEMQ;
  /* The solution Vector */
  ierr = DACreateGlobalVector(MySim->NewtonRaphsonDA, &(MySim->NRSolution));CHKERRQ(ierr);CHKMEMQ;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DefiantDestroySimulationVecs"
PetscErrorCode DefiantDestroySimulationVecs(BlackOilReservoirSimulation* MySim) {
  PetscErrorCode ierr;
  PetscFunctionBegin;

  ierr = DMMGDestroy(MySim->SimDMMG);CHKERRQ(ierr);CHKMEMQ;
  ierr = DADestroy(MySim->SimDA);CHKERRQ(ierr);CHKMEMQ;
  /* Pressures */
  ierr = VecDestroy(MySim->Po);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Pw);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Pg);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Pb);CHKERRQ(ierr);CHKMEMQ;
  /* Capillary pressures */
  ierr = VecDestroy(MySim->Pcow);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Pcog);CHKERRQ(ierr);CHKMEMQ;
  /* Saturations */
  ierr = VecDestroy(MySim->So);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Sw);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Sg);CHKERRQ(ierr);CHKMEMQ;
  /* Densities */
  ierr = VecDestroy(MySim->Rhoo);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Rhow);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Rhog);CHKERRQ(ierr);CHKMEMQ;
  /* Viscosities */
  ierr = VecDestroy(MySim->Muo);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Muw);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Mug);CHKERRQ(ierr);CHKMEMQ;
  /* Permeabilities */
  ierr = VecDestroy(MySim->K11);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->K22);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->K33);CHKERRQ(ierr);CHKMEMQ;
  /* Relative permeabilities */
  ierr = VecDestroy(MySim->Kro);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Krw);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Krg);CHKERRQ(ierr);CHKMEMQ;
  /* Transmissibility at the faces */
  ierr = VecDestroy(MySim->Tox1p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Tox2p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Tox3p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Tox1m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Tox2m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Tox3m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Twx1p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Twx2p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Twx3p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Twx1m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Twx2m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Twx3m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Tgx1p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Tgx2p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Tgx3p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Tgx1m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Tgx2m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Tgx3m);CHKERRQ(ierr);CHKMEMQ;
  /* Flow rates in standard conditions volumes at cell centers */
  ierr = VecDestroy(MySim->Qo);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Qw);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Qg);CHKERRQ(ierr);CHKMEMQ;
  /* Volume factors at cell centers */
  ierr = VecDestroy(MySim->Bo);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Bw);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Bg);CHKERRQ(ierr);CHKMEMQ;
  /* Volume factors at cell faces */
  ierr = VecDestroy(MySim->Box1p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Box2p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Box3p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Box1m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Box2m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Box3m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Bwx1p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Bwx2p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Bwx3p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Bwx1m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Bwx2m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Bwx3m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Bgx1p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Bgx2p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Bgx3p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Bgx1m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Bgx2m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Bgx3m);CHKERRQ(ierr);CHKMEMQ;
  /* The Gas solubility */
  ierr = VecDestroy(MySim->Rso);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Rsox1p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Rsox1m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Rsox2p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Rsox2m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Rsox3p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Rsox3m);CHKERRQ(ierr);CHKMEMQ;
  /* Area multiplied by permeability divided by grid size at the faces */
  ierr = VecDestroy(MySim->AKHx1m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->AKHx2m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->AKHx3m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->AKHx1p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->AKHx2p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->AKHx3p);CHKERRQ(ierr);CHKMEMQ;
  /* Density At the block interfaces */
  ierr = VecDestroy(MySim->Rhoox1m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Rhoox1p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Rhowx1m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Rhowx1p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Rhogx1m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Rhogx1p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Rhoox2m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Rhoox2p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Rhowx2m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Rhowx2p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Rhogx2m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Rhogx2p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Rhoox3m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Rhoox3p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Rhowx3m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Rhowx3p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Rhogx3m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Rhogx3p);CHKERRQ(ierr);CHKMEMQ;
  /* Viscosity At the block interfaces */
  ierr = VecDestroy(MySim->Muox1m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Muox1p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Muwx1m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Muwx1p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Mugx1m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Mugx1p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Muox2m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Muox2p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Muwx2m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Muwx2p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Mugx2m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Mugx2p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Muox3m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Muox3p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Muwx3m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Muwx3p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Mugx3m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Mugx3p);CHKERRQ(ierr);CHKMEMQ;
  /* Density divided by viscosity at the interfaces */
  ierr = VecDestroy(MySim->RhoByMuox1m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->RhoByMuox1p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->RhoByMuwx1m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->RhoByMuwx1p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->RhoByMugx1m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->RhoByMugx1p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->RhoByMuox2m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->RhoByMuox2p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->RhoByMuwx2m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->RhoByMuwx2p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->RhoByMugx2m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->RhoByMugx2p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->RhoByMuox3m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->RhoByMuox3p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->RhoByMuwx3m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->RhoByMuwx3p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->RhoByMugx3m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->RhoByMugx3p);CHKERRQ(ierr);CHKMEMQ;
  /* Relative Perms at the block interfaces */
  ierr = VecDestroy(MySim->RelPermox1m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->RelPermox1p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->RelPermox2m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->RelPermox2p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->RelPermox3m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->RelPermox3p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->RelPermwx1m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->RelPermwx1p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->RelPermwx2m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->RelPermwx2p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->RelPermwx3m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->RelPermwx3p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->RelPermgx1m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->RelPermgx1p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->RelPermgx2m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->RelPermgx2p);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->RelPermgx3m);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->RelPermgx3p);CHKERRQ(ierr);CHKMEMQ;

  ierr = VecDestroy(MySim->x1);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->x2);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->x3);CHKERRQ(ierr);CHKMEMQ;

  ierr = VecDestroy(MySim->A1);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->A2);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->A3);CHKERRQ(ierr);CHKMEMQ;

  ierr = VecDestroy(MySim->h1);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->h2);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->h3);CHKERRQ(ierr);CHKMEMQ;

  ierr = VecDestroy(MySim->Phi);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->FlowMask);CHKERRQ(ierr);CHKMEMQ;

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DefiantDestroy2PhIMPESSystem"
extern PetscErrorCode DefiantDestroy2PhIMPESSystem(BlackOilReservoirSimulation* MySim)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  /* Gravity collected terms */
  ierr = VecDestroy(MySim->Gravity);CHKERRQ(ierr);CHKMEMQ;
  /* Capillary Pressure collected term */
  ierr = VecDestroy(MySim->CapillaryPressure);CHKERRQ(ierr);CHKMEMQ;
  /* Matrix and RHS */
  ierr = MatDestroy(MySim->A);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->RHS);CHKERRQ(ierr);CHKMEMQ;
  ierr = KSPDestroy(MySim->ksp);CHKERRQ(ierr);CHKMEMQ;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DefiantDestroy3PhIMPESSystem"
extern PetscErrorCode DefiantDestroy3PhIMPESSystem(BlackOilReservoirSimulation* MySim)
{
  PetscFunctionBegin;

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DefiantDestroy2PhNewtonRaphsonSystem"
extern PetscErrorCode DefiantDestroy2PhNewtonRaphsonSystem(BlackOilReservoirSimulation* MySim)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  /* Matrix and RHS */
  ierr = MatDestroy(MySim->NewtonRaphsonF);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->NewtonRaphsonRHS);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->NewtonRaphsonSolution);CHKERRQ(ierr);CHKMEMQ;
  ierr = SNESDestroy(MySim->snes);CHKERRQ(ierr);CHKMEMQ;
  /* temporary vectors */
  ierr = VecDestroy(MySim->RHS0);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->RHS1);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->NRCapPressure);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Grav0);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Grav1);CHKERRQ(ierr);CHKMEMQ;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DefiantDestroy3PhNewtonRaphsonSystem"
extern PetscErrorCode DefiantDestroy3PhNewtonRaphsonSystem(BlackOilReservoirSimulation* MySim)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  /* Matrix and RHS */
  ierr = MatDestroy(MySim->NewtonRaphsonF);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->NewtonRaphsonRHS);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->NewtonRaphsonSolution);CHKERRQ(ierr);CHKMEMQ;
  ierr = SNESDestroy(MySim->snes);CHKERRQ(ierr);CHKMEMQ;
  /* temporary vectors */
  ierr = VecDestroy(MySim->RHS0);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->RHS1);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->RHS2);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->NRCapPressure);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->NRGasCapPressure);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Grav0);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Grav1);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->Grav2);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->NRSolution);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->PhiSoByBoOld);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->PhiSwByBwOld);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDestroy(MySim->PhiSgByBgPlusRsoSoByBoOld);CHKERRQ(ierr);CHKMEMQ;
  PetscFunctionReturn(0);
}
