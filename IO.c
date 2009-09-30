/*
 * DefiantIO.c
 *
 *  Created on: Sep 6, 2009
 *      Author: yye00
 */

#include "Defiant.h"

#undef __FUNCT__
#define __FUNCT__ "DefiantViewPressuresSTDOUT"
PetscErrorCode DefiantViewPressuresSTDOUT(BlackOilReservoirSimulation* MySim) {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = VecAssemblyBegin(MySim->Po);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->Pw);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->Pg);CHKERRQ(ierr);

  ierr = VecAssemblyEnd(MySim->Po);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->Pw);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->Pg);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "The value of Po is:\n");CHKERRQ(ierr);
  ierr = VecView(MySim->Po, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "The value of Pw is:\n");CHKERRQ(ierr);
  ierr = VecView(MySim->Pw, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "The value of Pg is:\n");CHKERRQ(ierr);
  ierr = VecView(MySim->Pg, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DefiantViewPressuresASCII"
PetscErrorCode DefiantViewPressuresASCII(BlackOilReservoirSimulation* MySim) {
  PetscErrorCode ierr;
  PetscViewer MyViewer;

  PetscFunctionBegin;
  ierr = VecAssemblyBegin(MySim->Po);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->Pw);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->Pg);CHKERRQ(ierr);

  ierr = VecAssemblyEnd(MySim->Po);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->Pw);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->Pg);CHKERRQ(ierr);

  ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, "Po.m", &MyViewer);CHKERRQ(ierr);
  ierr = PetscViewerSetFormat(MyViewer, PETSC_VIEWER_ASCII_MATLAB);CHKERRQ(ierr);
  ierr = VecView(MySim->Po, MyViewer);CHKERRQ(ierr);

  ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, "Pw.m", &MyViewer);CHKERRQ(ierr);
  ierr = PetscViewerSetFormat(MyViewer, PETSC_VIEWER_ASCII_MATLAB);CHKERRQ(ierr);
  ierr = VecView(MySim->Pw, MyViewer);CHKERRQ(ierr);

  ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, "Pg.m", &MyViewer);CHKERRQ(ierr);
  ierr = PetscViewerSetFormat(MyViewer, PETSC_VIEWER_ASCII_MATLAB);CHKERRQ(ierr);
  ierr = VecView(MySim->Pg, MyViewer);CHKERRQ(ierr);

  ierr = PetscViewerDestroy(MyViewer);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DefiantViewPressuresSTDOUT"
PetscErrorCode DefiantViewSaturationsSTDOUT(BlackOilReservoirSimulation* MySim) {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = VecAssemblyBegin(MySim->So);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->Sw);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->Sg);CHKERRQ(ierr);

  ierr = VecAssemblyEnd(MySim->So);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->Sw);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->Sg);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "The value of So is:\n");CHKERRQ(ierr);
  ierr = VecView(MySim->So, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "The value of Sw is:\n");CHKERRQ(ierr);
  ierr = VecView(MySim->Sw, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "The value of Sg is:\n");CHKERRQ(ierr);
  ierr = VecView(MySim->Sg, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DefiantViewSaturationsASCII"
PetscErrorCode DefiantViewSaturationsASCII(BlackOilReservoirSimulation* MySim) {
  PetscErrorCode ierr;
  PetscViewer MyViewer;

  PetscFunctionBegin;
  ierr = VecAssemblyBegin(MySim->So);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->Sw);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->Sg);CHKERRQ(ierr);

  ierr = VecAssemblyEnd(MySim->So);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->Sw);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->Sg);CHKERRQ(ierr);

  ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, "So.m", &MyViewer);CHKERRQ(ierr);
  ierr = PetscViewerSetFormat(MyViewer, PETSC_VIEWER_ASCII_MATLAB);CHKERRQ(ierr);
  ierr = VecView(MySim->So, MyViewer);CHKERRQ(ierr);

  ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, "Sw.m", &MyViewer);CHKERRQ(ierr);
  ierr = PetscViewerSetFormat(MyViewer, PETSC_VIEWER_ASCII_MATLAB);CHKERRQ(ierr);
  ierr = VecView(MySim->Sw, MyViewer);CHKERRQ(ierr);

  ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, "Sg.m", &MyViewer);CHKERRQ(ierr);
  ierr = PetscViewerSetFormat(MyViewer, PETSC_VIEWER_ASCII_MATLAB);CHKERRQ(ierr);
  ierr = VecView(MySim->Sg, MyViewer);CHKERRQ(ierr);

  ierr = PetscViewerDestroy(MyViewer);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "DefiantViewCoords"
PetscErrorCode DefiantViewCoords(BlackOilReservoirSimulation* MySim) {
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = VecAssemblyBegin(MySim->x1);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->x2);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->x3);CHKERRQ(ierr);

  ierr = VecAssemblyEnd(MySim->x1);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->x2);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->x3);CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD, "The value of X1 is:\n");CHKERRQ(ierr);
  ierr = VecView(MySim->x1, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "The value of X2 is:\n");CHKERRQ(ierr);
  ierr = VecView(MySim->x2, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "The value of X3 is:\n");CHKERRQ(ierr);
  ierr = VecView(MySim->x3, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
