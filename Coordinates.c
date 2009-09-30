/*
 * Coordinates.c
 *
 *  Created on: Sep 6, 2009
 *      Author: yye00
 */

#include "Defiant.h"

#undef __FUNCT__
#define __FUNCT__ "DefiantGetDACoords"
PetscErrorCode DefiantGetDACoords(BlackOilReservoirSimulation* MySim) {
  PetscErrorCode ierr;
  Vec Coords, vecs[3];

  PetscFunctionBegin;
  vecs[0] = MySim->x1;
  vecs[1] = MySim->x2;
  vecs[2] = MySim->x3;

  ierr = DAGetCoordinates(MySim->SimDA, &Coords);CHKERRQ(ierr);
  if (!Coords) {
    ierr = DASetUniformCoordinates(MySim->SimDA, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
    CHKERRQ(ierr);
    ierr = DAGetCoordinates(MySim->SimDA, &Coords);
    CHKERRQ(ierr);
  }
  ierr = VecStrideGatherAll(Coords, vecs, INSERT_VALUES);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DefiantSetDACoords"
PetscErrorCode DefiantSetDACoords(BlackOilReservoirSimulation* MySim) {
  PetscErrorCode ierr;
  Vec Coords, vecs[3];

  PetscFunctionBegin;

  /* Set the blocksize to one for the coordinate vectors */
  ierr = VecSetBlockSize(MySim->x1, 1);CHKERRQ(ierr);
  ierr = VecSetBlockSize(MySim->x2, 1);CHKERRQ(ierr);
  ierr = VecSetBlockSize(MySim->x3, 1);CHKERRQ(ierr);
  /* Set vecs from options */
  ierr = VecSetFromOptions(MySim->x1);CHKERRQ(ierr);
  ierr = VecSetFromOptions(MySim->x2);CHKERRQ(ierr);
  ierr = VecSetFromOptions(MySim->x3);CHKERRQ(ierr);

  vecs[0] = MySim->x1;
  vecs[1] = MySim->x2;
  vecs[2] = MySim->x3;

  ierr = DAGetCoordinates(MySim->SimDA, &Coords);CHKERRQ(ierr);
  ierr = VecStrideGatherAll(Coords, vecs, INSERT_VALUES);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
