/*
 * Defiant.c
 *
 *  Created on: Aug 24, 2009
 *      Author: yye00
 */

/*
 This is Defiant, a soon-to-be compositional reservoir simulator
 */

#include "Defiant.h"

extern PetscErrorCode TestIterate();
#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **argv) {
  PetscErrorCode ierr;

  /* Initialize PETSc */
  PetscInitialize(&argc, &argv, (char *) 0, help);

  ierr = DefiantIMPES2PhFivePoint();CHKERRQ(ierr);
  //ierr = TestIterate();CHKERRQ(ierr);

  ierr = PetscFinalize();CHKERRQ(ierr);

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "ComputeRHS"
PetscErrorCode ComputeRHS(DMMG dmmg, Vec b) {
  PetscErrorCode ierr;
  PetscInt mx, my, mz;
  PetscScalar h;

  PetscFunctionBegin;
  ierr = DAGetInfo((DA) dmmg->dm, 0, &mx, &my, &mz, 0, 0, 0, 0, 0, 0, 0);CHKERRQ(ierr);
  h = 1.0 / ((mx - 1) * (my - 1) * (mz - 1));
  ierr = VecSet(b, h);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ComputeMatrix"
PetscErrorCode ComputeMatrix(DMMG dmmg, Mat jac, Mat B) {
  DA da = (DA) dmmg->dm;
  PetscErrorCode ierr;
  PetscInt i, j, k, mx, my, mz, xm, ym, zm, xs, ys, zs;
  PetscScalar v[7], Hx, Hy, Hz, HxHydHz, HyHzdHx, HxHzdHy;
  MatStencil row, col[7];

  ierr = DAGetInfo(da, 0, &mx, &my, &mz, 0, 0, 0, 0, 0, 0, 0);CHKERRQ(ierr);
  Hx = 1.0 / (PetscReal) (mx - 1);
  Hy = 1.0 / (PetscReal) (my - 1);
  Hz = 1.0 / (PetscReal) (mz - 1);
  HxHydHz = Hx * Hy / Hz;
  HxHzdHy = Hx * Hz / Hy;
  HyHzdHx = Hy * Hz / Hx;
  ierr = DAGetCorners(da, &xs, &ys, &zs, &xm, &ym, &zm);CHKERRQ(ierr);

  for (k = zs; k < zs + zm; k++) {
    for (j = ys; j < ys + ym; j++) {
      for (i = xs; i < xs + xm; i++) {
        row.i = i;
        row.j = j;
        row.k = k;
        if (i == 0 || j == 0 || k == 0 || i == mx - 1 || j == my - 1 || k == mz
            - 1) {
          v[0] = 2.0 * (HxHydHz + HxHzdHy + HyHzdHx);
          ierr = MatSetValuesStencil(B, 1, &row, 1, &row, v, INSERT_VALUES);
          CHKERRQ(ierr);
        } else {
          v[0] = -HxHydHz;
          col[0].i = i;
          col[0].j = j;
          col[0].k = k - 1;
          v[1] = -HxHzdHy;
          col[1].i = i;
          col[1].j = j - 1;
          col[1].k = k;
          v[2] = -HyHzdHx;
          col[2].i = i - 1;
          col[2].j = j;
          col[2].k = k;
          v[3] = 2.0 * (HxHydHz + HxHzdHy + HyHzdHx);
          col[3].i = row.i;
          col[3].j = row.j;
          col[3].k = row.k;
          v[4] = -HyHzdHx;
          col[4].i = i + 1;
          col[4].j = j;
          col[4].k = k;
          v[5] = -HxHzdHy;
          col[5].i = i;
          col[5].j = j + 1;
          col[5].k = k;
          v[6] = -HxHydHz;
          col[6].i = i;
          col[6].j = j;
          col[6].k = k + 1;
          ierr = MatSetValuesStencil(B, 1, &row, 7, col, v, INSERT_VALUES);
          CHKERRQ(ierr);
        }
      }
    }
  }
  ierr = MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "TestIterate"
PetscErrorCode TestIterate() {
  PetscErrorCode ierr;
  PetscInt i, j, k, mx, my, mz, xm, ym, zm, xs, ys, zs, rank;
  PetscScalar ***array;
  Vec Pooo;
  Vec Mooo;
  DA da;

  PetscFunctionBegin;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  ierr = DACreate3d(PETSC_COMM_WORLD, DA_NONPERIODIC, DA_STENCIL_STAR, -3, -3,
        -3, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, 1, 1, 0, 0, 0,
        &(da));CHKERRQ(ierr);
  ierr = DACreateGlobalVector(da, &(Pooo));CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(Pooo,&Mooo);CHKERRQ(ierr);


  ierr = DAView(da, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  ierr = DAGetInfo(da, 0, &mx, &my, &mz, 0, 0, 0, 0, 0, 0, 0);CHKERRQ(ierr);

  ierr = DAGetCorners(da, &xs, &ys, &zs, &xm, &ym, &zm);CHKERRQ(ierr);


  ierr = DAVecGetArray(da,Pooo, &array);CHKERRQ(ierr);


  fprintf(stderr, "\n size in x=%d  in y=%d in z=%d\n", xm, ym, zm);
  fprintf(
      stderr,
      "\n For rank: %d starting indices are as follows in x=%d  in y=%d in z=%d\n",
      rank, xs, ys, zs);

  for (k = zs; k < zs + zm; k++) {
    for (j = ys; j < ys + ym; j++) {
      for (i = xs; i < xs + xm; i++) {
        if (i == 0 || j == 0 || k == 0 || i == mx - 1 || j == my - 1 || k == mz
            - 1) {
          array[k][j][i] = -3.1415912654;
          PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Proc rank: %d at boundary location i=%d j=%d k=%d, array is:%f\n", rank, i,j,k, array[k][j][i]);
          ierr = PetscSynchronizedFlush(PETSC_COMM_WORLD);CHKERRQ(ierr);
        } else {
          array[k][j][i] = rank;
          PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Proc rank: %d at location i=%d j=%d k=%d, array is:%f\n", rank, i,j,k, array[k][j][i]);
          ierr = PetscSynchronizedFlush(PETSC_COMM_WORLD);CHKERRQ(ierr);

        }
      }
    }
  }

  ierr = DAVecRestoreArray(da, Pooo, &array);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(Pooo);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(Pooo);CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD, "\n Printing Poo \n");CHKERRQ(ierr);
  ierr = VecView(Pooo,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;

  ierr = DAGetLocalVector(da,&Mooo);CHKERRQ(ierr);
  ierr = DAGlobalToLocalBegin(da,Pooo,INSERT_VALUES,Mooo);CHKERRQ(ierr);
  ierr = DAGlobalToLocalEnd(da,Pooo,INSERT_VALUES,Mooo);CHKERRQ(ierr);

  ierr = DAVecGetArray(da,Mooo,&array);CHKERRQ(ierr);



  if (rank==0){
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Proc rank: %d array is:%f\n", rank, array[1][1][5]);
    ierr = PetscSynchronizedFlush(PETSC_COMM_WORLD);CHKERRQ(ierr);
  }

  if (rank==1){
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Proc rank: %d array is:%f\n", rank, array[1][1][4]);
    ierr = PetscSynchronizedFlush(PETSC_COMM_WORLD);CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

