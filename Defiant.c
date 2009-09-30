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

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **argv) {
  PetscErrorCode ierr;

  /* Initialize PETSc */
  PetscInitialize(&argc, &argv, (char *) 0, help);

  ierr = DefiantIMPES2PhBuckleyLeverett();CHKERRQ(ierr);

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
PetscErrorCode TestIterate(BlackOilReservoirSimulation* MySim) {
  PetscErrorCode ierr;
  PetscInt i, j, k, mx, my, mz, xm, ym, zm, xs, ys, zs, rank;
  PetscInt BlockSize;
  PetscScalar ***array;
  PetscScalar *array2;
  Vec Coords;

  PetscFunctionBegin;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  ierr = DAView(MySim->SimDA, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  ierr = DAGetInfo(MySim->SimDA, 0, &mx, &my, &mz, 0, 0, 0, 0, 0, 0, 0);CHKERRQ(ierr);

  ierr = DAGetCorners(MySim->SimDA, &xs, &ys, &zs, &xm, &ym, &zm);CHKERRQ(ierr);

  ierr = DAVecGetArray(MySim->SimDA, MySim->Po, &array);CHKERRQ(ierr);

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
          array[k][j][i] = -1.0;
          //fprintf(stderr,"Proc rank: %d at boundary location i=%d j=%d k=%d, array is:%f\n", rank, i,j,k, array[k][j][i]);
        } else {
          array[k][j][i] = rank;
          //fprintf(stderr,"Proc rank: %d at location i=%d j=%d k=%d, array is:%f\n", rank, i,j,k, array[k][j][i]);
        }
      }
    }
  }

  ierr = DAVecRestoreArray(MySim->SimDA, MySim->Po, &array);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MySim->Po);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MySim->Po);CHKERRQ(ierr);

  ierr = DAGetCoordinates(MySim->SimDA, &Coords);CHKERRQ(ierr);
  if (!Coords) {
    ierr
        = DASetUniformCoordinates(MySim->SimDA, 0.0, 1.0, -2.0, -1.0, 3.0, 5.0);
    CHKERRQ(ierr);
    ierr = DAGetCoordinates(MySim->SimDA, &Coords);
    CHKERRQ(ierr);
  }

  ierr = VecGetBlockSize(Coords, &BlockSize);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "\n blocksize for DA is:%d", BlockSize);CHKERRQ(ierr);

  ierr = VecGetArray(Coords, &array2);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "\nX_COORDINATES %d double\n", xm);CHKERRQ(ierr);
  for (i = 0; i < xm; i++) {
    ierr = PetscPrintf(PETSC_COMM_WORLD, "\n%G ", PetscRealPart(array2[i*3]));
    CHKERRQ(ierr);
  }
  ierr = PetscPrintf(PETSC_COMM_WORLD, "\nY_COORDINATES %d double\n", ym);CHKERRQ(ierr);
  for (i = 0; i < ym; i++) {
    ierr = PetscPrintf(PETSC_COMM_WORLD, "\n%G ",
        PetscRealPart(array2[i*xm*3+1]));
    CHKERRQ(ierr);
  }
  ierr = PetscPrintf(PETSC_COMM_WORLD, "\nZ_COORDINATES %d double\n", zm);CHKERRQ(ierr);
  for (i = 0; i < zm; i++) {
    ierr = PetscPrintf(PETSC_COMM_WORLD, "\n%G ",
        PetscRealPart(array2[i*xm*ym*3+2]));
    CHKERRQ(ierr);
  }

  k = 0;
  ierr = PetscPrintf(PETSC_COMM_WORLD, "\n Alternate way \n");CHKERRQ(ierr);
  for (i = 0; i < mx * my * mz; i++) {
    //    ierr = PetscPrintf(PETSC_COMM_WORLD, "\nX_COORDINATES %G double\n", PetscRealPart(array2[k]));CHKERRQ(ierr);
    //    ierr = PetscPrintf(PETSC_COMM_WORLD, "\nY_COORDINATES %G double\n", PetscRealPart(array2[k+1]));CHKERRQ(ierr);
    //    ierr = PetscPrintf(PETSC_COMM_WORLD, "\nZ_COORDINATES %G double\n", PetscRealPart(array2[k+2]));CHKERRQ(ierr);
    k = k + 3;
  }

  ierr = PetscPrintf(PETSC_COMM_WORLD, "\n Alternate way \n");CHKERRQ(ierr);
  //  ierr = VecView(Coords, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "\n Alternate way \n");CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

