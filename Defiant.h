/*
 * Defiant.h
 *
 *  Created on: Aug 24, 2009
 *      Author: yye00
 */

#ifndef DEFIANT_H_
#define DEFIANT_H_

static char help[] = "Defiant, a soon-to-be compositional reservoir simulator.\n\n";

#include "petscda.h"
#include "petscksp.h"
#include "petscdmmg.h"

#define DEFIANT_DEBUG             0
#define DEFIANT_USE_X_VIEWER      0
#define DEFIANT_USE_ASCII_VIEWER  0
#define DEFIANT_USE_MATLAB_VIEWER 1

#if defined(PETSC_HAVE_MATLAB_ENGINE)
#undef DEFIANT_USE_MATLAB_VIEWER
#define DEFIANT_USE_MATLAB_VIEWER 1
#else
#undef DEFIANT_USE_MATLAB_VIEWER
#define DEFIANT_USE_MATLAB_VIEWER 0
#endif

/* FlowMask */
#define FLUID_FLOW    2063
#define NO_FLUID_FLOW 2064

/* Well orientation */
#define PERF_ORIENTATION_X1X1 1701
#define PERF_ORIENTATION_X2X2 1702
#define PERF_ORIENTATION_X3X3 1703
#define PERF_ORIENTATION_X1X2 1704
#define PERF_ORIENTATION_X2X3 1705
#define PERF_ORIENTATION_X3X1 1706


/* Some Math Stuff */
#define PI 3.14159265358979323846264338327
#define EPSILON 1e-30
#define  ABS(a) ((a) < 0.0 ? -(a) : (a))

/* Some tables */
typedef struct {
  PetscReal *X;    /* X values for the table */
  PetscReal *Y;    /* Y values for the table */
  PetscInt NumberOfEntries;
} DefiantTable;

/* Some table functions */
extern PetscErrorCode DefiantTableXFromY(DefiantTable *MyTable, PetscReal *X, PetscReal *Y);
extern PetscErrorCode DefiantTableYFromX(DefiantTable *MyTable, PetscReal *X, PetscReal *Y);

/* Well constraint type */
#define BHP_CONSTRAINT       1707
#define FLOW_RATE_CONSTRAINT 1708
#define MONITORING           1709

/* Well types */
#define OIL_INJECTOR         1710
#define WATER_INJECTOR       1711
#define GAS_INJECTOR         1712
#define OIL_PRODUCER         1713
#define WATER_PRODUCER       1714
#define GAS_PRODUCER         1715
#define PRODUCER             1716

typedef struct {
  /* Global DA indices */
  PetscInt I,J,K;
  /* The actual type of the well: injector or producer ? */
  PetscInt WellType;
  /* Is the perforation active */
  PetscTruth IsActive;
  /* Orientation of the perforation */
  PetscInt Orientation;
  /* What is the constraint at this perforation */
  PetscInt Constraint;
  /* What is the well-index at this perforation */
  PetscReal WellIndex;
  /* Perforation parameters */
  PetscReal Apo;
  PetscReal Apw;
  PetscReal Apg ;
  PetscReal Qo;    /* user set flow rate for flow rate constraint wells */
  PetscReal Qw;    /* user set flow rate for flow rate constraint wells */
  PetscReal Qg;    /* user set flow rate for flow rate constraint wells */
  PetscReal BHPo;  /* Perforation pressure for oil   */
  PetscReal BHPw;  /* Perforation pressure for water */
  PetscReal BHPg;  /* Perforation pressure for gas   */
  PetscReal zbh;   /* reference depth for above pressure */
  PetscReal Rw;    /* radius at perforation */
  PetscReal S;     /* Skin at perforation */
  /* Other perforation properties */
  PetscScalar So, Sw, Sg;
  PetscScalar Po, Pw, Pg;
  PetscScalar Kro, Krw, Krg;
  PetscScalar Muo, Muw, Mug;
  PetscScalar Rhoo, Rhow, Rhog;
  PetscScalar Pcow, Pcog;
  PetscScalar Bo, Bw, Bg;
  PetscScalar x1, x2, x3;
  PetscScalar h1, h2, h3;
  /* Owner core, need this for parallel */
  PetscInt OwnerRank;
} Perforation;

typedef struct {
  /* Number of Perforations */
  PetscInt NumberOfPerforations;
  /* Pointer to the Perforations array */
  Perforation *Perforations;
  /* Constraint type */
  PetscInt ConstraintType;
  /* Well parameters */
  PetscReal TotalQo;   /* cumulative oil flow rate, summed over all perfs user constraint   */
  PetscReal TotalQw;   /* cumulative water flow rate, summed over all perfs user constraint */
  PetscReal TotalQg;   /* cumulative gas flow rate, summed over all perfs user constraint   */
} Well;

typedef struct {
  /* The DMMG we will be using */
  DMMG  *SimDMMG;
  /* The distributed array structure we will be using */
  DA  SimDA;

  /* pressures */
  Vec Po, Pw, Pg, Pb;
  /* Capillary pressures */
  Vec Pcow, Pcog;
  /* Saturations */
  Vec So, Sw, Sg;
  /* Densities */
  Vec Rhoo, Rhow, Rhog;
  /* Viscosities */
  Vec Muo, Muw, Mug;
  /* Permeabilities */
  Vec K11, K22, K33;
  /* Relative permeabilities */
  Vec Kro, Krw, Krg;
  /* Transmissibility For Oil*/
  Vec Tox1p, Tox2p, Tox3p;
  Vec Tox1m, Tox2m, Tox3m;
  /* Transmissibility For Water*/
  Vec Twx1p, Twx2p, Twx3p;
  Vec Twx1m, Twx2m, Twx3m;
  /* Transmissibility For Gas*/
  Vec Tgx1p, Tgx2p, Tgx3p;
  Vec Tgx1m, Tgx2m, Tgx3m;
  /* Flow rates at cells */
  Vec Qo, Qw, Qg;
  /* Volume factors at cell center */
  Vec Bo, Bw, Bg;
  /* Volume vactors at cell faces for oil */
  Vec Box1p, Box2p, Box3p;
  Vec Box1m, Box2m, Box3m;
  /* Volume vactors at cell faces for water */
  Vec Bwx1p, Bwx2p, Bwx3p;
  Vec Bwx1m, Bwx2m, Bwx3m;
  /* Volume vactors at cell faces for gas */
  Vec Bgx1p, Bgx2p, Bgx3p;
  Vec Bgx1m, Bgx2m, Bgx3m;
  /* Gas solubility */
  Vec Rso;
  Vec Rsox1p, Rsox1m;
  Vec Rsox2p, Rsox2m;
  Vec Rsox3p, Rsox3m;


  /* Porosity */
  Vec Phi;
  /* Flow mask */
  Vec FlowMask;
  /* Coordinates, Areas and h123 */
  Vec x1, x2, x3;
  Vec A1, A2, A3;
  Vec h1, h2, h3;
  /* AKByH */
  Vec AKHx1m, AKHx1p;
  Vec AKHx2m, AKHx2p;
  Vec AKHx3m, AKHx3p;
  /* Density at the faces */
  Vec Rhoox1m, Rhoox1p, Rhowx1m, Rhowx1p, Rhogx1m, Rhogx1p;
  Vec Rhoox2m, Rhoox2p, Rhowx2m, Rhowx2p, Rhogx2m, Rhogx2p;
  Vec Rhoox3m, Rhoox3p, Rhowx3m, Rhowx3p, Rhogx3m, Rhogx3p;
  /* Viscosity at the faces */
  Vec Muox1m, Muox1p, Muwx1m, Muwx1p, Mugx1m, Mugx1p;
  Vec Muox2m, Muox2p, Muwx2m, Muwx2p, Mugx2m, Mugx2p;
  Vec Muox3m, Muox3p, Muwx3m, Muwx3p, Mugx3m, Mugx3p;
  /* Density By Viscosity at the faces */
  Vec RhoByMuox1m, RhoByMuox1p, RhoByMuwx1m, RhoByMuwx1p, RhoByMugx1m, RhoByMugx1p;
  Vec RhoByMuox2m, RhoByMuox2p, RhoByMuwx2m, RhoByMuwx2p, RhoByMugx2m, RhoByMugx2p;
  Vec RhoByMuox3m, RhoByMuox3p, RhoByMuwx3m, RhoByMuwx3p, RhoByMugx3m, RhoByMugx3p;
  /* Relative Permeability on the faces */
  Vec RelPermox1m, RelPermox1p, RelPermox2m, RelPermox2p, RelPermox3m, RelPermox3p;
  Vec RelPermwx1m, RelPermwx1p, RelPermwx2m, RelPermwx2p, RelPermwx3m, RelPermwx3p;
  Vec RelPermgx1m, RelPermgx1p, RelPermgx2m, RelPermgx2p, RelPermgx3m, RelPermgx3p;

  /* Computational tools */
  Mat A;
  Vec RHS;
  KSP ksp;
  PC  pc;
  /* Gravity collected term */
  Vec Gravity;
  /* Capillary Pressure collected term for both IMPES and Newton Raphson */
  Vec CapillaryPressure;

  /* Computational tools for Newton Raphson */
  SNES snes;
  DA   NewtonRaphsonDA; /* conventions DOF[0] -> Po, DOF[1] -> Sw, DOF[2] -> So */
  Mat  NewtonRaphsonF;
  Vec  NewtonRaphsonSolution, NewtonRaphsonRHS;
  Vec  RHS0, RHS1, RHS2;    /* temporary right hand side vectors for Newton Raphson */
  Vec  Grav0, Grav1, Grav2; /* temporary gravity vectors  for Newton Raphson */
  Vec  NRCapPressure, NRGasCapPressure;
  Vec  NRSolution;
  Vec  PhiSoByBoOld, PhiSwByBwOld, PhiSgByBgPlusRsoSoByBoOld;

  PetscReal StartTime,  EndTime;
  PetscReal CurrentTimeP, DeltaTP;
  PetscReal CurrentTimeS, DeltaTS;
  PetscReal DSDTmax, DSmax;
  PetscTruth AdaptiveTimeStep;

  /* Parameters */
  PetscReal RockCompressibility;
  PetscReal RockCompRefPressure;
  PetscReal RockCompRefPorosity;
  /* Fluid compressibilities */
  PetscReal OilCompressibility;
  PetscReal WaterCompressibility;
  PetscReal GasCompressibility;
  /* Fluid reference pressures for compressibility */
  PetscReal OilCompPressure;
  PetscReal WaterCompPressure;
  PetscReal GasCompPressure;
  /* Molecular Weights of oil, water and gas*/
  PetscReal WO, WW, WG;
  /* Standard conditions densities */
  PetscReal Rhoos, Rhows, Rhogs;
  /* Rock/Fluid properties */
  PetscReal Sor, Swc, Sgc;

  /* Well information */
  PetscInt NumberOfWells;
  /* Pointer to wells */
  Well *Wells;
  /* Flow Rate Collected Term */
  PetscReal SimQo, SimQw, SimQg;
  PetscReal SimWOR, SimGOR;

  /* Gravitational Acceleration */
  PetscScalar GravAcc;

  /* IO control */
  PetscReal OutputEvery; /* output frequency, lazy evaluation with improved impes */

  /* Table relationships for capillary pressures */
  DefiantTable PcowSw;
  DefiantTable PcogSg;
  /* Table relationships for relative perms from sats */
  DefiantTable KrwSw;
  DefiantTable KrowSw;
  DefiantTable KrgSg;
  DefiantTable KrogSg;

} BlackOilReservoirSimulation;

/* Defiant Simulation Creation Routines */
extern PetscErrorCode DefiantCreateSimulationVecs(BlackOilReservoirSimulation* MySim);
extern PetscErrorCode DefiantCreate2PhIMPESSystem(BlackOilReservoirSimulation* MySim);
extern PetscErrorCode DefiantCreate3PhIMPESSystem(BlackOilReservoirSimulation* MySim);
extern PetscErrorCode DefiantCreate2PhNewtonRaphsonSystem(BlackOilReservoirSimulation* MySim);
extern PetscErrorCode DefiantCreate3PhNewtonRaphsonSystem(BlackOilReservoirSimulation* MySim);
/* Defiant Simulation Destruction Routines */
extern PetscErrorCode DefiantDestroySimulationVecs(BlackOilReservoirSimulation* MySim);
extern PetscErrorCode DefiantDestroy2PhIMPESSystem(BlackOilReservoirSimulation* MySim);
extern PetscErrorCode DefiantDestroy3PhIMPESSystem(BlackOilReservoirSimulation* MySim);
extern PetscErrorCode DefiantDestroy2PhNewtonRaphsonSystem(BlackOilReservoirSimulation* MySim);
extern PetscErrorCode DefiantDestroy3PhNewtonRaphsonSystem(BlackOilReservoirSimulation* MySim);

/* Rock update functions, when we have frac mechanics there will be more here */
extern PetscErrorCode DefiantUpdatePorosity(BlackOilReservoirSimulation* MySim);

/* PVT functions */
extern PetscErrorCode DefiantUpdateDensity(BlackOilReservoirSimulation* MySim);
extern PetscErrorCode DefiantUpdateVolumeFactors(BlackOilReservoirSimulation* MySim);
extern PetscErrorCode DefiantUpdateViscosities(BlackOilReservoirSimulation* MySim);
extern PetscErrorCode DefiantUpdatePcow(BlackOilReservoirSimulation* MySim);
extern PetscErrorCode DefiantUpdatePcog(BlackOilReservoirSimulation* MySim);
extern PetscErrorCode DefiantUpdatePwFromPcow(BlackOilReservoirSimulation* MySim);
extern PetscErrorCode DefiantUpdatePgFromPcog(BlackOilReservoirSimulation* MySim);

/* Relative Permeability routines */
extern PetscErrorCode DefiantUpdateRelativePermeability(BlackOilReservoirSimulation* MySim);
extern PetscErrorCode DefiantComputeRelativePermsFromSatsTables(BlackOilReservoirSimulation* MySim);
extern PetscErrorCode DefiantComputeRelativePermsFromSatsCorey(BlackOilReservoirSimulation* MySim);
extern PetscErrorCode DefiantComputeRelativePermsFromSatsNaar(BlackOilReservoirSimulation* MySim);
extern PetscErrorCode DefiantComputeRelativePermsFromSatsWygal(BlackOilReservoirSimulation* MySim);
extern PetscErrorCode DefiantComputeRelativePermsFromSatsStonesI(BlackOilReservoirSimulation* MySim);
extern PetscErrorCode DefiantComputeRelativePermsFromSatsStonesII(BlackOilReservoirSimulation* MySim);

/* Interpolation routines */
extern PetscErrorCode DefiantComputeKAByHAtFaces(BlackOilReservoirSimulation* MySim);
extern PetscErrorCode DefiantComputeRhoAndMuAtFaces(BlackOilReservoirSimulation* MySim);
extern PetscErrorCode DefiantComputeRelativePermsAtFaces(BlackOilReservoirSimulation* MySim);
extern PetscErrorCode DefiantComputeVolumeFactorsAtFaces(BlackOilReservoirSimulation* MySim);
extern PetscErrorCode DefiantComputeRsoAtFaces(BlackOilReservoirSimulation* MySim);

/* Transmissibility calculation routines */
extern PetscErrorCode DefiantComputeTransmissibilities(BlackOilReservoirSimulation* MySim);



/* General 2Ph functions for various solvers */
extern PetscErrorCode DefiantPhGravity(BlackOilReservoirSimulation* MySim);
extern PetscErrorCode DefiantPhCapillaryPressure(BlackOilReservoirSimulation* MySim);
/* IMPES 2Ph functions */
extern PetscErrorCode DefiantIMPES2PhAssembleMatrix(BlackOilReservoirSimulation* MySim);
extern PetscErrorCode DefiantIMPES2PhAssembleRHS(BlackOilReservoirSimulation* MySim);
extern PetscErrorCode DefiantIMPES2PhUpdateSaturations(BlackOilReservoirSimulation* MySim);
extern PetscErrorCode DefiantIMPES2PhHandleWells(BlackOilReservoirSimulation* MySim);
extern PetscErrorCode DefiantIMPES2PhHandleWellsForSaturation(BlackOilReservoirSimulation* MySim);
extern PetscErrorCode DefiantIMPES2PhUpdateSaturations(BlackOilReservoirSimulation* MySim);
extern PetscErrorCode DefiantClearMatrixRHS(BlackOilReservoirSimulation* MySim);
extern PetscErrorCode DefiantIMPES2PhSolve(BlackOilReservoirSimulation* MySim);
extern PetscErrorCode DefiantIMPES2PhIterate(BlackOilReservoirSimulation* MySim);
/* 2Ph DMMG functions */
extern PetscErrorCode DefiantIMPES2PhDMMGComputeRHS(DMMG dmmg, Vec b);
extern PetscErrorCode DefiantIMPES2PhDMMGComputeMatrix(DMMG dmmg, Mat jac, Mat B);
extern PetscErrorCode DefiantIMPES2PhDMMGSolve(BlackOilReservoirSimulation* MySim);
extern PetscErrorCode DefiantIMPES2PhDMMGIterate(BlackOilReservoirSimulation* MySim);

/* NewtonRaphson 2Ph functions */
extern PetscErrorCode DefiantNewtonRaphson2PhGravity(BlackOilReservoirSimulation* MySim);
extern PetscErrorCode DefiantNewtonRaphson2PhCapillaryPressure(BlackOilReservoirSimulation* MySim);
extern PetscErrorCode DefiantNewtonRaphson2PhHandleWells(BlackOilReservoirSimulation* MySim);
extern PetscErrorCode DefiantNewtonRaphson2PhFormFunction(SNES snes,Vec X,Vec F,void *ptr);
extern PetscErrorCode DefiantNewtonRaphson2PhSolve(BlackOilReservoirSimulation* MySim);
extern PetscErrorCode DefiantNewtonRaphson2PhComputeLegacyTerms(BlackOilReservoirSimulation* MySim);
extern PetscErrorCode DefiantNewtonRaphson2PhIterate(BlackOilReservoirSimulation* MySim);

/* IMPES 3Ph functions */
extern PetscErrorCode DefiantIMPES3PhGravity(BlackOilReservoirSimulation* MySim);
extern PetscErrorCode DefiantIMPES3PhCapillaryPressure(BlackOilReservoirSimulation* MySim);
extern PetscErrorCode DefiantIMPES3PhUpdateSaturations(BlackOilReservoirSimulation* MySim);
extern PetscErrorCode DefiantIMPES3PhHandleWells(BlackOilReservoirSimulation* MySim);
extern PetscErrorCode DefiantIMPES3PhHandleWellsForSaturation(BlackOilReservoirSimulation* MySim);
extern PetscErrorCode DefiantIMPES3PhAssembleMatrix(BlackOilReservoirSimulation* MySim);
extern PetscErrorCode DefiantIMPES3PhAssembleRHS(BlackOilReservoirSimulation* MySim);
extern PetscErrorCode DefiantIMPES3PhSolve(BlackOilReservoirSimulation* MySim);
extern PetscErrorCode DefiantIMPES3PhUpdateSaturations(BlackOilReservoirSimulation* MySim);

extern PetscErrorCode DefiantIMPES3PhDMMGHandleWells(BlackOilReservoirSimulation* MySim);
extern PetscErrorCode DefiantIMPES3PhDMMGComputeRHS(DMMG dmmg, Vec b);
extern PetscErrorCode DefiantIMPES3PhDMMGComputeMatrix(DMMG dmmg, Mat jac, Mat B);
extern PetscErrorCode DefiantIMPES3PhDMMGSolve(BlackOilReservoirSimulation* MySim);
extern PetscErrorCode DefiantIMPES3PhDMMGIterate(BlackOilReservoirSimulation* MySim);

/* NewtonRaphson 3Ph functions */
extern PetscErrorCode DefiantNewtonRaphson3PhGravity(BlackOilReservoirSimulation* MySim);
extern PetscErrorCode DefiantNewtonRaphson3PhCapillaryPressure(BlackOilReservoirSimulation* MySim);
extern PetscErrorCode DefiantNewtonRaphson3PhHandleWells(BlackOilReservoirSimulation* MySim);
extern PetscErrorCode DefiantNewtonRaphson3PhFormFunction(SNES snes,Vec X,Vec F,void *ptr);
extern PetscErrorCode DefiantNewtonRaphson3PhSolve(BlackOilReservoirSimulation* MySim);
extern PetscErrorCode DefiantNewtonRaphson3PhComputeLegacyTerms(BlackOilReservoirSimulation* MySim);
extern PetscErrorCode DefiantNewtonRaphson3PhIterate(BlackOilReservoirSimulation* MySim);

/* Clean matrix and RHS */
extern PetscErrorCode DefiantClearMatrixRHS(BlackOilReservoirSimulation* MySim);

/* Relative permeability from saturation routines */
extern PetscErrorCode DefiantComputeKrw(PetscScalar  * Krw,  PetscScalar Sw);
extern PetscErrorCode DefiantComputeKrg(PetscScalar  * Krg,  PetscScalar Sg);
extern PetscErrorCode DefiantComputeKrow(PetscScalar * Krow, PetscScalar Sw);
extern PetscErrorCode DefiantComputeKrog(PetscScalar * Krog, PetscScalar Sg);

/* Production functions for BlackOil */
extern PetscErrorCode DefiantBlackOilSyncPerfOwners(BlackOilReservoirSimulation* MySim);
extern PetscErrorCode DefiantBlackOilComputePerfIndicesSyncPerfs(BlackOilReservoirSimulation* MySim);
extern PetscErrorCode DefiantBlackOil2PhProduction(BlackOilReservoirSimulation* MySim);
extern PetscErrorCode DefiantBlackOil3PhProduction(BlackOilReservoirSimulation* MySim);

/* Well handling and book-keeping functions */
extern PetscErrorCode DefiantAddPerforation(Well *MyWell);
extern PetscErrorCode DefiantAddWell(BlackOilReservoirSimulation* MySim);

/* These are global Statics that we sadly have to have for call-backs */
static BlackOilReservoirSimulation* CurrentSimulation;

/* Some examples */
extern PetscErrorCode DefiantIMPES2PhBuckleyLeverett();
extern PetscErrorCode DefiantIMPES3PhBuckleyLeverett();
extern PetscErrorCode DefiantNewton2PhBuckleyLeverett();
extern PetscErrorCode DefiantNewton3PhBuckleyLeverett();
extern PetscErrorCode DefiantIMPES2PhFivePoint();
extern PetscErrorCode DefiantIMPES3PhFivePoint();
extern PetscErrorCode DefiantNewton2PhFivePoint();
extern PetscErrorCode DefiantNewton3PhFivePoint();
extern PetscErrorCode DefiantIMPES2PhBenchmark();

/* Geometry routines */
extern PetscErrorCode DefiantSetDACoords(BlackOilReservoirSimulation* MySim);
extern PetscErrorCode DefiantGetDACoords(BlackOilReservoirSimulation* MySim);
extern PetscErrorCode DefiantSetDACoords(BlackOilReservoirSimulation* MySim);
extern PetscErrorCode DefiantViewCoords(BlackOilReservoirSimulation* MySim);

/* various viewing routines for pressure */
extern PetscErrorCode DefiantViewPressuresSTDOUT(BlackOilReservoirSimulation* MySim);
extern PetscErrorCode DefiantViewPressuresASCII(BlackOilReservoirSimulation* MySim);

/* various viewing routines for saturation */
extern PetscErrorCode DefiantViewSaturationsSTDOUT(BlackOilReservoirSimulation* MySim);
extern PetscErrorCode DefiantViewSaturationsASCII(BlackOilReservoirSimulation* MySim);


/* These routines are meant for reading and writing out variables to Valiant */
extern PetscErrorCode DefiantBlackOil2PhValiantLoadVecs(BlackOilReservoirSimulation* MySim);
extern PetscErrorCode DefiantBlackOil3PhValiantLoadVecs(BlackOilReservoirSimulation* MySim);
extern PetscErrorCode DefiantBlackOil2PhValiantWriteVecs(BlackOilReservoirSimulation* MySim);
extern PetscErrorCode DefiantBlackOil3PhValiantWriteVecs(BlackOilReservoirSimulation* MySim);

#endif /* DEFIANT_H_ */
