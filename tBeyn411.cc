#include <stdio.h>
#include <stdlib.h>

#include <string.h> 
#include <math.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>

#define HAVE_READLINE
#include <readline/readline.h>
#include <readline/history.h>

#include <libhrutil.h>
#include <libscuff.h>
#include <libSGJC.h>
#include "libBeyn.h"

#define II cdouble(0.0,1.0)

using namespace scuff;

/***************************************************************/
/* function passed to BeynSolver, called for each quadrature   */
/* point on the contour; its task is to replace the MxL matrix */
/* VHat with M(z) \ VHat, i.e. Inverse[ M(z) ] * VHat          */
/***************************************************************/
typedef struct BFData
 { 
   HMatrix *M;
 } BFData;

void BeynFunc(cdouble z, void *UserData, HMatrix *VHat)
{
  // unpack fields from data structure
  BFData *Data   = (BFData *)UserData;
  HMatrix *M     = Data->M;

  // assemble M matrix as described in Beyn section 4.11
  int NR = M->NR;
  double m = (double)NR;
  cdouble DE  =  2.0*m - 4.0*z/(6.0*m);   // diagonal entry
  cdouble ODE = -1.0*m -     z/(6.0*m);   // off-diagonal entry
  M->Zero();
  M->SetEntry(0,0,DE);  M->SetEntry(0,1,ODE); // top row
  for(int nr=1; nr<NR-1; nr++)                // rows 2...m-1
   { M->SetEntry(nr,nr,DE);
     M->SetEntry(nr,nr-1,ODE);
     M->SetEntry(nr,nr+1,ODE);
   };
  M->SetEntry(NR-1, NR-2, ODE);
  M->SetEntry(NR-1, NR-1, 0.5*DE + z/(z-1.0) ); // bottom row

  // replace VHat with M\VHat
  M->LUFactorize();
  M->LUSolve(VHat);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[]) 
{ 
  InstallHRSignalHandler();
  InitializeLog(argv[0]);

  /***************************************************************/
  /* process command-line arguments ******************************/
  /***************************************************************/
  cdouble z0  =  cdouble(150.0, 0.0);
  double Rx   = 148.0;
  double Ry   = 148.0;
  int L       = 10;
  int N       = 50;
  int Dim     = 400;
  /* name        type    #args  max_instances  storage    count  description*/
  OptStruct OSArray[]=
   { {"z0",                 PA_CDOUBLE, 1, 1, (void *)&z0,   0,   "center of elliptical contour ({150,0})"},
//
     {"Rx",                 PA_DOUBLE,  1, 1, (void *)&Rx,   0,   "horizontal radius of elliptical contour (148)"},
//
     {"Ry",                 PA_DOUBLE,  1, 1, (void *)&Ry,   0,   "vertical radius of elliptical contour (148)"},
//
     {"N",                  PA_INT,     1, 1, (void *)&N,    0,   "number of quadrature points (50)"},
//
     {"L",                  PA_INT,     1, 1, (void *)&L,    0,   "number of EVs expected in contour (10)"},
//
     {"Dim",                PA_INT,     1, 1, (void *)&Dim,  0,   "dimension of problem (400)"},
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);
  SetDefaultCD2SFormat("{%+.8e,%+.8e}");

  /***************************************************************/
  /* initialize BeynSolver data structure ************************/
  /***************************************************************/
  BeynSolver *Solver = CreateBeynSolver(Dim, L);

  /***************************************************************/
  /* initialize data structure needed for our BeynFunction  ******/
  /***************************************************************/
  struct BFData MyBFData, *BFD=&MyBFData;
  BFD->M = new HMatrix(Dim, Dim, LHM_COMPLEX);

  /***************************************************************/
  /* compute eigenvalues by Beyn method                          */
  /***************************************************************/
  int K=BeynSolve(Solver, BeynFunc, (void *)BFD, z0, Rx, Ry, N);

  printf("Found %i eigenvalues: \n",K);
  for(int k=0; k<K; k++)
   printf("%i: %s \n",k,CD2S(Solver->Lambda->GetEntry(k)));

}
