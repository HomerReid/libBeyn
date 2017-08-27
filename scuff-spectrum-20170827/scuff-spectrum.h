#ifndef SCUFF_SPECTRUM_H
#define SCUFF_SPECTRUM_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <libhrutil.h>
#include <libscuff.h>
#include <nlopt.h>

#define II cdouble(0.0,1.0)

using namespace scuff;

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct SSData
 { 
   RWGGeometry *G;
   HMatrix *M, *Mp, *Mm;
   HVector *Lambda;
   HVector *Sigma;
   cdouble Omega;
   double *kBloch;
   char *FileBase;
   double Delta;
   double RCond;
   double gradFFD[2], gradCFD[2];
   bool EstimateDerivative;
   bool Console;
   bool SecondOrderFD;
   bool DoOptimization;
   bool DBOptimization;
   bool DoRCond;
   bool Eigenvalues; 
   bool SVD;

 } SSData;

double Objective(unsigned n, const double* x,
                 double* grad, void* UserData);

#endif // SCUFF_SPECTRUM_H
