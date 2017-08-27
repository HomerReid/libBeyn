#ifndef BEYNMETHOD_H
#define BEYNMETHOD_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>
#include <fenv.h>

#include <libhmat.h>

/***************************************************************/
/* prototype for user-supplied function passed to BeynMethod.  */
/* The user's function should replace VHat with                */
/*  Inverse[ M(z) ] * VHat.                                    */
/***************************************************************/
typedef void (*BeynFunction)(cdouble z, void *UserData, HMatrix *VHat);

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct BeynData
{
   int M;   // dimension of matrices
   int L;   // number of columns of VHat matrix
   bool ShiftEVs;

   HMatrix *VHat;
   HVector *Sigma, *Lambda;
   void *Workspaces[4];
   
HMatrix *V0Full, *W0TFull; // FIXME
HMatrix *Eigenvectors;

 } BeynData;

BeynData *CreateBeynData(int M, int L);
void DestroyBeynData(BeynData *Data);
int BeynMethod(BeynData *Data, cdouble z0, double R,
               BeynFunction UserFunction, void *UserData, int N=25);

#endif
