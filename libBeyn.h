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
typedef struct BeynSolver
{
   int M;   // dimension of matrices
   int L;   // number of columns of VHat matrix

   HMatrix *VHat;
   HVector *Sigma, *Lambda;
   HMatrix *Eigenvectors;
   cdouble *Workspace;

 } BeynSolver;

// constructor, destructor
BeynSolver *CreateBeynSolver(int M, int L);
void DestroyBeynSolver(BeynSolver *Solver);

// reset the random matrix VHat used in the Beyn algorithm
// 
void ReRandomize(BeynSolver *Solver, unsigned int RandSeed=0);

// for both of the following routines,
// the return value is the number of eigenvalues found,
// and the eigenvalues and eigenvectors are stored in the
// Lambda and Eigenvectors fields of the BeynSolver structure

// Beyn method for circular contour of radius R,
// centered at z0, using N quadrature points
int BeynSolve(BeynSolver *Solver,
              BeynFunction UserFunction, void *UserData,
              cdouble z0, double R, int N=25);

// Beyn method for elliptical contour of horizontal, vertical
// radii Rx, Ry, centered at z0, using N quadrature points
int BeynSolve(BeynSolver *Solver,
              BeynFunction UserFunction, void *UserData,
              cdouble z0, double Rx, double Ry, int N=25);

#endif
