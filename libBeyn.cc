/* Copyright (C) 2005-2011 M. T. Homer Reid
 *
 * This file is part of SCUFF-EM.
 *
 * SCUFF-EM is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * SCUFF-EM is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/*
 * libBeyn.cc      -- implementation of Beyn's method for
 *                 -- nonlinear eigenvalue problems
 *
 * Homer Reid      -- 6/2016
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <libhrutil.h>
#include <libhmat.h>

#include <libBeyn.h>

#define II cdouble(0.0,1.0)

cdouble zrandN(double Mu=0.0, double Sigma=1.0)
{ return cdouble(randN(Mu,Sigma), randN(Mu,Sigma)); }

/***************************************************************/
/***************************************************************/
/***************************************************************/
BeynSolver *CreateBeynSolver(int M, int L)
{
  BeynSolver *Solver= (BeynSolver *)mallocEC(sizeof(*Solver));

  Solver->M = M;
  Solver->L = L;

  int MLMax = (M>L) ? M : L;
  int MLMin = (M<L) ? M : L;

  // storage for singular values, eigenvalues, and eigenvectors
  Solver->Sigma        = new HVector(MLMin);
  Solver->Lambda       = new HVector(L, LHM_COMPLEX);
  Solver->Eigenvectors = new HMatrix(M, L, LHM_COMPLEX);

  // random matrix used in algorithm
  Solver->VHat         = new HMatrix(M,L,LHM_COMPLEX);
  ReRandomize(Solver);

  // internal workspace: need storage for 3 MxL matrices
  // plus 3 LxL matrices
  #define MLBUFFERS 3
  #define LLBUFFERS 3
  size_t ML = MLMax*L, LL = L*L;
  size_t WorkspaceSize = (MLBUFFERS*ML + LLBUFFERS*LL)*sizeof(cdouble);
 
  Solver->Workspace = (cdouble *)mallocEC( WorkspaceSize );

  return Solver;
  
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void DestroyBeynSolver(BeynSolver *Solver)
{
  delete   Solver->Sigma;
  delete   Solver->Lambda;
  delete   Solver->Eigenvectors;
  delete   Solver->VHat;
  free(Solver->Workspace);

  delete Solver;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void ReRandomize(BeynSolver *Solver, unsigned int RandSeed)
{
  if (RandSeed==0) 
   RandSeed=time(0);
  srandom(RandSeed);
  HMatrix *VHat=Solver->VHat;
  for(int nr=0; nr<VHat->NR; nr++)
   for(int nc=0; nc<VHat->NC; nc++)
    VHat->SetEntry(nr,nc,zrandN());

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
int BeynSolve(BeynSolver *Solver,
              BeynFunction UserFunction, void *UserData,
              cdouble z0, double Rx, double Ry, int N)
{  
  if (Rx==Ry)
   Log("Applying Beyn method with z0=%s,R=%e,N=%i...",z2s(z0),Rx,N);
  else
   Log("Applying Beyn method with z0=%s,Rx=%e,Ry=%e,N=%i...",z2s(z0),Rx,Ry,N);

  int M                 = Solver->M;
  int L                 = Solver->L;
  HMatrix *VHat         = Solver->VHat;
  HVector *Sigma        = Solver->Sigma;
  HVector *Lambda       = Solver->Lambda;
  HMatrix *Eigenvectors = Solver->Eigenvectors;
  cdouble *Workspace    = Solver->Workspace;

  int MLMax = M>L ? M : L;
  size_t ML = MLMax*L, LL = L*L;
  cdouble *MLBuffers[3], *LLBuffers[3];
  MLBuffers[0] = Workspace;
  MLBuffers[1] = MLBuffers[0] + ML;
  MLBuffers[2] = MLBuffers[1] + ML;
  LLBuffers[0] = MLBuffers[2] + ML;
  LLBuffers[1] = LLBuffers[0] + LL;
  LLBuffers[2] = LLBuffers[1] + LL;

  /***************************************************************/
  /* evaluate contour integrals by numerical quadrature to get   */
  /* A0 and A1 matrices                                          */
  /***************************************************************/
  HMatrix       A0(M,L,LHM_COMPLEX,(void *)MLBuffers[0]);
  HMatrix       A1(M,L,LHM_COMPLEX,(void *)MLBuffers[1]);
  HMatrix MInvVHat(M,L,LHM_COMPLEX,(void *)MLBuffers[2]);
  A0.Zero();
  A1.Zero();
  double DeltaTheta = 2.0*M_PI / ((double)N);
  Log(" Evaluating contour integral (%i points)...",N);
  for(int n=0; n<N; n++)
   { 
     double Theta = ((double)n)*DeltaTheta;
     double CT    = cos(Theta), ST=sin(Theta);
     cdouble z1   = Rx*CT + II*Ry*ST;
     cdouble dz   = (II*Rx*ST + Ry*CT)/((double)N);
     MInvVHat.Copy(VHat);
     UserFunction(z0+z1, UserData, &MInvVHat);
     VecPlusEquals(A0.ZM, dz,    MInvVHat.ZM, M*L);
     VecPlusEquals(A1.ZM, z1*dz, MInvVHat.ZM, M*L);
   };

  /***************************************************************/
  /* perform linear algebra manipulations to get eigenvalues and */
  /* eigenvectors                                                */
  /***************************************************************/
  // A0 -> V0Full * Sigma * W0TFull'
  Log(" Computing SVD...");
  HMatrix V0Full(M,L,LHM_COMPLEX,(void *)MLBuffers[2]);
  HMatrix W0TFull(L,L,LHM_COMPLEX,(void *)LLBuffers[0]);
  A0.SVD(Sigma, &V0Full, &W0TFull);

  // compute effective rank K
  int K=0;
  double SigmaThreshold=1.0e-8;
  for(int k=0; k<Sigma->N; k++)
   if (Sigma->GetEntryD(k) > SigmaThreshold)
    K++;
  Log(" %i/%i relevant singular values",K,L);
  if (K==L)
   Warn("K=L=%i in BeynMethod (repeat with higher L)",K,L);
  if (K==0)
   { Warn("no eigenvalues found!");
     return 0;
   };

  // set V0, W0T = matrices of first K right, left singular vectors
  HMatrix V0(M,K,LHM_COMPLEX,(void *)MLBuffers[0]);
  HMatrix W0T(K,L,LHM_COMPLEX,(void *)LLBuffers[1]);
  for(int k=0; k<K; k++)
   { for(int m=0; m<M; m++) V0.SetEntry(m,k,V0Full.GetEntry(m,k));
     for(int l=0; l<L; l++) W0T.SetEntry(k,l,W0TFull.GetEntry(k,l));
   };

  // B <- V0' * A1 * W0 * Sigma^-1
  HMatrix TM2(K,L,LHM_COMPLEX,(void *)LLBuffers[0]);
  HMatrix B(K,K,LHM_COMPLEX,(void *)LLBuffers[2]);
  Log(" Multiplying V0*A1->TM...");
  V0.Multiply(&A1, &TM2, "--transA C");   // TM2 <- V0' * A1
  Log(" Multiplying TM*W0T->B...");
  TM2.Multiply(&W0T, &B, "--transB C");   //  B <- TM2 * W0
  for(int m=0; m<K; m++)                  //  B <- B * Sigma^{-1}
   for(int n=0; n<K; n++)
    B.ScaleEntry(m,n,1.0/Sigma->GetEntry(n));

  Log(" Eigensolving...");
  HVector LambdaK(K,LHM_COMPLEX,(void *)LLBuffers[0]);
  HMatrix S(K,K,LHM_COMPLEX,(void *)LLBuffers[1]);
  B.NSEig(&LambdaK, &S);

  Log(" Multiplying V0*S...");
  HMatrix EigenvectorsK(M,K,LHM_COMPLEX,(void *)MLBuffers[1]);
  V0.Multiply(&S, &EigenvectorsK);

  Solver->Lambda->Zero();
  Solver->Eigenvectors->Zero();
  for(int k=0; k<K; k++)
   { Solver->Lambda->SetEntry(k, LambdaK.GetEntry(k) + z0 );
     for(int m=0; m<M; m++)
      Solver->Eigenvectors->SetEntry(m,k, EigenvectorsK.GetEntry(m,k));
   };

  return K;
}

int BeynSolve(BeynSolver *Solver,
              BeynFunction UserFunction, void *UserData,
              cdouble z0, double R, int N)
{ return BeynSolve(Solver, UserFunction, UserData, z0, R, R, N); }
