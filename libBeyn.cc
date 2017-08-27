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

/***************************************************************/
/***************************************************************/
/***************************************************************/
cdouble zRandom()
{ 
  return 2.0*(drand48()-0.5) + II*2.0*(drand48()-0.5);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
BeynData *CreateBeynData(int M, int L)
{
  BeynData *Data = (BeynData *)mallocEC(sizeof(*Data));

  Data->M = M;
  Data->L = L;
  Data->ShiftEVs = true;

  int MLMax = (M>L) ? M : L;
  int MLMin = (M<L) ? M : L;

  size_t Size1 = MLMax*L*sizeof(cdouble);
  size_t Size2 = L*L*sizeof(cdouble);
  Data->Workspaces[0] = mallocEC(Size1);
  Data->Workspaces[1] = mallocEC(Size1);
  Data->Workspaces[2] = mallocEC(Size1);
  Data->Workspaces[3] = mallocEC(Size2);

  srand48(time(0));
  Data->VHat = new HMatrix(M,L,LHM_COMPLEX);
  for(int m=0; m<M; m++)
   for(int l=0; l<L; l++)
    Data->VHat->SetEntry(m, l, zRandom() );

  Data->Sigma  = new HVector(MLMin);
  Data->Lambda = new HVector(L, LHM_COMPLEX);
  Data->Eigenvectors = new HMatrix(M, L, LHM_COMPLEX);

Data->V0Full=new HMatrix(M, M, LHM_COMPLEX); // FIXME
Data->W0TFull=new HMatrix(L, L, LHM_COMPLEX); // FIXME

  return Data;
  
}

void DestroyBeynData(BeynData *Data)
{
  delete Data->VHat;
  for(int ntm=0; ntm<4; ntm++)
   free(Data->Workspaces[ntm]);
  delete Data->Sigma;
  delete Data->Lambda;

delete Data->V0Full; // FIXME
delete Data->W0TFull; // FIXME

  delete Data;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
int BeynMethod(BeynData *Data, cdouble z0, double R,
               BeynFunction UserFunction, void *UserData, int N)
{  

  Log("Applying Beyn method with z0=%s,R=%e,N=%i...",z2s(z0),R,N);

  int M                 = Data->M;
  int L                 = Data->L;
  bool ShiftEVs         = Data->ShiftEVs;
  HMatrix *VHat         = Data->VHat;  
  void **Workspaces     = Data->Workspaces;
  HVector *Sigma        = Data->Sigma;
  HVector *Lambda       = Data->Lambda;

HMatrix *V0Full = Data->V0Full;   // FIXME
HMatrix *W0TFull = Data->W0TFull; // FIXME

  /***************************************************************/
  /* evaluate contour integrals by numerical quadrature to get   */
  /* A0 and A1 matrices                                          */
  /***************************************************************/
  HMatrix TM(M,L,LHM_COMPLEX,Workspaces[0]);
  HMatrix A0(M,L,LHM_COMPLEX,Workspaces[1]);
  HMatrix A1(M,L,LHM_COMPLEX,Workspaces[2]);
  A0.Zero();
  A1.Zero();
  double DeltaTheta = 2.0*M_PI / ((double)N);
  Log(" Evaluating contour integral (%i points)...",N);
  for(int n=0; n<N; n++)
   { 
     cdouble Xi = exp(II*((double)n)*DeltaTheta);
     cdouble dz = R*Xi;
     cdouble z  = z0 + dz;
     TM.Copy(VHat);
     UserFunction(z, UserData, &TM);
     VecPlusEquals(A0.ZM, Xi, TM.ZM, M*L);
     VecPlusEquals(A1.ZM, Xi * (ShiftEVs ? dz : z), TM.ZM, M*L);
   };
  A0.Scale( R / ((double)N) );
  A1.Scale( R / ((double)N) );

  /***************************************************************/
  /* perform linear algebra manipulations to get eigenvalues and */
  /* eigenvectors                                                */
  /***************************************************************/
#if 0
V0Full: MxM
V0:     MxK
W0T:    KxK
TM:     KxK
B:      KxK
S:      KxK
Eigenvectors: MxK
#endif

  // A0 -> V0 * Sigma * W0'
  Log(" Computing SVD...");
  A0.SVD(Sigma, V0Full, W0TFull);

  // compute effective rank Q
  int K=0;
  double SigmaThreshold=1.0e-8;
  for(int k=0; k<Sigma->N; k++)
   if (Sigma->GetEntryD(k) > SigmaThreshold)
    K++;
  Log(" %i/%i relevant singular values",K,L);
  if (K==L)
   Warn("K=L=%i in BeynMethod (repeat with higher L)",K,L);

  HMatrix V0(M,K,LHM_COMPLEX,LHM_NORMAL,Workspaces[1]);
  V0Full->ExtractBlock(0, 0, &V0);

  HMatrix W0T(K,L,LHM_COMPLEX,LHM_NORMAL,Workspaces[0]);
  W0TFull->ExtractBlock(0,0, &W0T);

  // B <- V0' * A1 * W0 * Sigma^-1
  HMatrix TM2(K,L,LHM_COMPLEX,LHM_NORMAL,Workspaces[3]);
  HMatrix B(K,K,LHM_COMPLEX,LHM_NORMAL,Workspaces[2]);
  Log(" Multiplying V0*A1->TM...");
  V0.Multiply(&A1, &TM2, "--transA C");   // TM <- V0' * A1
  Log(" Multiplying TM*W0T->B...");
  TM2.Multiply(&W0T, &B, "--transB C");   //  B <- TM * W0
  for(int m=0; m<B.NR; m++)             //  B <- B * Sigma^{-1}
   for(int n=0; n<B.NC; n++)
    B.SetEntry(m,n,B.GetEntry(m,n)/Sigma->GetEntry(n));

  HMatrix S(K,K,LHM_COMPLEX,LHM_NORMAL,Workspaces[0]);
  HMatrix Eigenvectors(M,K,LHM_COMPLEX,LHM_NORMAL,Workspaces[2]);
  Lambda->Zero();
  HVector LambdaK(K,LHM_COMPLEX,((void *)(Lambda->ZV)));
  Log(" Eigensolving...");
  B.NSEig(&LambdaK, &S);
  if (ShiftEVs)
   for(int k=0; k<K; k++)
    Lambda->ZV[k] += z0;
  Log(" Multiplying V0*S...");
  V0.Multiply(&S, &Eigenvectors);

  Data->Eigenvectors->InsertBlock(&Eigenvectors, 0, 0);
  return K;
}
