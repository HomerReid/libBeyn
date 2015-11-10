#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <libhrutil.h>
#include <libscuff.h>
#include <libscuffInternals.h>
#include <libSGJC.h>
#include <scuff-spectrum.h>

#include <nlopt.h>

#define II cdouble(0.0,1.0)

using namespace scuff;

/***************************************************************/
/***************************************************************/
/***************************************************************/
void AssembleM(SSData *Data, cdouble Omega, double *kBloch, HMatrix *M)
{
  RWGGeometry *G = Data->G;
  int NS         = G->NumSurfaces;

  static void **Accelerators=0;
  if (!Accelerators)
   { Accelerators=(void **)mallocEC(NS*NS*sizeof(void *));
     for(int nsa=0; nsa<NS; nsa++)
      for(int nsb=0; nsb<NS; nsb++)
       Accelerators[nsa*NS + nsb]=G->CreateABMBAccelerator(nsa, nsb);
   };
 
  for(int nsa=0; nsa<NS; nsa++)
   for(int nsb=0; nsb<NS; nsb++)
    { int RowOffset=G->BFIndexOffset[nsa];
      int ColOffset=G->BFIndexOffset[nsb];
      G->AssembleBEMMatrixBlock(nsa, nsb, Omega, kBloch,
                                M, 0, RowOffset, ColOffset,
                                Accelerators[nsa*NS + nsb]);
    };
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
double GetRCond(HMatrix *M, cdouble *LogDet)
{  
  double Norm = M->GetNorm();
  M->LUFactorize();

  cdouble Sum=0.0;
  for (int n=0; n<M->NR; n++)
   Sum += log(M->GetEntry(n,n));
  *LogDet = Sum;

  return M->GetRCond( Norm );
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
double Objective(unsigned n, const double *x,
                 double* grad, void* UserData)
{
   SSData *Data = (SSData* )UserData;

   RWGGeometry *G     = Data->G;
   HMatrix *M         = Data->M;
   HMatrix *Mp        = Data->Mp;
   HMatrix *Mm        = Data->Mm;
   HVector *Lambda    = Data->Lambda;
   char *FileBase     = Data->FileBase;
   double Delta       = Data->Delta;
   bool SecondOrderFD = Data->SecondOrderFD;
   bool Eigenvalues   = Data->Eigenvalues;
   double *gradFFD    = Data->gradFFD;
   double *gradCFD    = Data->gradCFD;

   /***************************************************************/
   /***************************************************************/
   /***************************************************************/
   cdouble Omega;
   double *kBloch, kBlochBuffer[2];
   if (n==2)
    { Omega = cdouble(x[0],x[1]);
      kBloch  = Data->kBloch;
    }
   else if (n==1)
    { Omega = Data->Omega;
      kBloch=kBlochBuffer;
      kBloch[0]=x[0];
      kBloch[1]=0.0;
    }
   else
    { kBloch = Data->kBloch;
      Omega  = Data->Omega; 
    };

   /***************************************************************/
   /***************************************************************/
   /***************************************************************/
   Log("Assembling M at {w,k}={%e,%e} {%e,%e}...",
        real(Omega),imag(Omega),kBloch[0],kBloch[1]);
   AssembleM(Data, Omega, kBloch, M);

   /***************************************************************/
   /***************************************************************/
   /***************************************************************/
   cdouble LMinAbs=0.0, LMinRe=0.0, LMinIm=0.0;
   if (Eigenvalues)
    { Mp->Copy(M);
      Mp->NSEig(Lambda);
      double MinAbs=1.0e89, MinRe=1.0e89, MinIm=1.0e89;
      for(n=0; n<Lambda->N; n++)
       { cdouble L=Lambda->GetEntry(n);
         if (abs(L) < MinAbs) 
          { MinAbs=abs(L);      LMinAbs=L; };
         if (abs(real(L)) < MinRe) 
          { MinRe=abs(real(L)); LMinRe=L;  };
         if (imag(L) < MinIm) 
          { MinIm=abs(imag(L)); LMinIm=L;  };
       };
    };

   cdouble LogDetM;
   double MRCond = Data->RCond = GetRCond(M, &LogDetM);

   /***************************************************************/
   /***************************************************************/
   /***************************************************************/
   if (grad)
    for(int Mu=0; Mu<2; Mu++)
     { 
       double DeltaOmega = Delta*fabs(x[Mu]);
       if (DeltaOmega==0.0) DeltaOmega=Delta;

       cdouble dOmega = (Mu==0) ? DeltaOmega : II*DeltaOmega;

       AssembleM(Data, Omega + dOmega, kBloch, Mp);
       cdouble LogDetP;
       double MpRCond=GetRCond(Mp,&LogDetP);
       grad[Mu]=gradFFD[Mu] = ( MpRCond-MRCond ) / fabs(DeltaOmega);
       if (SecondOrderFD)
        { AssembleM(Data, Omega - dOmega, kBloch, Mm);
          cdouble LogDetm;
          double MmRCond=GetRCond(Mm,&LogDetm);
          gradCFD[Mu] = ( MpRCond-MmRCond ) / (2.0*fabs(DeltaOmega));
          grad[Mu]=gradCFD[Mu];
        };
     };

   /*--------------------------------------------------------------*/
   /*--------------------------------------------------------------*/
   /*--------------------------------------------------------------*/
   char FileName[100];
   snprintf(FileName, 100, "%s.spectrum",FileBase);
   static bool WroteHeader=false;
   if (!WroteHeader)
    { WroteHeader=true;
      FILE *ff=fopen(FileName,"a");
      fprintf(ff,"\n\n");
      fprintf(ff,"# scuff-spectrum running on %s at %s\n",GetHostName(),GetTimeString());
      fprintf(ff,"# 1   2 re, im Omega\n");
      fprintf(ff,"# 3   4 kx, ky \n");
      fprintf(ff,"# 5   6 re, im log det M\n");
      fprintf(ff,"# 7     MRCond\n");
      int nc=8;
      if (grad)
       { fprintf(ff,"# %i %i grad f (FFD)\n",nc,nc+1);
         fprintf(ff,"# %i %i grad f (CFD)\n",nc+2,nc+3);
         nc+=4;
       };
      if (Eigenvalues)
       { fprintf(ff,"# %i %i arg min abs|eigenvalue| \n",nc,nc+1);
         nc+=2;
         fprintf(ff,"# %i %i arg min re|eigenvalue| \n",nc,nc+1);
         nc+=2;
         fprintf(ff,"# %i %i arg min im|eigenvalue| \n",nc,nc+1);
         nc+=2;
       };
      fclose(ff);
    };

   /*--------------------------------------------------------------*/
   /*--------------------------------------------------------------*/
   /*--------------------------------------------------------------*/
   FILE *ff=fopen(FileName,"a");
   fprintf(ff,"%+.10e %+.10e ",real(Omega),imag(Omega));
   fprintf(ff,"%e %e ",kBloch[0],kBloch[1]);
   fprintf(ff,"%e %e ",real(LogDetM),imag(LogDetM));
   fprintf(ff,"%e    ",MRCond);
   if (grad)
    { fprintf(ff,"%.2e %.2e ",gradFFD[0],gradFFD[1]);
      fprintf(ff,"%.2e %.2e ",gradCFD[0],gradCFD[1]);
    };
   if (Eigenvalues)
    { fprintf(ff,"%+.8e %.8e \n",real(LMinAbs),imag(LMinAbs));
      fprintf(ff,"%+.8e %.8e \n",real(LMinRe),imag(LMinRe));
      fprintf(ff,"%+.8e %.8e \n",real(LMinIm),imag(LMinIm));
    };
   fprintf(ff,"\n");
   fclose(ff);

  return MRCond;
}
