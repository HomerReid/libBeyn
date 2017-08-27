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
#include <libscuffInternals.h>
#include <libSGJC.h>

#define II cdouble(0.0,1.0)

using namespace scuff;

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct MyData
 { 
   RWGGeometry *G;
   HMatrix *M, *Mp, *Mm;
   HVector *Lambda;
   char FileName[100];
   double Delta;
   double *kBloch;
   double Phi;
   bool EstimateDerivative;
   bool Console;

 } MyData;

/***************************************************************/
/***************************************************************/
/***************************************************************/
#if 0
void DoTheCalculation(MyData *Data, cdouble Omega)
{
   RWGGeometry *G  = Data->G;
   HMatrix *M      = Data->M;
   HMatrix *Mp     = Data->Mp;
   HMatrix *Mm     = Data->Mm;
   HVector *Lambda = Data->Lambda;
   char *FileName  = Data->FileName;
   double Delta    = Data->Delta;
   double Phi      = Data->Phi;

   Log("Assembling M at center frequency...");
   G->AssembleBEMMatrix(Omega, M);

   /*--------------------------------------------------------------*/
   /*- characterization 1: smallest-magnitude eigenvalue          -*/
   /*--------------------------------------------------------------*/
   Log("Computing eigenvalues...");
   Mp->Copy(M);
   Mp->NSEig(Lambda);
   double LambdaMin = abs(Lambda->GetEntry(0));
   double LambdaMax = LambdaMin;
   for(int n=1; n<Lambda->N; n++)
    { LambdaMin = fmin(LambdaMin, abs(Lambda->GetEntry(n)));
      LambdaMax = fmax(LambdaMax, abs(Lambda->GetEntry(n)));
    };

   /*--------------------------------------------------------------*/
   /*--------------------------------------------------------------*/
   /*--------------------------------------------------------------*/
   cdouble DeltaOmega = Delta*abs(Omega)*exp(II*Phi*M_PI/180.0);
   double MNorm = M->GetNorm();
   Log(" Norm of M is %e.",MNorm);
   Log("Assembling at tweaked frequencies pm %s...",z2s(DeltaOmega));
   G->AssembleBEMMatrix(Omega + DeltaOmega, Mp);
   G->AssembleBEMMatrix(Omega - DeltaOmega, Mm);

   /*---------------------------------------------------------------------*/
   /* set Mp = 1st-order (forward) finite-difference approx to dM/dOmega  */
   /* set Mn = 2nd-order (centered) finite-difference approx to dM/dOmega */
   /*---------------------------------------------------------------------*/
   Log("Finite-differencing...");
   for(int nr=0; nr<M->NR; nr++)
    for(int nc=0; nc<M->NC; nc++)
     {  
       cdouble  MEntry =  M->GetEntry(nr,nc);
       cdouble MpEntry = Mp->GetEntry(nr,nc);
       cdouble MmEntry = Mm->GetEntry(nr,nc);

       Mp->SetEntry(nr,nc, (MpEntry - MEntry) / DeltaOmega);
       Mm->SetEntry(nr,nc, (MpEntry - MmEntry) / (2.0*DeltaOmega));
     };
   HMatrix *dMdOmega1 = Mp;
   HMatrix *dMdOmega2 = Mm;

   /*--------------------------------------------------------------*/
   /*--------------------------------------------------------------*/
   /*--------------------------------------------------------------*/
   Log("LU-factorizing...");
   M->LUFactorize();

   double uMin=abs(M->GetEntry(0,0)), uMax=uMin;
   for(int n=1; n<M->NR; n++)
    { uMin = fmin(uMin, abs(M->GetEntry(n,n)) );
      uMax = fmax(uMax, abs(M->GetEntry(n,n)) );
    };

   Log("LU-solving and tracing...");
   M->LUSolve(dMdOmega1);
   M->LUSolve(dMdOmega2);
   cdouble NewtonStep1 = -dMdOmega1->GetTrace();
   cdouble NewtonStep2 = -dMdOmega2->GetTrace();
 
   Log("Fetching RCond(%e)......",MNorm);
   double MRCond = M->GetRCond(MNorm);
   Log("Fetching RCond=(%e).....",MRCond);

   /*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
   /*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
   /*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
   static bool WroteHeader=false;
   if (!WroteHeader)
    { WroteHeader = true;
      FILE *f=fopen(Data->FileName,"a");
      fprintf(f,"# FindBEMPoles running on %s at %s\n",GetHostName(),GetTimeString());
      fprintf(f,"# 1 2 re,imag Omega\n");
      fprintf(f,"# 3 4 Delta Phi\n");
      fprintf(f,"# 5 6 7  LambdaMin LambdaMax LMin/LMax\n");
      fprintf(f,"# 8 9 10 uMin uMax uMin/uMax\n");
      fprintf(f,"# 11     RCond\n");
      fprintf(f,"# 12 13  NewtonStep (1st order)\n");
      fprintf(f,"# 14 15  NewtonStep (2nd order)\n");
      fclose(f);
    };
  
   /*--------------------------------------------------------------*/
   /*--------------------------------------------------------------*/
   /*--------------------------------------------------------------*/
   FILE *f=fopen(FileName,"a");
   fprintf(f,"%e %e ",real(Omega),imag(Omega));
   fprintf(f,"%e %e ",Delta,Phi);
   fprintf(f,"%e %e %e ",LambdaMin,LambdaMax,LambdaMin/LambdaMax);
   fprintf(f,"%e %e %e ",uMin,uMax,uMin/uMax);
   fprintf(f,"%e ",MRCond);
   fprintf(f,"%e %e ",real(NewtonStep1),imag(NewtonStep1));
   fprintf(f,"%e %e ",real(NewtonStep2),imag(NewtonStep2));
   fprintf(f,"\n");
   fclose(f);
}
#endif

/***************************************************************/
/***************************************************************/
/***************************************************************/
void DoTheCalculation(MyData *Data, cdouble Omega)
{
   RWGGeometry *G  = Data->G;
   HMatrix *M      = Data->M;
   HMatrix *Mp     = Data->Mp;
   HMatrix *Mm     = Data->Mm;
   HVector *Lambda = Data->Lambda;
   double *kBloch  = Data->kBloch;
   char *FileName  = Data->FileName;
   double Delta    = Data->Delta;
   double Phi      = Data->Phi;
   bool EstimateDerivative = Data->EstimateDerivative;

   Log("Assembling M at center frequency...");
   G->AssembleBEMMatrix(Omega, kBloch, M);

   double MNorm = M->GetNorm();
   Log("LU-factorizing...");
   M->LUFactorize();

   cdouble LogDet=0.0;
   for(int n=0; n<M->NR; n++)
    LogDet += log(M->GetEntry(n,n));

   double MRCond = M->GetRCond(MNorm);

   cdouble dLogDetEstimated, dLogDetActual;
   if (EstimateDerivative)
    { 
      Log("Assembling M at tweaked frequencies...");
      cdouble DeltaOmega = Delta*abs(Omega)*exp(II*Phi*M_PI/180.0);
      G->AssembleBEMMatrix(Omega+DeltaOmega, kBloch, Mp);
      G->AssembleBEMMatrix(Omega-DeltaOmega, kBloch, Mm);
      for(int nr=0; nr<M->NR; nr++)
       for(int nc=0; nc<M->NC; nc++)
        Mm->SetEntry(nr, nc, (Mp->GetEntry(nr,nc) - Mm->GetEntry(nr,nc)) / (2.0*DeltaOmega));

      M->LUSolve(Mm);
      dLogDetEstimated = Mm->GetTrace();

      Mp->LUFactorize();
      cdouble LogDetp=0.0;
      for(int n=0; n<Mp->NR; n++)
       LogDetp += log(Mp->GetEntry(n,n));

      dLogDetActual=(LogDetp - LogDet) / DeltaOmega;
    };

   /*--------------------------------------------------------------*/
   /*--------------------------------------------------------------*/
   /*--------------------------------------------------------------*/
   static bool WroteHeader=false;
   if (!WroteHeader)
    { WroteHeader = true;
      FILE *f=fopen(Data->FileName,"a");
      fprintf(f,"# FindBEMPoles running on %s at %s\n",GetHostName(),GetTimeString());
      fprintf(f,"# 1 2 re,imag Omega\n");
      fprintf(f,"# 3,4 re, imag LogDet \n");
      fprintf(f,"# 5   RCond\n");
      if (EstimateDerivative)
       { fprintf(f,"# 6 7  8  real(dLogDetActual) imag(dLogDetActual) abs(dLogDetActual)\n");
         fprintf(f,"# 9 10 11 real(dLogDetEstimated) imag(dLogDetEstimated) abs(dLogDetEstimated)\n");
       };
      fclose(f);
    };

   /*--------------------------------------------------------------*/
   /*--------------------------------------------------------------*/
   /*--------------------------------------------------------------*/
   FILE *f=fopen(FileName,"a");
   fprintf(f,"%e %e ",real(Omega),imag(Omega));
   fprintf(f,"%e %e %e ",real(LogDet), imag(LogDet), MRCond);
   if (EstimateDerivative)
    { fprintf(f,"%.2e %.2e (%.2e) ",real(dLogDetActual), imag(dLogDetActual), abs(dLogDetActual));
      fprintf(f,"%.2e %.2e (%.2e) ",real(dLogDetEstimated), imag(dLogDetEstimated), abs(dLogDetEstimated));
    }
   fprintf(f,"\n");
   fclose(f);

   /*--------------------------------------------------------------*/
   /*--------------------------------------------------------------*/
   /*--------------------------------------------------------------*/
   if (Data->Console)
    { printf("At Omega=%s\n",CD2S(Omega));
      printf("LogDet = {%e,%e} \n",real(LogDet),imag(LogDet));
      printf("MRCond = {%e}    \n",MRCond);
      if (EstimateDerivative)
       { printf("dLogDet (actual): {%+.2e,%+.2e} (%.2e)\n",real(dLogDetActual), imag(dLogDetActual), abs(dLogDetActual));
         printf("dLogDet (est)   : {%+.2e,%+.2e} (%.2e)\n",real(dLogDetEstimated), imag(dLogDetEstimated), abs(dLogDetEstimated));
       };
    }
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[]) 
{ 
  InstallHRSignalHandler();
  SetLogFileName("FindBEMPoles.log");
  Log("FindBEMPoles running on %s",GetHostName());

  /***************************************************************/
  /* process command-line arguments ******************************/
  /***************************************************************/
  char *GeoFile=0;
  char *OmegaFile=0;
  char *FileBase=0;
  double Delta=1.0e-4;
  double Phi=45.0;
  double kBloch[2]={0.0,0.0};
  bool EstimateDerivative=false;
  bool Console=false;
  /* name        type    #args  max_instances  storage    count  description*/
  OptStruct OSArray[]=
   { {"geometry",   PA_STRING,  1, 1, (void *)&GeoFile,   0,  ".scuffgeo file"},
//
     {"OmegaFile",  PA_STRING,  1, 1, (void *)&OmegaFile, 0,  "Omega file"},
//
     {"FileBase",   PA_STRING,  1, 1, (void *)&FileBase,  0,  "base name of output file"},
//
     {"Delta",      PA_DOUBLE,  1, 1, (void *)&Delta,     0,  "delta parameter"},
//
     {"EstimateDerivative", PA_BOOL, 0, 1, (void *)&EstimateDerivative,     0,  ""},
     {"Console",    PA_BOOL,    0, 1, (void *)&Console,   0,  ""},
//
     {"Phi",        PA_DOUBLE,  1, 1, (void *)&Phi,       0,  "phi parameter in degrees"},
     {"kBloch",     PA_DOUBLE,  2, 1, (void *)kBloch,     0,  "kBloch"},
//
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);
  if (GeoFile==0)
   OSUsage(argv[0],OSArray,"--geometry option is mandatory");
  if (OmegaFile==0)
   OSUsage(argv[0],OSArray,"--OmegaFile option is mandatory");
  HVector *OmegaVector=new HVector(OmegaFile);
  if (OmegaVector->ErrMsg)
   ErrExit(OmegaVector->ErrMsg);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  SetLogFileName("FindBEMPoles.log","w");
  RWGGeometry *G=new RWGGeometry(GeoFile);
  HMatrix *M=G->AllocateBEMMatrix();
  HMatrix *dM=G->AllocateBEMMatrix();

  struct MyData MyMyData, *Data=&MyMyData;
  Data->G           = G;
  Data->M           = G->AllocateBEMMatrix();
  Data->Mp          = G->AllocateBEMMatrix();
  Data->Mm          = G->AllocateBEMMatrix();
  Data->Lambda      = new HVector(G->TotalBFs, LHM_COMPLEX);
  Data->Delta       = Delta;
  Data->kBloch      = (G->LDim>0) ? kBloch : 0;
  Data->Phi         = Phi;
  Data->EstimateDerivative = EstimateDerivative;
  Data->Console     = Console;
  snprintf(Data->FileName,100,"%s.FBP",FileBase ? FileBase : GetFileBase(GeoFile));

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  for (int no=0; no<OmegaVector->N; no++)
   { cdouble Omega = OmegaVector->GetEntry(no);
     DoTheCalculation(Data, Omega);
   };
  printf("Thank you for your support.\n");

}
