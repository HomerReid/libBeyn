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

#include <nlopt.h>

#define II cdouble(0.0,1.0)

using namespace scuff;

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct MyData
 { 
   RWGGeometry *G;
   HMatrix *M, *Mp, *Mm;
   double *kBloch;
   char FileName[100];
   double Delta;
   double Phi;
   double RCond;

 } MyData;

/***************************************************************/
/***************************************************************/
/***************************************************************/
void AssembleM(MyData *Data, cdouble Omega, double *kBloch, HMatrix *M)
{
  RWGGeometry *G = Data->G;
  int NS         = G->NumSurfaces;

  static void **Accelerators=0;
  if (!Accelerators)
   { Accelerators=(void **)mallocEC(NS*NS*sizeof(void *));
     for(int nsa=0; nsa<NS; nsa++)
      for(int nsb=0; nsb<NS; nsb++)
       Accelerators[nsa*NS + nsb]=CreateABMBAccelerator(nsa, nsb);
   };
 
  for(int nsa=0; nsa<NS; nsa++)
   for(int nsb=0; nsb<NS; nsb++)
    { int RowOffset=G->BFIndexOffset[nsa];
      int ColOffset=G->BFIndexOffset[nsb];
      G=AssembleBEMMatrixBlock(nsa, nsb, Omega, kBloch,
                               M, 0, RowOffset, ColOffset,
                               Accelerator[nsa*NS + nsb]);
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
double Objective(unsigned n, const double* x, double* grad, void* UserData)
{
   MyData *Data = (MyData* )UserData;

   RWGGeometry *G  = Data->G;
   HMatrix *M      = Data->M;
   HMatrix *Mp     = Data->Mp;
   HMatrix *Mm     = Data->Mm;
   double *kBloch  = Data->kBloch;
   char *FileName  = Data->FileName;
   double Delta    = Data->Delta;
   double Phi      = Data->Phi;

   Log("Assembling M at {w,k}={%e,%e} {%e,%e}...",x[0],x[1],kBloch[0],kBloch[1]);
   cdouble Omega = cdouble(x[0],x[1]);
   AssembleM(Data, Omega, kBloch, M);
   M->LUFactorize();
   cdouble LogDetM;
   double MRCond = Data->RCond = GetRCond(M, &LogDetM);
   
   double gradEst[2];
   if (grad)
    {  
      cdouble LogDetP, LogDetM;
      double DeltaOmega = Delta*fabs(x[0]);
      if (DeltaOmega==0.0) DeltaOmega=Delta;
      AssembleM(Data, Omega + DeltaOmega, kBloch, Mp);
      AssembleM(Data, Omega - DeltaOmega, kBloch, Mm);
      
      double MRCondP=GetRCond(Mp,&LogDetP);
      double MRCondM=GetRCond(Mm,&LogDetM);
      grad[0] =    ( MRCondP-MRCondM ) / (2.0*fabs(DeltaOmega));
      gradEst[0] = ( MRCondP-MRCond ) / (fabs(DeltaOmega));

      DeltaOmega = Delta*fabs(x[1]);
      if (DeltaOmega==0.0) DeltaOmega=Delta;
      AssembleM(Data, Omega + II*DeltaOmega, kBloch, Mp);
      AssembleM(Data, Omega - II*DeltaOmega, kBloch, Mm);
      MRCondP=GetRCond(Mp,&LogDetP);
      MRCondM=GetRCond(Mm,&LogDetM);
      grad[1] = ( MRCondP-MRCondM ) / (2.0*fabs(DeltaOmega));
      gradEst[1] = ( MRCondP-MRCond ) / (fabs(DeltaOmega));

    };

   /*--------------------------------------------------------------*/
   /*--------------------------------------------------------------*/
   /*--------------------------------------------------------------*/
   static bool WroteHeader=false;
   if (!WroteHeader)
    { WroteHeader=true;
      FILE *ff=fopen(FileName,"a");
      fprintf(ff,"# Optimize running on %s at %s\n",GetHostName(),GetTimeString());
      fprintf(ff,"# 1 2 re, im Omega\n");
      fprintf(ff,"# 3 4 re, im log det M\n");
      fprintf(ff,"# 5   MRCond\n");
      fprintf(ff,"# 6 7 grad f (coarse)\n");
      fprintf(ff,"# 8 9 grad f (accurate)\n");
      fclose(ff);
    };

   /*--------------------------------------------------------------*/
   /*--------------------------------------------------------------*/
   /*--------------------------------------------------------------*/
   FILE *ff=fopen(FileName,"a");
   fprintf(ff,"%e %e ",x[0],x[1]);
   fprintf(ff,"%e %e ",real(LogDetM),imag(LogDetM));
   fprintf(ff,"%e    ",MRCond);
   fprintf(ff,"%.2e %.2e ",gradEst[0],gradEst[1]);
   if (grad)
    fprintf(ff,"%.2e %.2e ",grad[0],grad[1]);
   else
    fprintf(ff,"0.0 0.0 ");
   fprintf(ff,"\n");
   fclose(ff);

  return MRCond;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[]) 
{ 
  InstallHRSignalHandler();
  SetLogFileName("Optimize.log");
  Log("Optimize running on %s",GetHostName());

  /***************************************************************/
  /* process command-line arguments ******************************/
  /***************************************************************/
  char *GeoFile=0;
  char *OmegaFile=0;
  char *FileBase=0;
  double Delta=1.0e-4;
  double Phi=45.0;
  int MaxEvals=100;
  double fAbsTol=1.0-10;
  double fRelTol=1.0e-2;
  double xRelTol=1.0-2;
  double kBloch[2]={0.0, 0.0};
  bool DerivativeBased=false;
  /* name        type    #args  max_instances  storage    count  description*/
  OptStruct OSArray[]=
   { {"geometry",   PA_STRING,  1, 1, (void *)&GeoFile,   0,  ".scuffgeo file"},
//
     {"OmegaFile",  PA_STRING,  1, 1, (void *)&OmegaFile, 0,  "Omega file"},
//
     {"kBloch",     PA_DOUBLE,  2, 1, (void *)kBloch,    0,   "kBloch"},
     {"FileBase",   PA_STRING,  1, 1, (void *)&FileBase,  0,  "base name of output file"},
//
     {"Delta",      PA_DOUBLE,  1, 1, (void *)&Delta,     0,  "delta parameter"},
     {"Phi",        PA_DOUBLE,  1, 1, (void *)&Phi,       0,  "phi parameter"},
//
     {"MaxEvals",   PA_INT,     1, 1, (void *)&MaxEvals,  0,  "max number of evals"},
//
     {"fAbsTol",    PA_DOUBLE,  1, 1, (void *)&fAbsTol,   0,  "fAbsTol"},
     {"fRelTol",    PA_DOUBLE,  1, 1, (void *)&fRelTol,   0,  "fRelTol"},
     {"xRelTol",    PA_DOUBLE,  1, 1, (void *)&xRelTol,   0,  "xRelTol"},
     {"DerivativeBased", PA_BOOL, 1,1, (void *)&DerivativeBased, 0, ""},
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
  SetLogFileName("Optimize.log","w");
  RWGGeometry *G=new RWGGeometry(GeoFile);
  HMatrix *M=G->AllocateBEMMatrix();
  HMatrix *dM=G->AllocateBEMMatrix();

  struct MyData MyMyData, *Data=&MyMyData;
  Data->G           = G;
  Data->M           = G->AllocateBEMMatrix();
  Data->Mp          = G->AllocateBEMMatrix();
  Data->Mm          = G->AllocateBEMMatrix();
  Data->Delta       = Delta;
  Data->Phi         = Phi;  
  Data->kBloch      = G->LDim > 0 ? kBloch : 0;
  snprintf(Data->FileName,100,"%s.FBP",FileBase ? FileBase : GetFileBase(GeoFile));
  
  nlopt_opt MyOpt;
  if (DerivativeBased)
   MyOpt = nlopt_create(NLOPT_LD_TNEWTON_PRECOND, 2);
  else
   MyOpt = nlopt_create(NLOPT_LN_SBPLX, 2);
  nlopt_set_min_objective(MyOpt, Objective, (void *)Data);
  nlopt_set_ftol_rel(MyOpt, fRelTol);
  nlopt_set_xtol_rel(MyOpt, xRelTol);

  double dx[2]={0.1,0.01};
  nlopt_set_initial_step(MyOpt, dx);

  double Lower[2]={0.0,  -10.0};
  double Upper[2]={10.0, +10.0};
  nlopt_set_lower_bounds(MyOpt, Lower);
  nlopt_set_upper_bounds(MyOpt, Upper);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  for (int no=0; no<OmegaVector->N; no++)
   { 
     cdouble Omega=OmegaVector->GetEntry(no);
     double x[2];
     x[0] = real(Omega);
     x[1] = imag(Omega);

     double fInitial=Objective(2, x, 0, Data);
     double RCInitial=Data->RCond;

     double fBest;
     nlopt_optimize(MyOpt, x, &fBest);
     printf("(%.2e,%.2e; %.2e) --> (%.6e, %.6e; %.2e)\n",
             real(Omega), imag(Omega), fInitial,
             x[0], x[1], fBest);

     double fFinal=Objective(2, x, 0, Data);
     double RCFinal=Data->RCond;
     FILE *f=vfopen("%s.Optimize","a",GetFileBase(GeoFile));
     fprintf(f,"%e %e ",x[0],x[1]);
     if (G->LDim==1)
      fprintf(f,"%e  ",kBloch[0]);
     if (G->LDim==2)
      fprintf(f,"%e %e ",kBloch[0],kBloch[1]);
     fprintf(f,"%e %e ",fInitial,RCInitial);
     fprintf(f,"%e %e ",fFinal, RCFinal);
     fprintf(f,"\n");
     fclose(f);

   };
  printf("Thank you for your support.\n");

}
