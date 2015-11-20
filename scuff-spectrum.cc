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
#include <scuff-spectrum.h>

#define II cdouble(0.0,1.0)

using namespace scuff;

/***************************************************************/
/***************************************************************/
/***************************************************************/
void AdjustkBloch(RWGGeometry *G, cdouble Omega, double *kBloch)
{
  if (imag(Omega)!=0.0) return;

  double Threshold=1.0e-4;
  if (char *s=getenv("SCUFF_KBLOCH_THRESHOLD"))
   sscanf(s,"%le",&Threshold);
  if (Threshold==0.0)
   return;
  double T2=Threshold*Threshold;

  double kB2 = kBloch[0]*kBloch[0];
  if (G->LDim==2) kB2+=kBloch[1]*kBloch[1];
  if (kB2==0.0) return;

  for(int nr=0; nr<G->NumRegions; nr++)
   { 
     cdouble k=Omega*G->RegionMPs[nr]->GetRefractiveIndex(Omega);
     double k2=norm(k);

     if ( fabs( k2 - kB2 ) > T2*kB2 )
      continue;

     Warn("|kBloch| very close to |k| in region %i; tweaking by 1 %c %e",
           nr, Threshold > 0.0 ? '+' : '-' , fabs(Threshold));
     if (G->LDim==1) 
      kBloch[0]*=(1.0+Threshold);
     else
      { kBloch[0]*=(1.0+Threshold) / M_SQRT2;
        kBloch[1]*=(1.0+Threshold) / M_SQRT2;
      };
     return;
   };

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void Optimize(SSData *Data, cdouble Omega0)
{
  /***************************************************************/ 
  /***************************************************************/ 
  /***************************************************************/ 
  double fRelTol=1.0e-2;
  double xRelTol=1.0e-2;
  static nlopt_opt MyOpt=0;
  if (MyOpt==0)
   { 
     if (Data->DBOptimization)
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
   };

  /***************************************************************/ 
  /***************************************************************/ 
  /***************************************************************/ 
  double x[2];
  x[0] = real(Omega0);
  x[1] = imag(Omega0);

  double fInitial=Objective(2, x, 0, Data);
  double RCInitial=Data->RCond;

  double fBest;
  nlopt_optimize(MyOpt, x, &fBest);
  printf("(%.2e,%.2e; %.2e) --> (%.6e, %.6e; %.2e)\n",
          real(Omega0), imag(Omega0), fInitial, x[0], x[1], fBest);

  double fFinal=Objective(2, x, 0, Data);
  double RCFinal=Data->RCond;

  FILE *f=vfopen("%s.Optimize","a",Data->FileBase);
  double *kBloch=Data->kBloch;
  fprintf(f,"%e %e ",x[0],x[1]);
  if (Data->G->LDim==1)
   fprintf(f,"%e  ",kBloch[0]);
  if (Data->G->LDim==2)
   fprintf(f,"%e %e ",kBloch[0],kBloch[1]);
  fprintf(f,"%e %e ",fInitial,RCInitial);
  fprintf(f,"%e %e ",fFinal, RCFinal);
  fprintf(f,"\n");
  fclose(f);

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[]) 
{ 
  InstallHRSignalHandler();
  SetLogFileName("scuff-spectrum.log");
  Log("scuff-spectrum running on %s",GetHostName());

  /***************************************************************/
  /* process command-line arguments ******************************/
  /***************************************************************/
  char *GeoFile=0;
  cdouble OmegaValues[10]; int nOmega=0;
  char *OmegaFile=0;
  char *OmegakBlochFile=0;
  char *FileBase=0;
  double Delta=1.0e-4;
  double Phi=45.0;
  double kBloch[2]={0.0,0.0};
  bool EstimateDerivative=false;
  bool Console=false;
  bool SecondOrderFD=false;
  bool DBOptimization=false;
  bool RCond=false;
  bool Eigenvalues=false;
  bool SVD=false;
  /* name        type    #args  max_instances  storage    count  description*/
  OptStruct OSArray[]=
   { {"geometry",           PA_STRING,  1, 1, (void *)&GeoFile,            0,   ".scuffgeo file"},
//
     {"kBloch",             PA_DOUBLE,  2, 1, (void *)kBloch,              0,  "kBloch"},
//
     {"Omega",              PA_CDOUBLE, 1, 9, (void *)OmegaValues,     &nOmega, ""},
     {"OmegaFile",          PA_STRING,  1, 1, (void *)&OmegaFile,          0,   ""},
     {"OmegakBlochFile",    PA_STRING,  1, 1, (void *)&OmegakBlochFile,    0,   ""},
//
     {"FileBase",           PA_STRING,  1, 1, (void *)&FileBase,           0,   ""},
//
     {"Delta",              PA_DOUBLE,  1, 1, (void *)&Delta,              0,   ""},
//
     {"EstimateDerivative", PA_BOOL,    0, 1, (void *)&EstimateDerivative, 0,   ""},
     {"SecondOrderFD",      PA_BOOL,    0, 1, (void *)&SecondOrderFD,      0,   ""},
     {"Console",            PA_BOOL,    0, 1, (void *)&Console,            0,   ""},
     {"DBOptimization",     PA_BOOL,    0, 1, (void *)&DBOptimization,     0,   ""},
     {"RCond",              PA_BOOL,    0, 1, (void *)&RCond,              0,   ""},
     {"Eigenvalues",        PA_BOOL,    0, 1, (void *)&Eigenvalues,        0,   ""},
     {"SVD",                PA_BOOL,    0, 1, (void *)&SVD,                0,   ""},
//
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);
  if (GeoFile==0)
   OSUsage(argv[0],OSArray,"--geometry option is mandatory");

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  HVector *OmegaVector=0;
  if (nOmega>0)
   { OmegaVector=new HVector(nOmega, LHM_COMPLEX);
     for(int no=0; no<nOmega; no++)
      OmegaVector->SetEntry(no,OmegaValues[no]);
   }
  else if (OmegaFile)
   { OmegaVector=new HVector(OmegaFile);
     if (OmegaVector->ErrMsg)
      ErrExit(OmegaVector->ErrMsg);
   };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  RWGGeometry *G=new RWGGeometry(GeoFile);

  struct SSData MyData, *Data=&MyData;
  Data->G           = G;
  Data->M           = G->AllocateBEMMatrix();
  Data->Mp          = G->AllocateBEMMatrix();
  Data->Mm          = G->AllocateBEMMatrix();
  Data->Lambda      = G->AllocateRHSVector();
  Data->Sigma       = G->AllocateRHSVector(true);
  Data->Delta       = Delta;
  Data->kBloch      = (G->LDim>0) ? kBloch : 0;
  Data->EstimateDerivative = EstimateDerivative;
  Data->Console     = Console;
  Data->SecondOrderFD = SecondOrderFD;
  Data->DBOptimization = DBOptimization;
  Data->DoRCond     = RCond;
  Data->Eigenvalues = Eigenvalues;
  Data->SVD         = SVD;
  Data->FileBase    = FileBase ? FileBase : strdup(GetFileBase(GeoFile));

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (OmegakBlochFile)
   { HMatrix *OkBMatrix=new HMatrix(OmegakBlochFile);
     if (OkBMatrix->ErrMsg)
      ErrExit(OkBMatrix->ErrMsg);
     Data->kBloch=kBloch;
     for (int n=0; n<OkBMatrix->NR; n++)
      { Data->Omega = OkBMatrix->GetEntry(n, 0);
        kBloch[0]   = OkBMatrix->GetEntryD(n, 1);
        if (OkBMatrix->NC>=3)
         kBloch[1]   = OkBMatrix->GetEntryD(n, 2);
        else
         kBloch[1]   = 0.0;

        AdjustkBloch(G, Data->Omega, kBloch);

        Objective(0, 0, 0, Data);
      };
   }
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  else if (OmegaVector)
   for (int n=0; n<OmegaVector->N; n++)
    Optimize(Data,OmegaVector->GetEntry(n));

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  printf("Thank you for your support.\n");
}

