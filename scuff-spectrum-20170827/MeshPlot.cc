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
#include <libhmat.h>

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[]) 
{ 
  InstallHRSignalHandler();
  SetLogFileName("MeshPlot.log");

  /***************************************************************/
  /* process command-line arguments ******************************/
  /***************************************************************/
  char *DataFile=0;
  int xCol  = 1;
  int yCol  = 2;
  int zCols[10];
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
     {"Phi",        PA_DOUBLE,  1, 1, (void *)&Phi,       0,  "phi parameter in degrees"},
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
  Data->Phi         = Phi;
  snprintf(Data->FileName,100,"%s.FBP",FileBase ? FileBase : GetFileBase(GeoFile));

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
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

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  for (int no=0; no<OmegaVector->N; no++)
   { cdouble Omega = OmegaVector->GetEntry(no);
     DoTheCalculation(Data, Omega);
   };
  printf("Thank you for your support.\n");

}
