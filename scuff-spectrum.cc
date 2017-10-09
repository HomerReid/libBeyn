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
#include "libBeyn.h"

#define II cdouble(0.0,1.0)

using namespace scuff;


/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct BFData
 { 
   RWGGeometry *G;
   HMatrix *M;
   double kBloch[2];
   FILE *LogFile;
 } BFData;

void BeynFunc(cdouble Omega, void *UserData, HMatrix *VHat)
{
  BFData *Data   = (BFData *)UserData;

  RWGGeometry *G = Data->G;
  HMatrix *M     = Data->M;
  double *kBloch = Data->kBloch;
  FILE *LogFile  = Data->LogFile;

  if (G->LDim==0)
   Log(" assembling BEM matrix at Omega=%s",CD2S(Omega));
  if (G->LDim==1)
   Log(" assembling BEM matrix at k={%e},Omega=%s", kBloch[0],CD2S(Omega));
  if (G->LDim==2)
   Log(" assembling BEM matrix at k={%e,%e},Omega=%s", kBloch[0],kBloch[1],CD2S(Omega));

  if (G->LDim==0)
   G->AssembleBEMMatrix(Omega, M);
  else
   G->AssembleBEMMatrix(Omega, kBloch, M);

  if (LogFile)
   fprintf(LogFile,"%e %e\n",real(Omega),imag(Omega));

  Log(" LUFactorizing...");
  M->LUFactorize();
  Log(" LUSolving...");
  M->LUSolve(VHat);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[]) 
{ 
  InstallHRSignalHandler();
  InitializeLog(argv[0]);

  /***************************************************************/
  /* process command-line arguments ******************************/
  /***************************************************************/
  char *GeoFile=0;
  int L = 10;
  char *ContourFile=0;
  char *FileBase=0;
  bool PlotContours=false;
  /* name        type    #args  max_instances  storage    count  description*/
  OptStruct OSArray[]=
   { {"geometry",           PA_STRING,  1, 1, (void *)&GeoFile,            0,   ".scuffgeo file"},
//
     {"L",                  PA_INT,     1, 1, (void *)&L,                  0,   "number of eigenvalues expected inside contour"},
//
     {"ContourFile",        PA_STRING,  1, 1, (void *)&ContourFile,        0,   "list of contours"},
//
     {"FileBase",           PA_STRING,  1, 1, (void *)&FileBase,           0,   "base name for output files"},
//
     {"PlotContours",       PA_BOOL,    0, 1, (void *)&PlotContours,       0,   "plot contours for visualization"},
//
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);
  if (GeoFile==0)
   OSUsage(argv[0],OSArray,"--geometry option is mandatory");
  if (ContourFile==0)
   OSUsage(argv[0],OSArray,"--ContourFileoption is mandatory");
  if (FileBase==0)
   FileBase=strdup(GetFileBase(GeoFile));

  RWGGeometry *G         = new RWGGeometry(GeoFile);
  HMatrix *M             = G->AllocateBEMMatrix();
  HMatrix *ContourMatrix = new HMatrix(ContourFile);

  int N = G->TotalBFs;
  BeynSolver *Solver = CreateBeynSolver(N, L);

  // open output file and write file header 
  FILE *f=vfopen("%s.Frequencies","a",FileBase);
  fprintf(f,"#%s running on %s (%s)\n",argv[0], GetHostName(), GetTimeString());
  int nc=0;
  for(int nc=0; nc<G->LDim; nc++)
   fprintf(f,"# %i k%c\n",nc+1,'x'+nc);
  fprintf(f,"#  %i,%i,%i,%i,%i,%i re,im omega0, Rx, Ry, N, L\n",nc+1,nc+2,nc+3,nc+4,nc+5,nc+6); nc+=6;
  fprintf(f,"#  %i,%i re,im omega1 \n", nc+1, nc+2); nc+=2;
  fprintf(f,"#  %i,%i re,im omega2 \n", nc+1, nc+2); nc+=2;
  fprintf(f,"# ... \n");

  SetDefaultCD2SFormat("%.8e %.8e");
  for(int nr=0; nr<ContourMatrix->NR; nr++)
   {
     struct BFData MyBFData = {G, M};
     double *kBloch = MyBFData.kBloch;

     for(int d=0; d<G->LDim; d++) 
      MyBFData.kBloch[d] = ContourMatrix->GetEntryD(nr, d);
     cdouble Omega0      = ContourMatrix->GetEntry (nr, G->LDim+0);
     double Rx           = ContourMatrix->GetEntryD(nr, G->LDim+1);
     double Ry           = ContourMatrix->GetEntryD(nr, G->LDim+2);
     int N               = (int)(ContourMatrix->GetEntryD(nr, G->LDim+3));
     if (Rx==Ry)
      Log("Looking for eigenvalues in contour (z0,R)={%s,%e}",CD2S(Omega0),Rx);
     else
      Log("Looking for eigenvalues in contour (z0,Rx,Ry)={%s,%e,%e}",CD2S(Omega0),Rx,Ry);

     if (PlotContours)
      { MyBFData.LogFile = vfopen("%s.contours","a",FileBase);
        fprintf(MyBFData.LogFile,"\n\n# (Omega0,Rx,Ry,N)=%s,%e,%e,%i\n",CD2S(Omega0),Rx,Ry,N);
      };

     int K=BeynSolve(Solver, BeynFunc, (void *)&MyBFData, Omega0, Rx, Ry, N);
     if (G->LDim>0)
      fprintVec(f,kBloch,G->LDim);
     fprintf(f,"%s %+e %+e %3i %3i ",CD2S(Omega0),Rx,Ry,N,L);
     fprintVecCR(f,Solver->Lambda->ZV,K);
     fflush(f);

     if (PlotContours)
      fclose(MyBFData.LogFile);

   }; 
  fclose(f);

}
