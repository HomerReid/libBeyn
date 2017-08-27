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
#include "MoreUtils.h"

#define II cdouble(0.0,1.0)

using namespace scuff;
namespace scuff {
cdouble GetG(double R[3], cdouble k, cdouble *dGBar, cdouble *ddGBar);
                }

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[]) 
{ 
  InstallHRSignalHandler();

  /***************************************************************/
  /* process command-line arguments ******************************/
  /***************************************************************/
  srand48(time(0));
  cdouble k=cdouble(drand48(), drand48());
  double R[3];
  double kBloch[2];
  R[0]=drand48();
  R[1]=drand48();
  R[2]=drand48();
  kBloch[0]=drand48();
  kBloch[1]=drand48();
  double Delta=1.0e-4;
  int LDim=2;
  /* name        type    #args  max_instances  storage    count  description*/
  OptStruct OSArray[]=
 {
     {"k",              PA_CDOUBLE, 1, 1, (void *)&k,              0, ""},
//
     {"R",                  PA_DOUBLE,  3, 1, (void *)R,                   0, ""},
//
     {"kBloch",             PA_DOUBLE,  2, 1, (void *)kBloch,              0, ""},
     {"Delta",              PA_DOUBLE,  1, 1, (void *)&Delta,              0, ""},
     {"LDim",               PA_INT,     1, 1, (void *)&LDim,              0, ""},
//
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);

  double LBV1[2]={1.0,0.0};
  double LBV2[2]={0.0,1.0};
  double *LBV[2]={LBV1, LBV2};

  GBarAccelerator *GBA=0;
  if (LDim>0)
   { GBA =new GBarAccelerator;
     GBA->k=k;
     GBA->kBloch=kBloch;
     GBA->ExcludeInnerCells=false;
     GBA->LDim=LDim;
     GBA->LBV1[0]=1.0; GBA->LBV1[1]=0.0; 
     GBA->LBV2[0]=0.0; GBA->LBV2[1]=1.0; 
     GBA->LBV[0]=GBA->LBV1;
     GBA->LBV[1]=GBA->LBV2;
     GBA->ForceFullEwald=true;
     GBA->I2D=0;
     GBA->I3D=0;
   };

  double r2=R[0]*R[0] + R[1]*R[1] + R[2]*R[2];
  double r=sqrt(r2);
  cdouble dG[3];
  cdouble G0 = GBA ? GetGBar(R, GBA, dG, 0, true) : GetG(R, k, dG, 0);

  cdouble dGdk[8];
  dGdk[0] = -II*r2*( (R[0]*dG[0]+R[1]*dG[1]+R[2]*dG[2])/r - II*k*G0);
  dGdk[1] = -k*G0*R[0];
  dGdk[2] = -k*G0*R[1];
  dGdk[3] = -k*G0*R[2];

  cdouble dk=Delta*k;
  cdouble GpdG, dGpdG[3], GmdG, dGmdG[3];
  if (GBA)
   { 
     GBA->k = k+dk;
     GpdG=GetGBar(R, GBA, dGpdG, 0, true);

     GBA->k = k-dk;
     GmdG=GetGBar(R, GBA, dGmdG, 0, true);
   }
  else
   { 
     GpdG=GetG(R, k+dk, dGpdG, 0);
     GmdG=GetG(R, k-dk, dGmdG, 0);
   };
  cdouble dGdkBF[8];
  dGdkBF[0] = (GpdG - GmdG) / (2.0*dk);
  for(int n=0; n<3; n++)
   dGdkBF[n+1] = (dGpdG[n] - dGmdG[n]) / (2.0*dk);

  Compare(dGdk, dGdkBF, 4, "HR", "BF");

}
