/****************************************************************************************
 Purpose: Handle the prior distribution.

 Author:
   Jinsong Chen and Susan S. Hubbard
   Lawrence Berkeley National Lab
   Earth Sciences Division
   Berkeley, CA 94720

 Written: January of 2003

------------------------------------------------------------------------------------------
 Copyright Notices:

 Bayesian Model (2013-100) (c) 2013, The Regents of the University of California, through
 Lawrence Berkeley National Laboratory (subject to receipt of any required approvals from
 the U.S. Dept. of Energy). All rights reserved.

 Redistribution and use in source and binary forms, with or without modification, are
 permitted provided that the following conditions are met:
 (1) Redistributions of source code must retain the above copyright notice, this list of
     conditions and the following disclaimer.
 (2) Redistributions in binary form must reproduce the above copyright notice, this list
     of conditions and the following disclaimer in the documentation and/or other
     materials provided with the distribution.
 (3) Neither the name of the University of California, Lawrence Berkeley National
     Laboratory, U.S. Dept. of Energy nor the names of its contributors may be used to
     endorse or promote products derived from this software without specific prior written
     permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
 EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
 SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
 OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
 TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 You are under no obligation whatsoever to provide any bug fixes, patches, or upgrades to
 the features, functionality or performance of the source code ("Enhancements") to anyone;
 however, if you choose to make your Enhancements available either publicly, or directly
 to Lawrence Berkeley National Laboratory, without imposing a separate written license
 agreement for such Enhancements, then you hereby grant the following license: a non-
 exclusive, royalty-free perpetual license to install, use, modify, prepare derivative
 works, incorporate into other computer software, distribute, and sublicense such 
 enhancements or derivative works thereof, in binary and source code form.
*****************************************************************************************/
#include "bayes.h"

PRIOR::PRIOR() {}
PRIOR::~PRIOR()
{
  free(x);
  free(z);
  free_dmalloc2(primean);
  free_dmalloc2(privar);
}

int PRIOR::ReadData(FILE *fp)
{
  int i,j,k,errID;
  double tmp1,tmp2,tmp3,tmp4;

  errID=fscanf(fp,"%d %d",&nx,&nz);
  if(errID==0) {
    printf("Prior data file does not include dimensions.\n");
    exit(1);
  };
  x=dmalloc(nx+1);
  z=dmalloc(nz+1);
  primean=dmalloc2(nz+1,nx+1);
  privar=dmalloc2(nz+1,nx+1);

  for(j=1;j<=nz;j++) 
  for(i=1;i<=nx;i++) {
     errID=fscanf(fp,"%lf %lf %lf %lf",&tmp1,&tmp2,&tmp3,&tmp4);
     if(errID==0) return(0);
     if(j==1) x[i]=tmp1;
     if(i==1) z[j]=tmp2;
     primean[i][j]=tmp3;
     privar[i][j]=tmp4;
  };
  return(1);
}

double PRIOR::NormalDist(double tmpx,double mu,double sigma)
{
  double const pi=3.1415926;
  double coef1,coef2,fx;

  coef1=1.0/(sqrt(2*pi)*sigma);
  coef2=(tmpx-mu)*(tmpx-mu)/(2.0*sigma*sigma);
  fx=coef1*exp(-coef2);
  return(fx);
}
