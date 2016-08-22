/****************************************************************************************
 Purpose: Likelihood functions

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

LKHD::LKHD(){}
LKHD::~LKHD()
{
  free(zz);
  free(xx);
  free_dmalloc2(gph);
  free(mu_std);
  free(fitdat);
}

int LKHD::ReadData(FILE *fp)
{ 
  int i,j,errID;
  double tmp1,tmp2,tmp3;

  errID=fscanf(fp, "%d %d",&nx,&nz);
  if(errID==0) {
    printf("The tomographic data file does not include dimensions\n");
    exit(1);
  };
  xx=dmalloc(nx+1);
  zz=dmalloc(nz+1);
  gph=dmalloc2(nz+1,nx+1);

  for(j=1;j<=nz;j++) 
  for(i=1;i<=nx;i++) {
     errID=fscanf(fp, "%lf %lf %lf",&tmp1,&tmp2,&tmp3);
     if(errID==0) return(0);
     if(j==1) xx[i]=tmp1;
     if(i==1) zz[j]=tmp2;
     gph[j][i]=tmp3;
  };

  lowx=xx[1];
  upx=xx[1];
  for(i=2;i<=nx;i++) {
     if(xx[i]<lowx) lowx=xx[i];
     if(xx[i]>upx) upx=xx[i];
  };

  lowz=zz[1];
  upz=zz[1];
  for(i=2;i<=nz;i++) {
     if(zz[i]<lowz) lowz=zz[i];
     if(zz[i]>upz) upz=zz[i];
  };
     
  return(1); 
}

int LKHD::ReadPetr(FILE *fp)
{
  int i,j,errID;
  char dummy_str[80];
  double tmp;

  /*
  errID=fscanf(fp,"%d",&nw);
  if(errID==0) {
    printf("Petrophysics model file does not include dimensions\n");
    exit(1);
  };
  */

  nw=1; //for HydroImage at this phase

  mu_std=dmalloc(9);
  fitdat=dmalloc(6);

  /*
  for(i=1;i<=nw;i++) {
     errID=fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf %lf %lf",
        tmp,tmp+1,tmp+2,tmp+3,tmp+4,tmp+5,tmp+6,tmp+7,tmp+8);
     if(errID==0) return(0);
     for(j=1;j<=8;j++) mu_std[i][j]=tmp[j];
  };
  */
  for(i=1;i<=8;i++) {
     fscanf(fp,"%s %lf",dummy_str,&tmp);
     mu_std[i]=tmp;
     //printf("%s %f\n",dummy_str,mu_std[i]);
  };

 fscanf(fp,"%s",dummy_str);//Regression formula
 fscanf(fp,"%s",dummy_str);//Regression formula

  /*Fitted data from regression models*/
  /*
  for(i=1;i<=nw;i++) {
     errID=fscanf(fp,"%lf %lf %lf %lf %lf %lf",
                  tmp,tmp+1,tmp+2,tmp+3,tmp+4,tmp+5);
     if(errID==0) return(0);
     for(j=1;j<=5;j++) fitdat[i][j]=tmp[j];
  };
 */

  for(i=1;i<=5;i++) {
     fscanf(fp,"%s %lf",dummy_str,&tmp);
     fitdat[i]=tmp;
     //printf("%d %s %lf\n",i,dummy_str,fitdat[i]);
  };

  return(1);
}

double LKHD::Likelihood(double gv,double hv)
{
  int j;
  double const pi=3.1415926;
  double hv_mu,hv_std,hv_norm,gv_mu,gv_std,gv_mean;
  double beta[4],S_square,mu,sq_sigma,sigma,coef1,coef2,fvk;
 
  hv_mu=mu_std[1];
  hv_std=mu_std[2];
  gv_mu=mu_std[3];
  gv_std=mu_std[4];

  for(j=0;j<4;j++) beta[j]=fitdat[j+1];
  S_square=fitdat[5];

  hv_norm=(hv-hv_mu)/hv_std;
  gv_mean=beta[0]+beta[1]*hv_norm+beta[2]*hv_norm*hv_norm+
          beta[3]*hv_norm*hv_norm*hv_norm;

  /*Transfer to original scale*/
  mu=gv_mu+gv_mean*gv_std;
  sq_sigma=S_square*gv_std*gv_std;

  sigma=sqrt(sq_sigma);
  coef1=1.0/(sqrt(2*pi)*sigma);
  coef2=(gv-mu)*(gv-mu)/(2.0*sq_sigma);
  fvk=coef1*exp(-coef2);
  return(fvk);
}

int LKHD::Interp2(double xloc,double zloc,double *rtn)
{
  int lftdwn,rghtdwn,lftup,rghtup;
  int i,j,ki,kj;
  double ri,rj,tmpki1,tmpki2,tmpkj1,tmpkj2,tmp;

  /*whether (xloc,zloc) is out of the bounds?*/
  if(xloc>upx || xloc<lowx || zloc>upz || zloc<lowz) return(0);

  /*Here (ki,kj) is the index of the left-down corner of each
    pixel. Be careful: it is not common order but consistent
    with my own datasets.*/
    
  for(ki=1,i=2;i<=nz-1;i++) if(zloc>zz[i]) ki=i;
  for(kj=1,j=2;j<=nx-1;j++) if(xloc>xx[j]) kj=j;

  ri=(zloc-zz[ki])/(zz[ki+1]-zz[ki]);
  rj=(xloc-xx[kj])/(xx[kj+1]-xx[kj]);
  tmpki1=gph[ki][kj]+rj*(gph[ki][kj+1]-gph[ki][kj]);
  tmpki2=gph[ki+1][kj]+rj*(gph[ki+1][kj+1]-gph[ki+1][kj]);
  tmpkj1=gph[ki][kj]+ri*(gph[ki+1][kj]-gph[ki][kj]);
  tmpkj2=gph[ki][kj+1]+ri*(gph[ki+1][kj+1]-gph[ki][kj+1]);
  tmp=tmpki1+ri*(tmpki2-tmpki1);

  /*Check validity of values at each corner*/
  lftdwn=(gph[ki][kj]!=0.0? 1:0);
  lftup=(gph[ki+1][kj]!=0.0? 1:0);
  rghtdwn=(gph[ki][kj+1]!=0.0? 1:0);
  rghtup=(gph[ki+1][kj+1]!=0.0? 1:0);

  /*If values at all four corner are valid, interpolate directly*/
  if(lftdwn==1 && lftup==1 && rghtdwn==1 && rghtup==1) {
    (*rtn)=tmp; //Valid at four corners
    return(1);
  }; 

  /*If values at all four corner are invalid*/
  if(lftdwn==0 && lftup==0 && rghtdwn==0 && rghtup==0) {
    return(0); 
  };

  /*Only one value is valid*/
  if(lftdwn==1 && lftup==0 && rghtdwn==0 && rghtup==0) {
    (*rtn)=gph[ki][kj];
    return(1);
  }
  else if(lftdwn==0 && lftup==1 && rghtdwn==0 && rghtup==0) {
    (*rtn)=gph[ki+1][kj];
    return(1);
  } 
  else if(lftdwn==0 && lftup==0 && rghtdwn==1 && rghtup==0) {
    (*rtn)=gph[ki][kj+1];
    return(1);
  }
  else if(lftdwn==0 && lftup==0 && rghtdwn==0 && rghtup==1) {
    (*rtn)=gph[ki+1][kj+1];
    return(1);
  };

  /*Only two values are valid*/
  if(lftdwn==1 && lftup==1 && rghtdwn==0 && rghtup==0) {
    (*rtn)=tmpkj1;
    return(1);
  }
  else if(lftdwn==1 && lftup==0 && rghtdwn==1 && rghtup==0) {
    (*rtn)=tmpki1;
    return(1);
  }
  else if(lftdwn==1 && lftup==0 && rghtdwn==0 && rghtup==1) {
    (*rtn)=0.5*(gph[ki][kj]+gph[ki+1][kj+1]);
    return(1);
  }
  else if(lftdwn==0 && lftup==1 && rghtdwn==1 && rghtup==0) {
    (*rtn)=0.5*(gph[ki+1][kj]+gph[ki][kj+1]);
    return(1);
  }
  else if(lftdwn==0 && lftup==1 && rghtdwn==0 && rghtup==1) {
    (*rtn)=tmpki2;
    return(1);
  } 
  else if(lftdwn==0 && lftup==0 && rghtdwn==1 && rghtup==1) {
    (*rtn)=tmpkj2;
    return(1);
  };

  /*Three values are valid*/
  if(lftdwn==1 && lftup==1 && rghtdwn==1 && rghtup==0) {
    (*rtn)=(gph[ki][kj]+gph[ki+1][kj]+gph[ki][kj+1])/3;
    return(1);
  } 
  else if(lftdwn==1 && lftup==1 && rghtdwn==0 && rghtup==1) {
    (*rtn)=(gph[ki][kj]+gph[ki+1][kj]+gph[ki+1][kj+1])/3;
    return(1);
  } 
  else if(lftdwn==1 && lftup==0 && rghtdwn==1 && rghtup==1) {
    (*rtn)=(gph[ki][kj]+gph[ki][kj+1]+gph[ki+1][kj+1])/3;
    return(1);
  } 
  else if(lftdwn==0 && lftup==1 && rghtdwn==1 && rghtup==1) {
    (*rtn)=(gph[ki+1][kj]+gph[ki][kj+1]+gph[ki+1][kj+1])/3;
    return(1);
  }; 
}
