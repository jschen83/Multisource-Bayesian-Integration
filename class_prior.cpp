#include "bayes.h"

PRIOR::PRIOR() {}
PRIOR::~PRIOR()
{
  free(x);
  free(z);
  free_dmalloc2(primean);
  free_dmalloc2(privar);
}

/*
int PRIOR::ReadPara(int tmpnx,int tmpnz)
{
  nx=tmpnx;
  nz=tmpnz;
  x=dmalloc(nx+1);
  z=dmalloc(nz+1);
  primean=dmalloc2(nz+1,nx+1);
  privar=dmalloc2(nz+1,nx+1);
  return(1);
}
*/

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
