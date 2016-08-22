
/****************************************************************
   Purpose: Header file of Bayesian integration      
   Written by: Jinsong Chen
   Date: JAN-04-2003 
******************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <string.h>
#include <errno.h>

#define ABS(x) ((x)>=0 ? (x) : -(x))
#define min(a,b) ((a)<=(b) ? (a):(b))
#define max(a,b) ((a)>=(b) ? (a):(b))

#define EPSILON 1e-8

FILE *myfopen(const char*, const char*);

double *dmalloc(long);
double **dmalloc2(long,long);
void free_dmalloc2(double**);

class PRIOR
{
  private:

  public:
    PRIOR();
    ~PRIOR();
    int nx,nz;
    double *x,*z;
    double **primean,**privar;
    int ReadData(FILE*);
    double NormalDist(double,double,double);
};

class LKHD
{
  private:
    int nx,nz,nw;
    double *xx,*zz,**gph,*mu_std,*fitdat;
    double lowx,upx,lowz,upz;

  public:
    LKHD();
    ~LKHD();
    int ReadData(FILE*);
    int ReadPetr(FILE*);
    double Likelihood(double,double);
    int Interp2(double,double,double*);
};

class POST
{
  private:
    LKHD lkhd;
    PRIOR prior; 

  public:  
    POST();
    ~POST();
    int InputData(FILE*,FILE*,FILE*);
    int EstimatePost(double,double,double,double*,double*);
    int BayesIntergration(FILE*);
};
