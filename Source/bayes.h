/*****************************************************************************************
 Purpose: Header file of Bayesian integration      

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
