/*****************************************************************************************
 Purpose:
   Provide various sub-functions for the program.

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

FILE *myfopen(const char *filename, const char *mode)
{
  FILE *tmp_fp;

  tmp_fp=fopen(filename, mode);
  if(tmp_fp==NULL) {
    fprintf(stderr, "The file %s cannot be opened.\n",filename);
    exit(1);
    };
  return(tmp_fp);
}

double *dmalloc(long num)
{
  double *dptr;

  if((dptr=(double*)malloc(num*sizeof(double)))==NULL) {
    fprintf(stderr, "No enough memories for double pointers.\n");
    exit(1);
    };
  return(dptr);
}

double **dmalloc2(long nrow, long ncol)
{
  long i;
  double **x, *xtemp;
  double *dmalloc(long);

  //Allocate (nrow) double pointers for recording the address of first element at each row
  if((x=(double**)malloc(nrow*sizeof(double*)))==NULL) {
    fprintf(stderr, "No enough memories for double matrix.\n");
    exit(1);
    };

  // Allocate (nrow*ncol) doubles for storing entire element 
  xtemp=dmalloc(nrow*ncol);

  // Connect the x-pointers to the allocated memeories
  for(i=0;i<nrow;i++,xtemp+=ncol) x[i]=xtemp;
  return(x);
}

void free_dmalloc2(double **a)
{
  free(a[0]);
  free(a);
}
