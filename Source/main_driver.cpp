/****************************************************************************************
 Purpose: Main function for Bayesian integration of geophysical and hydrogeological data.

 Author:
   Jinsong Chen and Susan S. Hubbard
   Lawrence Berkeley National Lab
   Earth Sciences Division
   Berkeley, CA 94720

 Written: January of 2003
 Updated: September of 2005

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

int main(void)
{
  FILE *fp1,*fp2,*fp3;
  char fname[80];
  POST post;

  printf("Please input the file name of prior data:");
  scanf("%s",fname);
  fp1=myfopen(fname,"r");

  printf("Input the file name of tomographic data:");
  scanf("%s",fname);
  fp2=myfopen(fname,"r");

  printf("Input the file name of petrophysical models:");
  scanf("%s",fname);
  fp3=myfopen(fname,"r");

  post.InputData(fp1,fp2,fp3);
  fclose(fp1);
  fclose(fp2);
  fclose(fp3);
  
  printf("Input the file name of output:");
  scanf("%s",fname);
  fp1=myfopen(fname,"w");

  post.BayesIntergration(fp1);
  fclose(fp1);
}
