/*****************************************************************
  Purpose: Main function for Bayesian integration of geophysical
           and hydrogeological data.
  Written: Sept. 2, 2005                                        
  Author: J. Chen                                             
******************************************************************/ 
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
