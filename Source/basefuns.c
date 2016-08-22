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

  /*allocate (nrow) double pointers for recording the
  address of the first element at each row*/
  if((x=(double**)malloc(nrow*sizeof(double*)))==NULL) {
    fprintf(stderr, "No enough memories for double matrix.\n");
    exit(1);
    };

  /*allocate (nrow*ncol) doubles for storing entire element*/ 
  xtemp=dmalloc(nrow*ncol);

  /*connect the x-pointers to the allocated memeories*/ 
  for(i=0;i<nrow;i++,xtemp+=ncol) x[i]=xtemp;
  return(x);
}

void free_dmalloc2(double **a)
{
  free(a[0]);
  free(a);
}
