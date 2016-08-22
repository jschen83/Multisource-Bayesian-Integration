#include "bayes.h"

POST::POST() {}
POST::~POST() {}

int POST::InputData(FILE *fp1,FILE *fp2,FILE *fp3)
{
  int errID;

  errID=prior.ReadData(fp1);
  if(errID==0) {
    printf("Unsuccessful when reading the PRIOR data file\n");
    return(0);
  };

  errID=lkhd.ReadData(fp2);
  if(errID==0) {
    printf("Unsuccessful when reading the tomographic data file\n");
    return(0);
  };

  errID=lkhd.ReadPetr(fp3);
  if(errID==0) {
    printf("Unsuccessful when reading the petrophysics file\n");
    return(0);
  };

  return(1);
}

int POST::EstimatePost(double v,double prior_mu,
                       double prior_std,double *post_mu,
                       double *post_std)
{
  int i,n;
  double lx,hx,h,tmpstd,priorfx,likelihood,*x,*post;
  double area,moment1,moment2;

  /*Suppose a wide range which is enough to cover
    most interested range*/

  n=1001;
  tmpstd=prior_std;
  lx=prior_mu-3*tmpstd;
  hx=prior_mu+3*tmpstd;
  x=dmalloc(n+1);
  post=dmalloc(n+1);
  h=(hx-lx)/(n-1);
  for(i=1;i<=n;i++) x[i]=lx+(i-1)*h;

  for(i=1;i<=n;i++) {
     priorfx=prior.NormalDist(x[i],prior_mu,prior_std);
     likelihood=lkhd.Likelihood(v,x[i]);
     post[i]=priorfx*likelihood; //Need normalizing
  };

  /*Calculate psuedo poster distribution
     mean and second moments*/
   area=0.5*h*(post[1]+post[n]);
   moment1=0.5*h*(post[1]*x[1]+post[n]*x[n]);
   moment2=0.5*h*(post[1]*x[1]*x[1]+post[n]*x[n]*x[n]);

  for(i=2;i<=n-1;i++) {
     area += h*post[i];
     moment1 += h*post[i]*x[i];
     moment2 += h*post[i]*x[i]*x[i];
  };

  /*Calculate posterior mean and standard deviation*/
   (*post_mu)=moment1/area;
   (*post_std)=sqrt(moment2/area-moment1*moment1/(area*area));

  free(x);
  free(post);
  return(1);
} 

int POST::BayesIntergration(FILE *fp)
{
  int gph_ok,errID,i,j,pnz,pnx;
  double xloc,zloc,prior_mu,prior_std,post_mu,post_std;
  double gph_val,prior_var;

  pnz=prior.nz;
  pnx=prior.nx;

  for(j=1;j<=pnz;j++) 
  for(i=1;i<=pnx;i++) {
     xloc=prior.x[i];
     zloc=prior.z[j];
     gph_ok=lkhd.Interp2(xloc,zloc,&gph_val);
     prior_mu=prior.primean[i][j];
     prior_var=prior.privar[i][j];

     if(prior_var<1e-6) 
       prior_std=0.0;
     else
       prior_std=sqrt(prior_var);

     if(gph_ok==0 || prior_std<1e-4) {
       post_mu=prior_mu;
       post_std=prior_std;
     }
     else {
       errID=EstimatePost(gph_val,prior_mu,
                          prior_std,&post_mu,&post_std);
       if(errID==0) {
         printf("Unsuccessful in calling EstimatePost\n");
         return(0);
       };  
     };
     fprintf(fp,"%8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n",
                xloc,zloc,prior_mu,prior_std,post_mu,post_std);

  };
  return(1);
}
