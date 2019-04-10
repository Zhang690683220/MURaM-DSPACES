#ifndef __ATOM__   // __ATOM__
#define __ATOM__

#include<string.h>
#define Hperg 4.407e23;
#define k 1.3806503e-16
#define amu 1.66053892173E-24
#define ev 1.60218e-12
#define ev_kj_moli 96.485
#define inz_H2 15.422086
#define inz_H 13.59844
#define inz_Hm1 0.754
#define D0_h2 4.478
#define D0_h2p 2.654
  

using std::string;

class atom{
 public:
  
      int na;
      int static const nH=9;
      char elm[2];
      string *ATM,*ATM2;
      int *Z;
      double *chi1,*chi2;
      double *W,*mass; 
      double *abu, *perg;
      double *EH,*gH;
      double muavg;
      double abutot;

      double *uu1,*uu2,*uu3;
      atom(int,double*);
      ~atom();
     void partf(int,double);
};

class molecb{
 public:
    int nm; 
    double *Y,*dY;
    molecb(int,double);
   ~molecb();
};

int * sort(double * list, int size);

#endif                // __ATOM__

