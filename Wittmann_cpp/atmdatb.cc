/* ATMDATB is the modified ATMDAT routine so that it outputs not only U1, U2, U3, but also the derivatives of
the logarithms of said energy levels with respect to temperature DU1, DU2, DU3 (except items 39 to 77) plus 56.
the rest are taken as 0 at the moment.
SUPPLY ATOMIC PARAMETERS FOR 92 ELEMENTS (HYDROGEN TO URANIUM)
A.D. WITTMANN, GOETTINGEN (1975) */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <iostream>

#include "atom.h"

int * sort(double * list, int size)
{
  int * ind = new int [size];
  for(int i=0; i<size; i++)
    ind[i] = i;

  for(int i=0; i<size; i++)
  {
    for(int j=size-1; j>i; j--)
    {
      if(list[j]>list[j-1])
      {
        double swap=list[j-1];
        int swapint = ind[j-1];

        list[j-1]=list[j];
        ind[j-1] = ind[j];

        list[j]= swap;
        ind[j] = swapint;
      }
    }
  }
  return ind;
}


using std::string;
  /* Names of the elements in capital letters */
  string ATM0[]={"H","HE","LI","BE","B","C","N","O","F","NE",
    "NA","MG","AL","SI","P","S","CL","AR","K","CA","SC","TI","V","CR",
    "MN","FE","CO","NI","CU","ZN","GA","GE","AS","SE","BR","KR",
    "RB","SR","Y","ZR","NB","MO","TC","RU","RH","PD","AG","CD","IN",
    "SN","SB","TE","I","XE","CS","BA","LA","CE","PR","ND","PM",
    "SM","EU","GD","TB","DY","HO","ER","TM","YB","LU","HF","TA","W",
    "RE","OS","IR","PT","AU","HG","TL","PB","BI","PO","AT","RN",
    "FR","RA","AC","TH","PA","U"};

/*     Names of the elements in small letters  */
      string ATM20[]={"h","he","li","be","b","c","n","o","f","ne",
     "na","mg","al","si","p","s","cl","ar","k","ca","sc","ti","v","cr",
     "mn","fe","co","ni","cu","zn","ga","ge","as","se","br","kr",
     "rb","sr","y","zr","nb","mo","tc","ru","rh","pd","ag","cd","in",
     "sn","sb","te","i","xe","cs","ba","la","ce","pr","nd","pm",
     "sm","eu","gd","tb","dy","ho","er","tm","yb","lu","hf","ta","w",
     "re","os","ir","pt","au","hg","tl","pb","bi","po","at","rn",
     "fr","ra","ac","th","pa","u"};

//     Here E1 and E2 are just dimension statements and W is the atomic mass of 92 elements
      double W0[]={1.008,4.003,6.939,9.012,10.811,12.011,14.007,16.,18.998,
     20.183,22.99,24.312,26.982,28.086,30.974,32.064,35.453,39.948,39.102,
     40.08,44.956,47.90,50.942,51.996,54.938,55.847,58.933,58.71,
     63.54,65.37,69.72,72.59,74.92,78.96,79.91,83.80,85.47,87.62,88.905,
     91.22,92.906,95.94, 99.00,101.07,102.9,106.4,107.87,112.40,114.82,
     118.69,121.75,127.6,126.9,131.3,132.9,137.34,138.91,140.12,140.91,
     144.24,147.00,150.35,151.96,157.25,158.92,162.50,164.93,
     167.26,168.93,173.04,174.97,178.49,180.95,183.85,186.2,190.2,192.2,
     195.09,196.97,200.59,204.37,207.19,208.98,210.,211.,222.,
     223.,226.1,227.1,232.04,231.,238.03};

//     Here E1 and E2 are just dimension statements and W is the atomic mass of 92 elements
      int Z0[]={1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,
        24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,
        47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,
        70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92};
     

//   Ionization energies of 92 elements in eV for outermost electron
      double E10[]={13.595,24.58,5.39,9.32,8.298,11.256,14.529,13.614,17.418,
     21.559,5.138,7.644,5.984,8.149,10.474,10.357,13.012,15.755,4.339,
     6.111,6.538,6.825,6.738,6.763,7.432,7.896,7.863,7.633,7.724,9.391,
     5.997,7.88,9.81,9.75,11.840,13.996,4.176,5.692,6.377,6.838,6.881,
     7.10,7.28,7.36,7.46,8.33,7.574,8.991,5.785,7.34,8.64,9.01,10.454,
     12.127,3.893,5.210,5.577,5.466,5.422,5.489,5.554,5.631,5.666,6.141,
     5.852,5.927,6.018,6.101,6.184,6.254,5.426,6.650,7.879,7.980,7.870,
     8.70,9.10,9.00,9.22,10.43,6.105,7.415,7.287,8.43,9.30,10.745,
     4.,5.276,6.9,6.,6.,6.};

//    Ionization energies of 92 elements in eV for 2nd outermost electron
      double E20[]={0.,54.403,75.62,18.21,25.15,24.376,29.59,35.11,34.98,41.07,
     47.290,15.03,18.823,16.34,19.72,23.405,23.798,27.62,31.81,11.868,
     12.891,13.63,14.205,16.493,15.636,16.178,17.052,18.15,20.286,17.96,
     20.509,15.93,18.63,21.50,21.60,24.565,27.50,11.027,12.233,13.13,
     14.316,16.15,15.26,16.76,18.07,19.42,21.48,16.904,18.86,14.63,
     16.50,18.60,19.09,21.20,25.10,10.001,11.060,10.850,10.550,10.730,
     10.899,11.069,11.241,12.090,11.519,11.670,11.800,11.930,12.050,
     12.184,13.900,14.900,16.2,17.7,16.60,17.00,20.00,18.56,20.50,18.75,
     20.42,15.03,16.68,19.,20.,20.,22.,10.144,12.1,12.,12.,12.};

//    Excitation energies of H for 9 levels in eV 
     double EH0[]={0.000,82258.214,97491.217,102822.768,105290.515,106631.019,
       107439.301,107963.906,108323.575,108580.843,108771.193,108915.969,109028.639,
       109118.040,109190.163,109249.191,109298.112,109339.108,109373.803,109677.619};

//    statistical weight g for H
     double gH0[]={2.0,8.0,18.0,32.0,50.0,72.0,98.0,128.0,162.00,200.00,242.00,
       288.00,338.00,392.00,450.00,512.00,578.00,648.00,722.00,1.00};
      
atom::atom(int na, double* eps0){

     ATM = new string [na];        
     ATM2 = new string [na];        

     Z = new int [na];

     W = new double [na];        
     mass = new double [na];
     perg = new double [na];
     abu = new double [na];
     chi1 = new double [na];        
     chi2 = new double [na];        

     uu1 = new double [na];
     uu3 = new double [na];
     uu2 = new double [na];
     
     EH = new double [nlevelH];
     gH = new double [nlevelH];

     /* Sort abundances */
     int * na2 = sort(eps0,92);


     muavg=0.0;
     abutot=0.0;
     for (int i=0;i<na;i++){
       int j = na2[i];

       ATM[i]=ATM0[j];
       ATM2[i]=ATM20[j];
       Z[i] = Z0[j];
       W[i] = W0[j];
       chi1[i] = E10[j];
       chi2[i] = E20[j];
       
       abu[i]=pow(10.,(eps0[i]-12.));
       mass[i] = W[i]*amu;
       abutot+=abu[i];
       muavg+=abu[i]*W[i];
     }
     
     muavg/=abutot;


     double invmuav = 1.0/(abutot*muavg*amu);
     for (int i=0;i<na;i++)
     {
       perg[i] = abu[i]*invmuav;
     }

     for (int i=0;i<na;i++)
       fprintf(stdout,"atom %d: Z = %d, abundance = %e, perg %e \n",i,Z[i], abu[i],perg[i]);

     fprintf(stdout, "-------------------------------\n");
     fprintf(stdout, "Mean mu = %e \n", muavg);

     for(int i=0;i<nlevelH;i++)
     {
       EH[i]=EH0[i]*1.2398e-4;
       gH[i]=gH0[i];
     }
}
atom::~atom(){
fprintf(stderr,"atom destructor called\n");
delete [] ATM;
delete [] ATM2;
delete [] W;
delete [] chi1;
delete [] chi2;
delete [] EH;
delete [] gH;
}
