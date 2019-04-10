#ifndef __IONS_WITTMANN__   // __IONS_WITTMANN__
#define __IONS_WITTMANN__

#include "atom.h"

      int const ncontr=27;
      int const nmol=2;
      
      double Hm_1,Hp_1,Hn_1,H2p_1,H2n_1,mI_1[ncontr];
      double phtot;

      double f_ji[ncontr][2];
      double *cmol, *gmol, *dcmol;

      double pe_pg10(const atom & at,double,double,double);
      double saha(double,double,double,double,double);
      double acota(double,double,double);
      double acotasig(double,double,double);
      double* leeabun(int);
      double * Heps(double tt);
      double * invert_eos_newton(double p0,double s0,double rho, double eps,double pmin,double pmax,double smin,double smax,double emin,double emax,double rmin,double rmax);
      double bilinear(int n1, int n2,double *x,double *y,float **fxy,double xx,double yy);

      double T_interp(double ee, double dd);
      double p_interp(double ee, double dd);
      double s_interp(double ee, double dd);
      double ne_interp(double ee, double dd);
      double rhoi_interp(double ee, double dd);
      double amb_interp(double ee, double dd);
      double d3_interp(double pp, double ss);
      double eps3_interp(double pp, double ss);

#endif                // __IONS_WITTMANN__

