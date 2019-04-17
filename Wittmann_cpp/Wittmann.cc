#include<stdio.h>
#include<string.h>
#include<math.h>
#include<iostream>
#include <stdlib.h>
#include <pwd.h>
#include <unistd.h>
#include <stdarg.h>
#include <errno.h>
#include <time.h>
#include <ctype.h>
#include <ctime>
#include <cstdlib>
#include <cfloat>

#include "Wittmann.h"
#include "atom.h"
       
const double pe_tol = 1.0e-8;
const double search_tol = 1.0e-8;
const int itermax = 1e7;

const double eos_gamma = 1.65;
const double eos_mu = 0.62;
const double eos_mu_n = 1.3;
const double Rgas=8.314e7;
const double c_temp= (eos_gamma-1.)*eos_mu/Rgas;
const double c_pres= (eos_gamma-1.);
const double mh = 1.67262158e-24;
static int debug_verbose = 0;

float ** ttbl, ** ptbl, ** netbl, ** rhoitbl, ** ambtbl, ** stbl, ** rhotbl, ** epstbl;

double * eps_grid, * r_grid, *s_grid, *p_grid;
int Ns,Nr,Neps,Np;

double delta_eps,delta_r,delta_s,delta_p;
double inv_delta_r,inv_delta_eps,inv_delta_s,inv_delta_p;

/* internal energy offset from MURaM EOS, if required, will need to calculate more accurately */
const double eps_off = 0.0; // 0.8e11;
const double ss_off = 0.0; // -2.08e9;

float** array_2d_contiguous(int nx,int ny){
  float **array = new float*[nx];

  array[0] = new float[nx*ny];

  for(int i=0;i<nx;i++)
    array[i] = array[0]+i*ny;

  return array;
}

using namespace std;
  int main(){
    /* Calls Wittman routines, provided from Borrero, then down the chain:
    Robert Cameron (cameron@mps.mpg.de)
    L. S. Anusha (bhasari@mps.mpg.de)
    D.P. Version to make 3D EOS for Hion MURaM (przybylski@mps.mpg.de) */
        
    /* Minimum and Maximum T */
    double T_min, T_max;
    
    double pg,rhon,rhoi,T,E, pe,ptot;
    int nlevelH=9;
    int ilevelH,iexi;
    double nH_jk[nlevelH];
    double nI_1,nI_2,xi_1,xi_2;
    double E_min, E_max,pe_min,pe_max,pe_a,pe_b;
    double E_search, E_old, E_a, E_b;
    double T_search, T_old, T_a, T_b;
    double xnonhtot,xntot,xnhtot,xnhtot1;

    /* density grid */
    Nr = 3000;
    r_grid = new double [Nr];
    double rmin=log(1.0e-20);
    double rmax=log(2.5e-3);
    
    /* energy grid */
    Neps = 2000;
    eps_grid = new double [Neps];
    double emin=log(1.0e12);
    double emax=log(2.0e15);
    
    /* inverse grids */
    Np = 300;
    p_grid = new double [Np];
    /* Inverse grids only used for bottom boundary, so low pressures not needed */
    double pmin = 4.0;
    double pmax = 20.0;

    Ns = 500;
    s_grid = new double [Ns];
    double smin = 5e8;
    double smax = 3.0e8;

    double n_i[ncontr];
    double mntot;
    double xy1,sumi,ski,sumexi,sumnexi;
    
    ttbl = array_2d_contiguous(Nr,Neps);
    ptbl = array_2d_contiguous(Nr,Neps);
    netbl = array_2d_contiguous(Nr,Neps);
    rhoitbl = array_2d_contiguous(Nr,Neps);
    ambtbl = array_2d_contiguous(Nr,Neps);
    stbl = array_2d_contiguous(Nr,Neps);
    rhotbl = array_2d_contiguous(Ns,Np);
    epstbl = array_2d_contiguous(Ns,Np);


    /* create density grid */

    r_grid[0]=rmin; 
    delta_r=(rmax-rmin)/(Nr-1);
    inv_delta_r = 1.0/delta_r;

    for(int ir=1;ir<Nr-1;ir++){
      r_grid[ir]=r_grid[ir-1]+delta_r;
    }
    r_grid[Nr-1]=rmax;
    
    /* create energy grid */

    eps_grid[0]=emin;
    delta_eps=(emax-emin)/(Neps-1);
    inv_delta_eps = 1.0/delta_eps;

    for(int ieps=1;ieps<Neps-1;ieps++){
      eps_grid[ieps]=eps_grid[ieps-1] + delta_eps;
    }
    eps_grid[Neps-1]=emax;
    
    /* We read the abundances */
    double *eps0=leeabun(92);
    
    /* create atom object */
    atom at(ncontr,eps0);
    
    delete [] eps0;

    int non_converge_count=0;
    int t_out_of_range_count=0;

    for(int ir=0;ir<Nr;ir++){
      double rh=exp(r_grid[ir]);

      /* computer Ntot for all included elements */
      xntot=0.0;
      for (int i=0;i<ncontr;i++)
      {
        n_i[i]=rh*at.perg[i];
        xntot+=n_i[i];
      }
      xnhtot = n_i[0];
      xnonhtot = xntot-xnhtot;

      for(int ieps=0;ieps<Neps;ieps++){
        double ep = exp(eps_grid[ieps]);
      
        fprintf(stdout,"-------------------------------------\n");
        fprintf(stdout,"ir=%d rho=%e | ieps=%d eps=%e \n",ir,rh,ieps,ep);

        double factor = 0.2;
        int converged = 0;
        int minit=0;
        while ((!converged)&&(minit<1000))
        { 
          if (ieps > 1)
            T_min = exp(ttbl[ir][ieps-1]) - (1+factor)*fabs(exp(ttbl[ir][ieps-1])-exp(ttbl[ir][ieps-2]));
          else if (ir > 1)
            T_min = exp(ttbl[ir-1][ieps]) - (1+factor)*fabs(exp(ttbl[ir-1][ieps])-exp(ttbl[ir-2][ieps]));
          else if (ieps > 0)
            T_min=0.95*T_search;
          else if (ir > 0)
            T_min = exp(ttbl[ir-1][ieps])*0.95;
          else
            T_min = 1.0e3;

          T_min = max(1.0e3,T_min);

          T=T_min;
          T_a=T_min;
          fprintf(stdout,"------------  Tmin, Temp=%e -------------------------\n",T);

          for(int i=0;i<ncontr;i++)
            at.partf(i,T);

          /* initial guess for pe */
          if (ieps > 1)
            pe = exp(netbl[ir][ieps-1]) - (1+factor)*fabs(exp(netbl[ir][ieps-1])-exp(netbl[ir][ieps-2]));
          else if (ir > 1)
            pe = exp(netbl[ir-1][ieps]) - (1+factor)*fabs(exp(netbl[ir-1][ieps])-exp(netbl[ir-2][ieps]));
          else 
            pe = xntot;

          pe*=k*T;

          double err_pe = 1.0e8;
          int pe_iter=0;
          int max_pe_iter = 1e6;

          while ((err_pe > pe_tol)&&(pe_iter<max_pe_iter))
          {
            double pe1=pe_pg10(at, T,pe,xnhtot);

            xnhtot1=Hm_1+Hp_1+Hn_1+2.0*H2n_1+2.0*H2p_1;

            pe = pe/xnhtot1*xnhtot;
            Hm_1 = Hm_1/xnhtot1*xnhtot;
            Hp_1 = Hp_1/xnhtot1*xnhtot;
            Hn_1 = Hn_1/xnhtot1*xnhtot;
            H2n_1 = H2n_1/xnhtot1*xnhtot;
            H2p_1 = H2p_1/xnhtot1*xnhtot;

            err_pe = fabs(pe-pe1)/pe;
            pe = pe1;
            pe_iter+=1;
          }

          pe_a = pe;
          pe_min = pe;

          fprintf(stdout,"Hn_1=%e,Hp_1=%e,Hm_1=%e H2n_1=%e H2p_1=%e ne=%e| xnhtot=%e xnhtot1=%e \n",Hn_1,Hp_1,Hm_1,H2n_1,H2p_1,pe*k*T,xnhtot, xnhtot1);

          sumexi=0.;
          /* except the last level */
          for(iexi=0;iexi<nlevelH-1;iexi++){
            nH_jk[iexi]=
              (Hn_1)*at.gH[iexi]*exp(-at.EH[iexi]*ev/(k*T))/at.uu1[0];
            sumexi+=nH_jk[iexi]*at.EH[iexi]*ev;
          }

          /* molecular H2, H2+ and H- energies */
          double * Hmol_eps = Heps(T);
          sumi = Hmol_eps[0]*H2n_1 + Hmol_eps[1]*H2p_1 + Hmol_eps[2]*Hm_1 + Hmol_eps[3]*Hn_1 + Hmol_eps[4]*Hp_1;

          /* Excitation energies */
          for(int i=1;i<ncontr;i++){
            double c = n_i[i]/(1.0+f_ji[i][0]+f_ji[i][1]*f_ji[i][0]);
            nI_1=f_ji[i][0]*c;
            nI_2=nI_1*f_ji[i][1];
            /* First ionisation energy */
            sumi+=(nI_1*at.chi1[i]+nI_2*(at.chi1[i]+at.chi2[i]))*ev;
          }

          ski = 1.5*((Hm_1+Hp_1+Hn_1+H2n_1+H2p_1+xnonhtot)*k*T+pe);

          E_min=ski + sumi+sumexi;
          E_a=ep-E_min/rh;

          if (E_a < 0)
            factor*=2;
          else if ( E_a!=E_a)
            factor*=0.9;
          else
            converged=1;

          minit+=1;
        }
        
        factor = 0.5;
        converged = 0;
        int maxit = 0;
        while ((!converged)&&(maxit<1000))
        {
          if (ieps > 1)
            T_max = exp(ttbl[ir][ieps-1]) + (1.0+factor)*(exp(ttbl[ir][ieps-1])-exp(ttbl[ir][ieps-2]));
          else if (ir > 1)
            T_max = exp(ttbl[ir-1][ieps]) + (1.0+factor)*(exp(ttbl[ir-1][ieps])-exp(ttbl[ir-2][ieps]));
          else if (ieps > 0)
            T_max=1.5*T_search;
          else if (ir > 0)
            T_max = exp(ttbl[ir-1][ieps])*1.25;
          else
            T_max = 1.0e7;

          T_max = min(1.0e7,T_max);

          T=T_max;
          T_b=T_max;
          fprintf(stdout,"------------  Tmax, Temp=%e -------------------------\n",T);

          for(int i=0;i<ncontr;i++)
            at.partf(i,T);

          /* initial guess for pe */
          if (ieps > 1)
            pe = exp(netbl[ir][ieps-1]) + (1+factor)*(exp(netbl[ir][ieps-1])-exp(netbl[ir][ieps-2]));
          else if (ir > 1)
            pe = exp(netbl[ir-1][ieps]) + (1+factor)*(exp(netbl[ir-1][ieps])-exp(netbl[ir-2][ieps]));
          else 
            pe = xntot;

          pe*=k*T;

          double err_pe = 1.0e8;
          int pe_iter=0;
          int max_pe_iter = 1e6;

          while ((err_pe > pe_tol)&&(pe_iter<max_pe_iter))
          {
            double pe1=pe_pg10(at, T,pe,xnhtot);

            xnhtot1=Hm_1+Hp_1+Hn_1+2.0*H2n_1+2.0*H2p_1;

            pe = pe/xnhtot1*xnhtot;
            Hm_1 = Hm_1/xnhtot1*xnhtot;
            Hp_1 = Hp_1/xnhtot1*xnhtot;
            Hn_1 = Hn_1/xnhtot1*xnhtot;
            H2n_1 = H2n_1/xnhtot1*xnhtot;
            H2p_1 = H2p_1/xnhtot1*xnhtot;

            err_pe = fabs(pe-pe1)/pe;
            pe = pe1;
            pe_iter+=1;
          }

          pe_b = pe;
          pe_max = pe;

          xnhtot1=Hm_1+Hp_1+Hn_1+2.0*H2n_1+2.0*H2p_1;
          fprintf(stdout,"Hn_1=%e,Hp_1=%e,Hm_1=%e H2n_1=%e H2p_1=%e ne=%e| xnhtot=%e xnhtot1=%e \n",Hn_1,Hp_1,Hm_1,H2n_1,H2p_1,pe*k*T,xnhtot, xnhtot1);

          sumexi=0.;
          /* except the last level */
          for(iexi=0;iexi<nlevelH-1;iexi++){
            nH_jk[iexi]=
              (Hn_1)*at.gH[iexi]*exp(-at.EH[iexi]*ev/(k*T))/at.uu1[0];
            sumexi+=nH_jk[iexi]*at.EH[iexi]*ev;
          }

          /* molecular H2, H2+ and H- energies  */
          double *Hmol_eps = Heps(T);
          sumi = Hmol_eps[0]*H2n_1 + Hmol_eps[1]*H2p_1 + Hmol_eps[2]*Hm_1 + Hmol_eps[3]*Hn_1 + Hmol_eps[4]*Hp_1;

          /* Excitation energies */
          for(int i=1;i<ncontr;i++){
            double c = n_i[i]/(1.0+f_ji[i][0]+f_ji[i][1]*f_ji[i][0]);
            nI_1=f_ji[i][0]*c;
            nI_2=nI_1*f_ji[i][1];
            /* First ionisation energy */
            sumi+=(nI_1*at.chi1[i]+nI_2*(at.chi1[i]+at.chi2[i]))*ev;
          }

          ski = 1.5*((Hm_1+Hp_1+Hn_1+H2n_1+H2p_1+xnonhtot)*k*T+pe);

          E_max=ski + sumi +sumexi;
          E_b=ep-E_max/rh;

          if (E_b > 0)
            factor*=2;
          else if (E_b!=E_b)
            factor*=0.9;
          else 
            converged=1;

          maxit+=1;
        }

        T = -E_a*(T_b-T_a)/(E_b-E_a)+T_a;

        fprintf(stdout,"------------  Initial guess Temp=%e -------------------------\n",T);
        fprintf(stdout,"E_min=%e eps_min=%e T_min=%e pe_min=%e\n",E_min,E_min/rh,T_min,pe_a);
        fprintf(stdout,"E_max=%e eps_max=%e T_max=%e pe_max=%e\n",E_max,E_max/rh,T_max,pe_b);

        /* Initialize */
        E=E_min; 

        fprintf(stdout,"------------  Begin iterations -------------------------\n");
        if((E_a > 0. && E_b > 0.)||(E_a < 0. && E_b < 0.))
        {
          fprintf(stdout,"E_a and E_b out of bounds\n");
          fprintf(stdout,"Use ideal gas for T, pg and PV=NRT for electron number\n");
          t_out_of_range_count+=1;
          T = c_temp*ep;
          T_search = T;
          pg = c_pres*ep/rh;

          /* pV=nRT */
          pe = pg-rh/(eos_mu_n*mh)*T*k; 
        }
        else 
        {
          /* Temperature Search */

          double rel_err=1.e+19;
          int iterTn=0;
          int itermax=100000;

          while(rel_err > search_tol && iterTn < itermax)
          {
            iterTn=iterTn+1;

            for(int i=0;i<ncontr;i++)
              at.partf(i,T);

            /* initial guess for pe */
            pe = -E_a*(pe_b-pe_a)/(E_b-E_a)+pe_a;
        
            double err_pe = 1.0e8;
            int pe_iter = 0;
            int max_pe_iter = 1e3;
            while ((err_pe > pe_tol)&&(pe_iter<max_pe_iter))
            {
              double pe1=pe_pg10(at, T,pe,xnhtot);

              xnhtot1=Hm_1+Hp_1+Hn_1+2.0*H2n_1+2.0*H2p_1;

              pe = pe/xnhtot1*xnhtot;
              Hm_1 = Hm_1/xnhtot1*xnhtot;
              Hp_1 = Hp_1/xnhtot1*xnhtot;
              Hn_1 = Hn_1/xnhtot1*xnhtot;
              H2n_1 = H2n_1/xnhtot1*xnhtot;
              H2p_1 = H2p_1/xnhtot1*xnhtot;

              err_pe = fabs(pe-pe1)/pe;
              pe = pe1;
              pe_iter+=1;
            }
              
            xnhtot1=Hm_1+Hp_1+Hn_1+2.0*H2n_1+2.0*H2p_1;

            ptot=(Hm_1+Hp_1+Hn_1+H2n_1+H2p_1+xnonhtot)*k*T;
            pg = pe + ptot;

            sumexi=0.;
            sumnexi=0.;
            for(iexi=0;iexi<nlevelH-1;iexi++){ ////except the last level 
              nH_jk[iexi]=
                (Hn_1)*at.gH[iexi]*exp(-at.EH[iexi]*ev/(k*T))/at.uu1[0];
              sumexi+=nH_jk[iexi]*at.EH[iexi]*ev;
              sumnexi+=nH_jk[iexi];
            }
        
            /* Ion density */
            rhoi=(Hp_1+2*H2p_1+Hm_1)*at.W[0];
            rhon=(Hn_1+2*H2n_1)*at.W[0];

            /* molecular H2, H2+ and H- energies */
            double *Hmol_eps = Heps(T);
            sumi = Hmol_eps[0]*H2n_1 + Hmol_eps[1]*H2p_1 + Hmol_eps[2]*Hm_1 + Hmol_eps[3]*Hn_1 + Hmol_eps[4]*Hp_1;

            /* Excitation energies */
            for(int i=1;i<ncontr;i++){
              
              double nN_0 = n_i[i]/(1.0+f_ji[i][0]+f_ji[i][1]*f_ji[i][0]);

              nI_1=f_ji[i][0]*nN_0;
              nI_2=nI_1*f_ji[i][1];

              /* First ionisation energy */
              sumi+=(nI_1*at.chi1[i]+nI_2*(at.chi1[i]+at.chi2[i]))*ev;

              /* Ion density */
              rhoi+=(nI_1+nI_2)*at.W[i];
              rhon+=(nN_0)*at.W[i];
            }
            
            rhoi*=amu;
            rhon*=amu;

            ski=1.5*pg;

            E_old=E;
            E=ski+sumi+sumexi; /* Actual energy */

            /* Difference in energy */
            E_search=ep-E/rh;          

            T_old = T;
            int sw = (int) (E_search*E_a < 0.0);
            T_b = sw*T + (1-sw)*T_b;
            E_b = sw*E_search + (1-sw)*E_b; 
            pe_b = sw*pe + (1-sw)*pe_b;

            int sw2 = (int) (E_search*E_b < 0.0);
            T_a = sw2*T + (1-sw2)*T_a;
            E_a = sw2*E_search + (1-sw2)*E_a;
            pe_a = sw2*pe + (1-sw2)*pe_a;

            //T_search = 0.5*(T_a+T_b);
            T_search = -E_a*(T_b-T_a)/(E_b-E_a)+T_a;
            rel_err=fabs((T_search-T_old)/T_old);
            T=T_search;

            if (debug_verbose==1) {
              fprintf(stdout,"sumi=%e sumexi=%e ski=%e ptot=%e pe=%e\n",sumi,sumexi,ski,ptot,pe);
              fprintf(stdout,"T_a=%e T_search=%e T_b=%e | E_a=%e E_search=%e E_b=%e\n",T_a,T_search,T_b,E_a,E_search,E_b);
              fprintf(stdout,"T_old=%e T_new=%e | eps_old=%e eps_new=%e\n",T_old,T,E_old/rh,E/rh);
              fprintf(stdout,"iterTn=%d rel_err=%e\n",iterTn,rel_err);
              fprintf(stdout,"-------------------------------------------------\n");
            }
          }
          if(iterTn==itermax){
            fprintf(stdout,"Bisection method did not find a root\n");
            fprintf(stdout,"T=%e Hp_1=%e Hn_1=%e\n",T,Hp_1,Hn_1);
            fprintf(stdout,"H2p_1=%e H2n_1=%e\n",H2p_1,H2n_1);
            non_converge_count+=1;
          }
        }
        ttbl[ir][ieps]= (float) log(T);

        ptbl[ir][ieps]= (float) log(pg);
        rhoitbl[ir][ieps]=(float) log(rhoi);
        ambtbl[ir][ieps]=(float) 0.0;
        netbl[ir][ieps]=(float) log(pe/k/T);
        

        fprintf(stdout,"Hp_1=%e Hn_1=%e Hm_1=%e H2p_1=%e H2n_1=%e| xnhtot=%e xnhtot1=%e \n",Hp_1,Hn_1,Hm_1,H2p_1,H2n_1,xnhtot,xnhtot1);
        fprintf(stdout,"T=%lf pg=%e Rho=%e rhoi=%e rhon=%e eps=%e pe =%e ne=%e\n",T,pg,rh,rhoi, rhon,ep,pe,pe/k/T);
        fprintf(stdout,"End iteration\n");
        fprintf(stdout," \n");
      }
    }
      
    for (int ir=0;ir<Nr;ir++)
    {
      double ss1=0.0;
      stbl[ir][0] = 0.0;

      double hon2 = delta_eps/2.0;
      double hon3 = delta_eps/3.0;
      double h3on8 = 3.0*delta_eps/8.0;

      stbl[ir][1] = hon2*(exp(eps_grid[0])/exp(ttbl[ir][0])+exp(eps_grid[1])/exp(ttbl[ir][1]));

      ss1=hon3*(exp(eps_grid[0])/exp(ttbl[ir][0])+4.0*exp(eps_grid[1])/exp(ttbl[ir][1]));

      for (int ieps=2;ieps<Neps;ieps+=2)
      {
        double ss2=ss1+hon3*exp(eps_grid[ieps])/exp(ttbl[ir][ieps]);
        stbl[ir][ieps]=ss2;
        ss2+=hon2*(exp(eps_grid[ieps])/exp(ttbl[ir][ieps])+
            exp(eps_grid[ieps+1])/exp(ttbl[ir][ieps+1]));
        stbl[ir][ieps+1]=ss2;
        ss1+=hon3*(2.0*exp(eps_grid[ieps])/exp(ttbl[ir][ieps])+4.0*exp(eps_grid[ieps+1])/exp(ttbl[ir][ieps+1]));
      }
    }

    double ss1=0.0;
    for (int ieps=0;ieps<Neps;ieps++)
    {
      stbl[0][0]+=ss1;
    }

    double hon2 = delta_r/2.0;
    double hon3 = delta_r/3.0;
    double h3on8 = 3.0*delta_r/8.0;

    double ss2 = -hon2*(exp(ptbl[0][0])/(exp(ttbl[0][0])*exp(r_grid[0]))+exp(ptbl[1][0])/(exp(ttbl[1][0])*exp(r_grid[1])));

    for (int ieps=0;ieps<Neps;ieps++)
    {
      stbl[1][ieps]+=ss2;
    }

    ss1-=hon3*(exp(ptbl[0][0])/(exp(ttbl[0][0])*exp(r_grid[0]))+4.0*exp(ptbl[1][0])/(exp(ttbl[1][0])*exp(r_grid[1])));

    for (int ir=2;ir<Nr;ir+=2)
    {
      ss2=ss1-hon3*exp(ptbl[ir][0])/(exp(ttbl[ir][0])*exp(r_grid[ir]));
      double ss3=ss2-hon2*(exp(ptbl[ir][0])/(exp(ttbl[ir][0])*exp(r_grid[ir]))+exp(ptbl[ir+1][0])/(exp(ttbl[ir+1][0])*exp(r_grid[ir+1])));
      ss1-=hon3*(2.0*exp(ptbl[ir][0])/(exp(ttbl[ir][0])*exp(r_grid[ir]))+4.0*exp(ptbl[ir+1][0])/(exp(ttbl[ir+1][0])*exp(r_grid[ir+1])));
      for (int ieps=0;ieps<Neps;ieps++)
      {
        stbl[ir][ieps]+=ss2;
        stbl[ir+1][ieps]+=ss3;
      }
    }

    fprintf(stdout,"non converge count = %i | t out of range count %i \n",non_converge_count, t_out_of_range_count);
    fprintf(stdout, "setting up inverse tables, smin=%e, smax=%e, pmin=%e,pmax=%e",smin,smax,pmin,pmax);


    p_grid[0]=pmin;
    delta_p=(pmax-pmin)/(Np-1);
    inv_delta_p=1.0/delta_p;

    for(int ip=1;ip<Np-1;ip++){
      p_grid[ip]=p_grid[ip-1]+delta_p;
    }
    p_grid[Np-1]=pmax;

    s_grid[0]=smin;
    delta_s=(smax-smin)/(Ns-1);
    inv_delta_s = 1.0/delta_s;
    for(int is=1;is<Ns-1;is++){
      s_grid[is]=s_grid[is-1]+delta_s;
    }
    s_grid[Ns-1]=smax;
    
    /* Invert for rho and eps tables */
    /* NB Inverse tables are only currently used for the lower boundary conditions.
     * As such we set the min/max rho/eps based on what is expected between a range of depths from -1 Mm
     * the tables currently do not work to a high enough density to do the the 20 Mm deep boxes. Around
     * 3.0e-3 would be needed */


    double ee0 = 2.5e12;
    double ee1 = 1.0e14;
    double ss0 = smin;
    ss1 = smax;
    double rr0 = 1.0e-10;
    double rr1 = 2.5e-3;
    double pp0 = exp(pmin);
    double pp1 = exp(pmax);

    fprintf(stdout, "Starting Inversion \n ");
    fprintf(stdout, " -------------------------------------- \n ");
      

    double max_err_p = 0.0;
    double max_err_s = 0.0;

    for(int is=0;is<Ns;is++)
    {
      int ip=0;

      double p0 = exp(p_grid[ip]);
      double s0 = s_grid[is];

      double rh = log(1e-6);
      double ep = log(5.0e12);

      double * output = invert_eos_newton(p0,s0,rh,ep,pp0,pp1,ss0,ss1,ee0,ee1,rr0,rr1);

      rhotbl[is][ip] = output[0];
      epstbl[is][ip] = output[1];
      max_err_p = max(max_err_p,output[2]);
      max_err_s = max(max_err_s,output[3]);
    }
    fprintf(stdout, "First iteration, ip= 0, max_err_p %e, max_err_s %e \n",max_err_p,max_err_s);
    /* First iteration */
    for(int ip=1;ip<Np;ip++)
    {
      max_err_p = 0.0;
      max_err_s = 0.0;
      double p0 = exp(p_grid[ip]);
      for(int is=0;is<Ns;is++)
      {
        double s0 = s_grid[is];
        
        double rh = rhotbl[is][ip-1];
        double ep = epstbl[is][ip-1];
      
        double * output = invert_eos_newton(p0,s0,rh,ep,pp0,pp1,ss0,ss1,ee0,ee1,rr0,rr1);

        rhotbl[is][ip] = output[0];
        epstbl[is][ip] = output[1];
        max_err_p = max(max_err_p,output[2]);
        max_err_s = max(max_err_s,output[3]);
      }
    fprintf(stdout, "First iteration, ip= %d, max_err_p %e, max_err_s %e \n",ip,max_err_p,max_err_s);
    }
    fprintf(stdout, "First Iteration Complete Starting Inversion Iteration 2\n ");
    fprintf(stdout, " -------------------------------------- \n ");
    /* Second iteration */
    for(int ip=1;ip<Np-1;ip++)
    {
      max_err_p = 0.0;
      max_err_s = 0.0;
      double p0 = exp(p_grid[ip]);
      for(int is=0;is<Ns;is++)
      {
        double s0 = s_grid[is];

        double rh = (rhotbl[is][ip+1]+rhotbl[is][ip-1])/2.0;
        double ep = (epstbl[is][ip+1]+epstbl[is][ip-1])/2.0;
        
        double * output = invert_eos_newton(p0,s0,rh,ep,pp0,pp1,ss0,ss1,ee0,ee1,rr0,rr1);

        rhotbl[is][ip] = output[0];
        epstbl[is][ip] = output[1];
        max_err_p = max(max_err_p,output[2]);
        max_err_s = max(max_err_s,output[3]);
      }
      fprintf(stdout, "Second iteration, ip= %d, max_err_p %e, max_err_s %e \n",ip,max_err_p,max_err_s);
    }

    FILE *outfp;
    float header[14];

    header[0] = (float) Neps;
    header[1] = (float) Nr;
    header[2] = (float) Np;
    header[3] = (float) Ns;
    header[4] = (float) eps_grid[0];
    header[5] = (float) eps_grid[Neps-1];
    header[6] = (float) r_grid[0];
    header[7] = (float) r_grid[Nr-1];
    header[8] = (float) p_grid[0];
    header[9] = (float) p_grid[Np-1];
    header[10] = (float) s_grid[0];
    header[11] = (float) s_grid[Ns-1];
    header[12] = (float) 0.8e11;
    header[13] = (float) -2.08e9;

    outfp = fopen("./eostest.dat","w");
    fwrite(header,sizeof(float),14,outfp);
    fwrite(&ptbl[0][0],sizeof(float),Neps*Nr,outfp);
    fwrite(&ttbl[0][0],sizeof(float),Neps*Nr,outfp);
    fwrite(&stbl[0][0],sizeof(float),Neps*Nr,outfp);
    fwrite(&netbl[0][0],sizeof(float),Neps*Nr,outfp);
    fwrite(&rhoitbl[0][0],sizeof(float),Neps*Nr,outfp);
    fwrite(&ambtbl[0][0],sizeof(float),Neps*Nr,outfp);
    fwrite(&epstbl[0][0],sizeof(float),Np*Ns,outfp);
    fwrite(&rhotbl[0][0],sizeof(float),Np*Ns,outfp);
    fclose(outfp);
  }

/* PE_PG10 (included in EQUISUBMU) evaluates the electonic pressure from the
 * the gaseous pressure and an estimate of the electronic pressure
 * calculates the electronic pressure from the pg and from an estimate of the pe1 */

double pe_pg10(const atom &at, double t,double pe,double xnhtot)
{
  double g1,g2,g3,g4,g5;
  double f1,f2,f3,f4,f5,fe;

  if(t<500.){
    fprintf(stdout,"pe_pg10: temperature < 500 K; temperature = 500 K\n");
    t=500.;
  }
  double theta=5040./t;

  molecb mol(nmol,theta);
  cmol=new double [nmol]();
  dcmol=new double [nmol]();

  for (int i=0;i<nmol;i++){
    cmol[i]=mol.Y[i+1];
    dcmol[i]=mol.dY[i+1];
  }

  g4=0.;
  g5=0.;

  /* to set min vlaue for pe */
  if((pe<=1.0d-42)||(pe!=pe))
  {
    pe=1.0d-42;
    g4=0.;
    g5=0.;
  } 
  else
  {
    g4 = pe*pow(10.,cmol[0]);
    g5 = pe*pow(10.,cmol[1]);
  }

  double a,b,c,d,e;

  g1=0.;

  /* for all other elements except hydrogen, assume LTE and call saha */
  /* Determine how much of each is ionised */
  for(int i=1;i<ncontr;i++)
  {
    a=saha(theta,at.chi1[i],at.uu1[i],at.uu2[i],pe);
    f_ji[i][0]=a;

    b=saha(theta,at.chi2[i],at.uu2[i],at.uu3[i],pe);
    f_ji[i][1]=b;
    c=at.abu[i]/(1.+a*(1.+b));
    g1+=c*a*(1.+2.*b);
  }

  /* phtot */
  phtot=xnhtot*k*t;

  /* p(h+)/p(h) */
  g2=saha(theta,at.chi1[0],at.uu1[0],at.uu2[0],pe);
  f_ji[0][0]=g2;

  /* p(h)/p(h-) */
  g3 = saha(theta,inz_Hm1,1.0,at.uu1[0],pe);
  
  if (g3 > 0.0)
    g3=1.0/g3;

  /* step 1: LTE f1 and f2:later replace by NE */
  a=1.+g2+g3;
  b=2.*(1.+g2/g5*g4);
  c=g5;
  d=g2-g3;
  e=g2/g5*g4;

  double c1=c*b*b+a*d*b-e*a*a;
  double c2=2.0*a*e-d*b+a*b*g1;
  double c3=-(e+b*g1);

  f1=0.5*c2/c1;
  double sgn=c1>0.0?1.00:-1.00;
  f1=-f1+sgn*sqrtl(f1*f1-c3/c1);
  f2=g2*f1;

  g2 = f2/f1;
  f_ji[0][0]=g2;
  e=g2/g5*g4;

  f3=g3*f1;
  f5=(1.0-a*f1)/b;
  f4=e*f5;
  fe=f2-f3+f4+g1;
 
  if(f5<1.e-4L)
  {
    long double const6=g5/pe*f1*f1;
    long double const7=f2-f3+g1;
    for(int i=0;i<50;i++)
    {
      f5=phtot*const6;
      f4=e*f5;
      fe=const7+f4;
    }
  }


  pe=fe*phtot;

  if(pe <= 0.0)
    pe=1.e-42;
 
  /* Scale by hydrogen pressure to conserve hydrogen number
   * and convert to number densities */

  double pretonum = phtot/(k*t);

  Hm_1=f3*pretonum;
  Hp_1=f2*pretonum;
  Hn_1=f1*pretonum;
  H2p_1=f4*pretonum;
  H2n_1=f5*pretonum;

  pe = pe;

  delete [] cmol;
  delete [] dcmol;

  return pe;
}

double acota(double x,double x0,double x1){
  double xx=x;
  if(x<x0)xx=x0;
  if(x>x1)xx=x1;
  return xx;
}

double acotasig(double x,double x0,double x1){
  double xx=x;
  if(x<0.){
    x=-x;
    xx=acota(x,x0,x1);
    xx=-xx;
  }
  else{
    xx=acota(x,x0,x1);
  }
  return xx;
}

double * Heps(double tt)
{
  double* E = new double [5];

  double theta = 5040.39/tt;

  /* Vardya Polynomial fits to molecular energies */
  double dE_H2 = 2.6757-1.4772*theta+0.60602*theta*theta-0.12427*theta*theta*theta
    +0.0097503*theta*theta*theta*theta;
  dE_H2*=k*tt;

  double dE_H2p = 2.9216-2.0036*theta+1.7231*theta*theta-0.82685*theta*theta*theta
    +0.15253*theta*theta*theta*theta;
  dE_H2p*=k*tt;

  /* H2 */
  E[0] = dE_H2;
  /* H2p */
  E[1] = (dE_H2p+D0_h2-D0_h2p+inz_H)*ev;
  /* Hm */
  E[2] = (0.5*D0_h2-inz_Hm1)*ev;
  /* Hn */
  E[3] = (0.5*D0_h2)*ev;
  /* Hp */
  E[4] = (0.5*D0_h2+inz_H)*ev;

  return E;
}

double * invert_eos_newton(double p0,double s0,double rho, double eps,double pmin,double pmax,double smin,double smax,double emin,double emax,double rmin,double rmax)
{

  double tol = 1.0e-4;
  double err_p = 1.0e10;
  double err_s = 1.0e10;
  double fct = 1.01;
  
  int maxiter = 1e4;
  int it=0;
   
  double rho1 = exp(rho);
  double eps1 = exp(eps);
    
  double p1 = p_interp(eps1,rho1);
  double s1 = s_interp(eps1,rho1);

  while(((err_p > tol)&&(err_s > tol))&&(it<maxiter))
  {
    double p_r = (p_interp(eps1,fct*rho1)-p1)/((fct-1)* rho1);
    double p_e = (p_interp(fct*eps1,rho1)-p1)/((fct-1)* eps1);
    double s_r = (s_interp(eps1,fct*rho1)-s1)/((fct-1)* rho1);
    double s_e = (s_interp(fct*eps1,rho1)-s1)/((fct-1)* eps1);

    double det = p_r*s_e-p_e*s_r;

    double r1 = ( s_e*(p0-p1)-p_e*(s0-s1))/det;
    double e1 = (-s_r*(p0-p1)+p_r*(s0-s1))/det;

    rho1 = min(max(rmin,rho1+r1),rmax);
    eps1 = min(max(emin,eps1+e1),emax);

    p1 = p_interp(eps1,rho1);
    s1 = s_interp(eps1,rho1);

    err_p = fabs(p1-p0)/p0;
    err_s = fabs(s1-s0)/s0;

    it+=1;
  }

  double * output = new double [4];
  output[0] = log(rho1);
  output[1] = log(eps1);
  output[2] = err_p;
  output[3] = err_s;

  return output;
}

double bilinear(int n1, int n2,double *x,double *y,float **fxy,double xx,double yy)
{

  double c00,c01,c10,c11;
  double ff;

  for (int i=0; i<n1-1; i++)
  {
    for (int j=0; j<n2-1; j++)
    {

      if(((xx>=x[i] && xx<x[i+1])||(xx>x[i] && xx<=x[i+1]))&&
          ((yy>=y[j] && yy<y[j+1])||(yy>y[j] && yy<=y[j+1])))
      {
        c00=(x[i+1]-xx)*(y[j+1]-yy)/((x[i+1]-x[i])*(y[j+1]-y[j]));
        c01=(xx-x[i])*(y[j+1]-yy)/((x[i+1]-x[i])*(y[j+1]-y[j]));
        c10=(x[i+1]-xx)*(yy-y[j])/((x[i+1]-x[i])*(y[j+1]-y[j]));
        c11=(xx-x[i])*(yy-y[j])/((x[i+1]-x[i])*(y[j+1]-y[j]));
        ff=c00*fxy[i][j]+c01*fxy[i+1][j]+c10*fxy[i][j+1]+c11*fxy[i+1][j+1];
      }

      else if(xx < x[0] || xx > x[n1-1] || yy < y[0] || yy > y[n2-1])
      {
        fprintf(stdout,"Out of range %e %e %e %e %e %e\n",xx,x[0],x[n1-1],yy,y[0],y[n2-1]);
        ff=0.;
      }
    }
  }
  return ff;
}

double s_interp(double ee, double dd) {

  int i,j;
  double logr, ss,ee1;

   ee1 = log(ee);
   logr = log(dd);

  i = (int) ((ee1-eps_grid[0])*inv_delta_eps );
  if (i < 0)         i=0;
  if (i > Neps - 2) i=Neps-2;

  j = (int) ( (logr-r_grid[0])*inv_delta_r );
  if (j < 0)        j=0;
  if (j > Nr - 2) j=Nr - 2;

  ss =  ((logr-r_grid[j])   * (ee1-eps_grid[i])   * stbl[j+1][i+1]
       +(logr-r_grid[j])   * (eps_grid[i+1]-ee1) * stbl[j+1][i]
       +(r_grid[j+1]-logr) * (ee1-eps_grid[i])   * stbl[j][i+1]
       +(r_grid[j+1]-logr) * (eps_grid[i+1]-ee1) * stbl[j][i])*inv_delta_eps*inv_delta_r;

  return ss = ss;
}

double p_interp(double ee, double dd) {

  int i,j;
  double logr, logp,ee1;

  ee1 = log(ee);
  logr = log(dd);

  i = (int) ((ee1-eps_grid[0])*inv_delta_eps );
  if (i < 0)         i=0;
  if (i > Neps - 2) i=Neps-2;

  j = (int) ( (logr-r_grid[0])*inv_delta_r );
  if (j < 0)        j=0;
  if (j > Nr - 2) j=Nr - 2;

  logp =  ((logr-r_grid[j])   * (ee1-eps_grid[i])   * ptbl[j+1][i+1]
       +(logr-r_grid[j])   * (eps_grid[i+1]-ee1) * ptbl[j+1][i]
       +(r_grid[j+1]-logr) * (ee1-eps_grid[i])   * ptbl[j][i+1]
       +(r_grid[j+1]-logr) * (eps_grid[i+1]-ee1) * ptbl[j][i])*inv_delta_eps*inv_delta_r;

  logp = exp(logp);

  return logp;
}
