#include <errno.h>
#include <string.h>
#include <mpi.h>
#include <time.h>
#include "grid.H"
#include "run.H"
#include "physics.H"
#include "io.h"
#include "mem.h"
#include "flags.h"
#include "cfg.h"
#include "atom.h"
#include "rt.h"
#include "rt_nlte.h"
#include "angle.h"

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <complex>
#include <cmath>

#include "precision.h"
#include "tool.h"
#include <sys/time.h>
#include <sys/resource.h>

#include "rt_nlte.h"
#include "atom.h"
//#define print_ions 1
#define kbl  1.3806503E-16L            // erg/K

void RTS_NLTE::chemical_equilibrium(long double rho_in,long double eps_in)     
     {
//L. S. Anusha (bhasari@mps.mpg.de)
//Calls routines implementing Wittmann provided to cameron by Borrero; to Anusha by Robert
//STOPS AT 3000 degrees K.
#ifdef print_ions
      fprintf(stderr,"In chemical equilibrium: lte_switch=%d\n",lte_switch);
#endif
      long double T,T_prev,N_el_prev;
      long double Temp=0.L;
      int ninzH=2;
      itermax=1000;
//      int nlevelHe=17;
//      int ninzHe=3;
//      int ilevelH,ilevelHe;
      long double E_min, E_max, E_old=0.L, E_current=0.L;
      long double E_search, E_a, E_b;
      long double T_search, T_a, T_b;
//      long double rel_err,crec=1e-20L;
      long double rel_err,crec=1e-8L;
      long double n_i[ncontr];
      int iterTn;
      int iter_en=1;
//////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
       long double TK=1000.L;

            for(int i=0;i<ncontr;i++)n_i[i]=0.L;
           
// We read the abundances
/////////////Construct Pg grids//////////////////////////////////////////
/*
//       long double Tmin=1000.L
       long double Tmin;
//       if(lte_switch==1) Tmin=1000.L;
//       if(lte_switch==0) Tmin=Temp_in-1500.L;
//       if(lte_switch==0) 
          Tmin=3000.L;
       long double Tmax=5.e+8L;
*/
       long double Tmin;
       long double Tmax;
       Tmin=3000.L;
       Tmax=5.e+8L;
//       if(rho_in>=2.e-7) 
//        {
//        Tmin=Temp_in-1000.L;
//        Tmax=Temp_in+1000.L;
//        }



//         fprintf(stderr,"Temp_in=%Le Tmin=%Le\n",Temp_in,Tmin);
//       long double Pg_min=3.5707081952196196e-08;
//     long double Pg_max=4859227.1531350948;
//       long double Pg_max=8.3178167846129988e+21;
////////////////////////////////////////////////////////////////////////
//////////////////create an atom object/////////////////////////////////
//already read in rt_nlte.cc
            chi1=new long double[ncontr];
            chi2=new long double[ncontr];
            for(int i=0;i<ncontr;i++)
            {
            chi1[i]=(long double)at->chi1[i];
            chi2[i]=(long double)at->chi2[i];
//            fprintf(stderr,"i=%d chi1=%Le chi2=%Le\n",i,chi1[i],chi2[i]);
            }
/////////////////////Compute mass of all elements///////////////////
            mass[0]=(long double)(at->W[0]*amu);  ////atomic mass unit
            for(int i=0;i<ncontr;i++)mass[i]=(long double)(at->W[i]*amu);

            long double xmntot=mass[0]; //H
            for(int i=1;i<ncontr;i++)
            {
            xmntot=xmntot+mass[i]*(abu[i]); //mass is the atomic mass
            }
//            for(int i=0;i<ncontr;i++)
//            fprintf(stderr,"mass=%Le abu=%Le\n",mass[i],abu[i]);
 

////////////Compute N_tot for a given rho///////////////////////////////
           long double rh=rho_in;
/////////Compute N_tot for a given rho//////////////////////////////////
//            fprintf(stderr,"xmntot=%Le\n",xmntot);
            n_i[0]=rh/xmntot; ////Get n_H
//            fprintf(stderr,"nH=%Le\n",n_i[0]);
            xnhtot=n_i[0];

            xntot=0.L;
            for (int i=0;i<ncontr;i++)
            {
            n_i[i]=(abu[i])*xnhtot; //// n_H times abundance 
            xntot=xntot+n_i[i];
//            fprintf(stderr,"i=%d n_i=%Le\n",i,n_i[i]);
            }
#ifdef print_ions
            fprintf(stderr,"rho=%Le eps=%Le\n",rho_in,eps_in);
#endif
             
            uu1=new long double[ncontr];
            uu2=new long double[ncontr];
            uu3=new long double[ncontr];
///////////////NR METHOD FOR FINDING T////////////////////////////

///////////////////////////////////////////////////////////////////
////////////////////////////////////T-search////////////////////////
          rel_err=1.e+19;
          iterTn=0;
          T=Temp_in; //initialization
////////INNER ITERATION STARTS///////////////////////////////////////
           do 
            {
            iterTn=iterTn+1;
            T=T>TK?T:TK;
            E_current=compute_energy(T,xnhtot,xntot);
            long double Tu=T+10.;// T perturbation
            long double Eu_current=compute_energy(Tu,xnhtot,xntot);
            long double E_prime_T=-(Eu_current-E_current)/(rh*10.);

            E_search=eps_in-(E_current)/rh; ////Difference in energy
            T_prev=T;
             
            T_search=T_prev-E_search/E_prime_T;

//            rel_err=fabs(eps_in-(E_current/rh))/eps_in;
            rel_err=fabs(T_search-T_prev)/T_search;

//            rel_err=fabs(T_search-T_prev);
            T=T_search; ////Update T

#ifdef print_ions
            fprintf(stderr,"T_prev=%Le T_search=%Le diff=%Le\n",T_prev,T_search,T_search-T_prev);
            fprintf(stderr,"E_current=%Le Eu_current=%Le E_search=%Le E_prime_T=%Le\n",
                            E_current,Eu_current,E_search,E_prime_T);
            fprintf(stderr,"iterTn=%d rel_err=%Le\n",iterTn,rel_err);
#endif

#ifdef print_ions
         fprintf(stderr,"Hp_1=%Le Hn_1=%Le Hm_1=%Le H2p_1=%Le H2n_1=%Le\n",Hp_1,Hn_1,Hm_1,H2p_1,H2n_1);
         fprintf(stderr," \n");
#endif
            }while(rel_err > crec && iterTn < itermax);

            if(iterTn==itermax) 
            {
            fprintf(stderr,"Newton's method did not find a root: lte_switch=%d\n",lte_switch);
            fprintf(stderr,"T=%Le Hp_1=%Le Hn_1=%Le\n",T,Hp_1,Hn_1);
            fprintf(stderr,"H2p_1=%Le H2n_1=%Le\n",H2p_1,H2n_1);
             exit(0);
            }
//////////////////////////////////////////////////////////////
//            T=0.5L*(T_a+T_b);
//            T=0.5L*(T+T_prev);
//////////////////////////////////////////////////////////////////////
            pe=inner_iteration(T,xnhtot,xntot);
            N_el=pe/(kbl*T);
/////////////////////////////////////////////////////////////////////
            Temp=T;
            P_g=(N_el+ntot)*kbl*T; //update when computed
//////////////////output/////////////////////////////////////
            T_out=T; 
            Pg_out=P_g;
            Pe_out=pe;
            Nel_out=N_el;
            nHt=Hp_1+Hn_1; //atomic hydrogen
//////////////////output/////////////////////////////////////
#ifdef print_ions
         fprintf(stderr,"Hp_1=%Le Hn_1=%Le Hm_1=%Le H2p_1=%Le H2n_1=%Le\n",Hp_1,Hn_1,Hm_1,H2p_1,H2n_1);
         fprintf(stderr,"charge cons: Hp_1=%Le Hm_1=%Le H2p_1=%Le Nel=%Le g1*nh=%Le ch-n=%Le ch-p=%Le\n",
                 Hp_1,Hm_1,H2p_1,Nel_out,g1*phtot/(kbl*T),Hp_1+H2p_1+g1*phtot/(kbl*T),Hm_1+Nel_out);
         fprintf(stderr," \n");
        
         fprintf(stderr,"T=%Le Pg=%Le Pe=%Le Nel=%Le nHt=%Le\n",T_out,Pg_out,Pe_out,Nel_out,nHt);
         fprintf(stderr,"P_g=%Le Rho=%Le Eps=%Le\n",P_g,rh,eps_in);
         fprintf(stderr,"End iteration\n");
         fprintf(stderr," \n");
#endif
         delete [] uu1;
         delete [] uu2;
         delete [] uu3;
         delete [] chi1;
         delete [] chi2;
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
       
}

long double RTS_NLTE::compute_energy(long double T,long double xnhtot,long double xntot)
{
//            fprintf(stderr,"Temp=%Le\n",T);
       long double summ,sumn,sumnexi,sumi;
       long double n_i1[ncontr];
       long double nI_1=0.L,nI_2=0.L,xi_1=0.L,xi_2=0.L;
       int nlevelH=nlevel;

            for(int i=0;i<ncontr;i++)
              n_i1[i]=0.L;

            for(int i=0;i<ncontr;i++)
             {
             int ii=i+1;
             float64_t Tdd=(double)T;
             at->partfun(ii,Tdd);
             uu1[i]=(long double)at->u1;
             uu2[i]=(long double)at->u2;
             uu3[i]=(long double)at->u3;
             }
/////////////////////////////////////////////////////////////////////////////
            pe=inner_iteration(T,xnhtot,xntot);
#ifdef print_ions
         fprintf(stderr,"Hp_1=%Le Hn_1=%Le Hm_1=%Le H2p_1=%Le H2n_1=%Le\n",Hp_1,Hn_1,Hm_1,H2p_1,H2n_1);
         fprintf(stderr," \n");
#endif
/////////////////////////////////////////////////////////////////////////////
            N_el=pe/(kbl*T);

            long double xnhtot1=Hm_1+Hp_1+Hn_1+2.0L*H2n_1+2.0L*H2p_1;
//          fprintf(stderr,"Hm_1=%Le H2n_1=%Le H2p_1=%Le\n",Hm_1,H2n_1,H2p_1);
//          fprintf(stderr,"xnhtot=%Le xnhtot1=%Le\n",xnhtot,xnhtot1);
            sumn=0.L;
            summ=0.L;
            for(int i=0;i<ncontr;i++)
            {
            n_i1[i]=(abu[i])*xnhtot1;
            summ=summ+mass[i]*n_i1[i];
            sumn=sumn+n_i1[i];
//            fprintf(stderr,"abu=%Le mass=%Le summ=%Le sumn=%Le n_i1=%Le\n",abu[i],mass[i],summ,sumn,n_i1[i]);
            }
//           fprintf(stderr,"mu_a=%Le\n",summ/sumn); //////nearly 1.29, I get 1.3
            ntot=sumn;

            P_g=(N_el+ntot)*kbl*T; //update when computed

            if(lte_switch==1)
            {
            sumexi=0.L;
//////////with electrons and protons: <nlevel+1
//////////remove electrons: <nlevel
//////////remove protons: <nlevel-1
            for(int iexi=0;iexi<nlevelH-1;iexi++) ////except the last level 
             {
            long double gH=(long double)(at->gH[iexi]);
            long double EH0=(long double)(at->EH[iexi]*ev);
            long double EHkT=((long double)(at->EH[iexi]*ev))/(kbl*T);

            nH_jk[iexi][0]=
            (Hn_1)*gH*expl(-EHkT)/uu1[0];
//            fprintf(stderr,"gH=%Le uu1=%Le exp=%Le\n",gH,uu1[0],expl(-EHkT));
            sumexi=sumexi+nH_jk[iexi][0]*EH0;
             }
            }
            if(lte_switch==0)
            {
            sumexi=0.L;
            for(int iexi=0;iexi<nlevelH-1;iexi++) ////except the last level 
             {
            long double EH0=(long double)(at->EH[iexi]*ev);

            sumexi=sumexi+nH_jk[iexi][0]*EH0;
             }
            }
            long double diff_sumexi=(sumexi-sumexi_old);
///////////////////////////////////////////////////////////
            sumHpi=(Hp_1)*(chi1[0]*(long double)(ev)+0.5L*(long double)((D0_h2)*ev)); ////compute sumi
            long double diff_sumHpi=(sumHpi-sumHpi_old);

/////////H2 formation energy in the neutrals//////////////////
/////////Mihals et al.(1988) table 2//////////////////////////
             sumHni=Hn_1*0.5L*((long double)((D0_h2)*ev));
            long double diff_sumHni=(sumHni-sumHni_old);
//////////////////////////////////////////////////////////////
             
///////////Excitation energies of other elements////////////////////////////////             
            sumexi_o=0.L;
           for(int i=1;i<ncontr;i++)
             {
              xi_1=f_ji[i][0]/(1.0L+f_ji[i][0]+f_ji[i][1]*f_ji[i][0]);
              xi_2=f_ji[i][1]*f_ji[i][0]/(1.0L+f_ji[i][0]+f_ji[i][1]*f_ji[i][0]);
              nI_1=xi_1*n_i1[i];
              nI_2=xi_2*n_i1[i];
           sumexi_o=sumexi_o+(nI_1*chi1[i]+nI_2*chi2[i])*(long double)(ev); ////First ionization energy
              }
            long double diff_sumexi_o=(sumexi_o-sumexi_o_old);

            ski=(3.0L/2.0L)*kbl*T*(N_el)+(3.0L/2.0L)*kbl*T*(ntot);
             long double ldeltat=(long double)deltat;
             long double dNel=(N_el-Nel_old)/ldeltat;
             long double dTemp=(T-Temp_old)/ldeltat;

            long double diff_ski=(ski-ski_old);
             
            sum_en_all=Hm_1*(long double)((0.5*D0_h2*ev-Hm1_ien*ev));
            long double diff_sum_en_all=(sum_en_all-sum_en_all_old); 
           
            long double E_tot=sumHpi
                             +sumHni
                             +sumexi
                             +sumexi_o
                             +ski
                             +sum_en_all;
return E_tot;
}
//////////////////////////////////INNER ITERATION FOR GUSESSED T AND GIVEN RHO////////////////////////////////////
     long double RTS_NLTE::inner_iteration(long double t,long double xnhtot,long double xntot)
      {
//         fprintf(stderr,"lte_switch=%d\n",lte_switch);
        long double pe1,pe2;
//        if(lte_switch==1) 
         pe1=xntot*kbl*t; ////Initialize
//        if(lte_switch==0) pe1=Pe_in; 
/////////////////////////////////////////////////////////////
        pe2=pefrompg10(t,xnhtot,pe1,xntot);
//        fprintf(stderr,"pe1=%Le xntot=%Le k=%Le t=%Le\n",pe1,xntot,k,t);
/////////////////////////////////////////////////////////////

      return pe2;
      }
///////////////////////INNER ITERATION//////////////////////////
////////////////////////////////////////////////////////////////
// PEFROMPG10 (included in EQUISUBMU) evaluates the electonic pressure from the
// temperature and gas pressure
// pefrompg10 evaluates the electonic pressure p corresponding to t1 and pg

      long double RTS_NLTE::pefrompg10(long double t,long double xnhtot,long double p,long double xntot)
      {
        long double prec=1e-8L;
        long double dif=1.L;
        int n2=0;
        long double p1=p;
        long double p2;
        do 
         {
	   n2+=1;
           p=(p+p1)/2.L; //If p is small enough, this step will make pe=0.
           p1=p;

           p2=pe_pg10(t,p,xnhtot);
           p=p2;
           dif=fabsl(fabsl( (p-p1))/(p+p1) );

//         fprintf(stderr,"n2=%d p_new=%Le p_old=%Le dif=%Le\n",n2,p,p1,dif);
         } while (dif>prec && n2<itermax);
          
          if(n2==itermax) 
          {
          fprintf(stderr,"Iteration did not converge\n");
           exit(0);
          }
	return p2;
       }
/////////////////////////////////////////////////////////////////////
//PE_PG10 (included in EQUISUBMU) evaluates the electonic pressure from the
//the gaseous pressure and an estimate of the electronic pressure
//calculates the electronic pressure from the pg and from an estimate of the pe1

      long double RTS_NLTE::pe_pg10(long double t,long double pe1,long double xnhtot)
      {
      long double g2,g3,g4,g5;
      long double f1,f2,f3,f4,f5,fe;
      long double xmin,xmax;
      long double pe2;

      if(t<500.L)
      {
      fprintf(stderr,"pe_pg10: temperature < 500 K; temperature = 500 K\n");
      t=500.L;
      }
      long double theta=5040.L/t;
      float64_t thetad=(double)theta;

        molecb mol(nmol,thetad);

      cmol=new long double [nmol];
      gmol=new long double [nmol];
      dcmol=new long double [nmol];
      long double du0,du1,du2;
//      fprintf(stderr,"nmol=%d \n",nmol);

      for (int i=0;i<nmol;i++)
         {
         cmol[i]=(long double)(mol.Y[i+1]);
         dcmol[i]=(long double)(mol.dY[i+1]);
//       fprintf(stderr,"i=%d Y=%Le dY=%Le\n",i,cmol[i],dcmol[i]);
         }

 
      g4=0.L;
      g5=0.L;
      pe2=pe1;
///     to set min vlaue for pe2
      if(pe2<0.L || fabs(pe2) <= 1.e-15L)
       {
         pe2=1.e-15L;
         g4=0.L;
         g5=0.L;
       }
      else
       {
        for(int i=0;i<2;i++)
        {
          xmin=-30.L;
          xmax=+30.L;
           long double cmol1=cmol[i];
//          gmol[i]=pe2*powl(10.L,(cmol[i]));
        }
         g4=pe2*powl(10.L,(cmol[0]));
         g5=pe2*powl(10.L,(cmol[1]));
//         fprintf(stderr,"g4=%Le g5=%Le\n",g4,g5);
       }

// now I calculate the levels u0, u1, u2 and their derivatives

       

      g1=0.L;
// for all other elements except hydrogen, assume LTE and call saha
       long double a,b,c,d,e,c1,c2,c3;
       for(int i=1;i<ncontr;i++)
       {
        long double chi=chi1[i];
        long double u1i=uu1[i];
        long double u2i=uu2[i];
//         fprintf(stderr,"i=%d chi=%Le u1i=%Le u2i=%Le\n",i,chi,u1i,u2i);

        a=saha(theta,chi,u1i,u2i,pe2);
        f_ji[i][0]=a;

        long double chii=chi2[i];
        long double u3i=uu3[i];
        b=saha(theta,chii,u2i,u3i,pe2);
        f_ji[i][1]=b;

        c=1.L+a*(1.L+b);
        g1=g1+((abu[i]))/c*a*(1.L+2.L*b);
//        fprintf(stderr,"i=%d a=%Le b=%Le c=%Le g1=%Le\n",i,a,b,c,g1);
       }
//       fprintf(stderr,"g1=%Le\n",g1);
///////////COMPUTE phtot from density////////////////////
        phtot=xnhtot*kbl*t;
///////////////////compute g2 for step 1 and g3////////////////////////////
       long double chi0=chi1[0];
       long double u10=uu1[0];
       long double u11=uu2[0];
       if(lte_switch==1)
       {
       g2=saha(theta,chi0,u10,u11,pe2);   // p(h+)/p(h)      
       f_ji[0][0]=g2;
       }
       else
       {
       g2=Hp_1/Hn_1;
       f_ji[0][0]=g2;
       }

       long double cc=Hm1_ien*1.L;
       long double dd=1.L;
       g3=saha(theta,cc,dd,u10,pe2);        // p(h)/p(h-) 
       xmin=1.e-30L;
       xmax=1.e+30L;
       long double g31=g3;
//       g3=acota(g31,xmin,xmax);
       if(g3>0.L) g3=1.0L/g3;
/////////////////////////////////////////////////////////////////
/////////////step 1: LTE f1 and f2:later replace by NE///////////
//        fprintf(stderr,"g2=%Le g3=%Le g4=%Le g5=%Le\n",g2,g3,g4,g5);
      
        a=1.L+g2+g3;
        b=2.L*(1.L+g2/g5*g4);
        c=g5;
        d=g2-g3;
        e=g2/g5*g4;
        xmin=1.e-15L;
        xmax=1.e+15L;
        long double a1=a;
//        a=acotasig(a1,xmin,xmax);
        long double d1=d;
//        d=acotasig(d1,xmin,xmax);
//        fprintf(stderr,"a=%Le b=%Le d=%Le e=%Le\n",a,b,d,e);

        c1=c*b*b+a*d*b-e*a*a;
        c2=2.L*a*e-d*b+a*b*g1;
        c3=-(e+b*g1);
//        fprintf(stderr,"c1=%Le c2=%Le c3=%Le\n",c1,c2,c3);

        if(lte_switch==1)
        {
        f1=0.5L*c2/c1;
        long double sgn=c1>0.L?1.0L:-1.0L;
        f1=-f1+sgn*sqrtl(f1*f1-c3/c1);
        f2=g2*f1;
        Hn_1=f1*xnhtot;
        Hp_1=f2*xnhtot;
//        fprintf(stderr,"f1=%Le t1=%Le t2=%Le sgn=%Le\n",f1,f1*f1-c3/c1,sqrtl(f1*f1-c3/c1),sgn);
        }
//////////////////////////////////////////////////////////////////
     


/////////////////step 2///////////////////////////////////////////
         if(lte_switch==0)
         {
        f1=Hn_1*kbl*t/phtot; //fixed from input
        f2=Hp_1*kbl*t/phtot; //fixed from input
         }
//////////////////////////////////////////////////////////////////

         
////////////STEP 2:My method//////////////////////////////////////
        if(f1 !=0.)
        {
        g2=Hp_1/Hn_1;
        f_ji[0][0]=g2;
        e=g2/g5*g4;
        f3=g3*f1;
        f5=(1.L-a*f1)/b;
        f4=e*f5;
        fe=f2-f3+f4+g1;
        xmin=1.e-30L;
        xmax=1.e+30L;
        long double fe1=fe;
//      fe=acota(fe1,xmin,xmax);

//      fprintf(stderr,"f1=%Le f2=%Le f3=%Le f4=%Le f5=%Le fe=%Le\n",
//                       f1,f2,f3,f4,f5,fe);
        }
        else
        {
          f_ji[0][0]=0.L;
          f5=1.L/b;
          f4=e*f5;
          f3=0.L;
          fe=f2-f3+f4+g1;
          xmin=1.e-30L;
          xmax=1.e+30L;
          long double fe1=fe;
//          fe=acota(fe1,xmin,xmax);
        }
      
///*
         if(f5<1.e-4L) 
         {
          long double const6=g5/pe2*f1*f1;
          long double const7=f2-f3+g1;
          for(int i=0;i<5;i++)
          {
           f5=phtot*const6;
           f4=e*f5;
           fe=const7+f4;
          long double fe1=fe;
//          fe=acota(fe1,xmin,xmax);
          }
         }
//*/
      pe2=fe*phtot;
///  added condition to avoid nan
      if(pe2<.0) pe2=1.e-15L;

//now notation is number density.
        Hm_1=f3*phtot/(kbl*t);
        Hp_1=f2*phtot/(kbl*t);
        Hn_1=f1*phtot/(kbl*t);
        H2p_1=f4*phtot/(kbl*t);
        H2n_1=f5*phtot/(kbl*t);
//       fprintf(stderr,"f1=%Le f2=%Le f3=%Le f4=%Le f5=%Le fe=%Le pe2=%Le T=%Le\n",
//                       f1,f2,f3,f4,f5,fe,pe2,t);
      delete [] cmol;
      delete [] gmol;
      delete [] dcmol;
       
      return pe2;
      }
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

        long double RTS_NLTE::acota(long double x,long double x0,long double x1)
        {
//      restricts the input between x0 ans x1
        long double xx=x;
        if(x<x0)xx=x0;
        if(x>x1)xx=x1;
   
        return xx;
        }

        long double RTS_NLTE::acotasig(long double x,long double x0,long double x1)
        {
        long double xx=x;
        if(x<0.)
         {
           x=-x;
           xx=acota(x,x0,x1);
           xx=-xx;
         }
         else
         {
           xx=acota(x,x0,x1);
         }
        return xx;
        }


void RTS_NLTE::chemical_equilibrium_bisection(long double rho_in,long double eps_in)     
     {
//L. S. Anusha (bhasari@mps.mpg.de)
//Calls routines implementing Wittmann provided to cameron by Borrero; to Anusha by Robert
//STOPS AT 3000 degrees K.
#ifdef print_ions
      fprintf(stderr,"In chemical equilibrium: lte_switch=%d\n",lte_switch);
#endif
      long double T,T_prev,N_el_prev;
      long double Temp=0.L;
      int ninzH=2;
      itermax=10000;
//      int nlevelHe=17;
//      int ninzHe=3;
//      int ilevelH,ilevelHe;
      long double E_min, E_max, E_old=0.L, E_current=0.L;
      long double E_search, E_a, E_b;
      long double T_search, T_a, T_b;
      long double rel_err,crec=1e-20L;
//      long double rel_err,crec=1e-8L;
      long double n_i[ncontr];
      int iterTn;
      int iter_en=1;
//////////////////////////////////////////////////////////
            sumexi_old=sumexi;
            sumHpi_old=sumHpi;
            sumHni_old=sumHni;
            sumexi_o_old=sumexi_o;
            ski_old=ski;
            sum_en_all_old=sum_en_all;
////////////////////////////////////////////////////////////
       long double TK=1000.L;

            for(int i=0;i<ncontr;i++)n_i[i]=0.L;
           
// We read the abundances
/////////////Construct Pg grids//////////////////////////////////////////
/*
//       long double Tmin=1000.L
       long double Tmin;
//       if(lte_switch==1) Tmin=1000.L;
//       if(lte_switch==0) Tmin=Temp_in-1500.L;
//       if(lte_switch==0) 
          Tmin=3000.L;
       long double Tmax=5.e+8L;
*/
       long double Tmin;
       long double Tmax;
       Tmin=3000.L;
       Tmax=5.e+8L;

//Whatever temperature we use, output temperature differs slightly from Temp_in.
//Therefore when lte_switch=1 its confusing if we use this 
       if(lte_switch==0)
       {
       if(rho_in>=2.e-7) 
//       if(rho_in>=1.e-7) 
        {
        Tmin=Temp_in-1000.L;
        Tmax=Temp_in+1000.L;
        }
       }


//         fprintf(stderr,"Temp_in=%Le Tmin=%Le\n",Temp_in,Tmin);
//       long double Pg_min=3.5707081952196196e-08;
//     long double Pg_max=4859227.1531350948;
//       long double Pg_max=8.3178167846129988e+21;
////////////////////////////////////////////////////////////////////////
//////////////////create an atom object/////////////////////////////////
//already read in rt_nlte.cc
            chi1=new long double[ncontr];
            chi2=new long double[ncontr];
            for(int i=0;i<ncontr;i++)
            {
            chi1[i]=(long double)(at->chi1[i]);
            chi2[i]=(long double)(at->chi2[i]);
//            fprintf(stderr,"i=%d chi1=%Le chi2=%Le\n",i,chi1[i],chi2[i]);
            }
/////////////////////Compute mass of all elements///////////////////
            mass[0]=(long double)(at->W[0]*amu);  ////atomic mass unit
            for(int i=0;i<ncontr;i++)mass[i]=(long double)(at->W[i]*amu);

            long double xmntot=mass[0]; //H
            for(int i=1;i<ncontr;i++)
            {
            xmntot=xmntot+mass[i]*(abu[i]); //mass is the atomic mass
            }
//            for(int i=0;i<ncontr;i++)
//            fprintf(stderr,"mass=%Le abu=%Le\n",mass[i],abu[i]);
 

////////////Compute N_tot for a given rho///////////////////////////////
           long double rh=rho_in;
/////////Compute N_tot for a given rho//////////////////////////////////
//            fprintf(stderr,"xmntot=%Le\n",xmntot);
            n_i[0]=rh/xmntot; ////Get n_H
//            fprintf(stderr,"nH=%Le\n",n_i[0]);
            xnhtot=n_i[0];

            xntot=0.L;
            for (int i=0;i<ncontr;i++)
            {
            n_i[i]=(abu[i])*xnhtot; //// n_H times abundance 
            xntot=xntot+n_i[i];
//            fprintf(stderr,"i=%d n_i=%Le\n",i,n_i[i]);
            }
#ifdef print_ions
            fprintf(stderr,"rho=%Le eps=%Le\n",rho_in,eps_in);
#endif
             
            uu1=new long double[ncontr];
            uu2=new long double[ncontr];
            uu3=new long double[ncontr];
///////////////BISECTION METHOD FOR FINDING T////////////////////////////
//////////////////////////TMIN///////////////////////////////////////////
            T=Tmin;
            T=T>TK?T:TK;
            T_a=T;
            
            
            E_min=compute_energy(T,xnhtot,xntot);
            E_a=eps_in-E_min/rh;

#ifdef print_ions
            fprintf(stderr,"T_min=%Le\n",Tmin);
            fprintf(stderr,"E_a=%Le E_min/rh=%Le\n",E_a,E_min/rh);
//            fprintf(stderr,"sumi=%Le sumexi=%Le ski=%Le N_el=%Le ntot=%Le pe=%Le\n",
//              sumi,sumexi,ski,N_el,ntot,pe);
#endif
////////////////////////////////////////////////////////
//////////////////TMAX//////////////////////////////////
            T=Tmax;
            T=T>TK?T:TK;
            T_b=T;
//            fprintf(stderr,"Temp=%Le\n",T);
            E_max=compute_energy(T,xnhtot,xntot);
            E_b=eps_in-E_max/rh;

#ifdef print_ions
            fprintf(stderr,"T_max=%Le\n",Tmax);
            fprintf(stderr,"E_b=%Le E_max/rh=%Le\n",E_b,E_max/rh);
//            fprintf(stderr,"sumi=%Le sumexi=%Le ski=%Le N_el=%Le ntot=%Le\n",sumi,sumexi,ski,N_el,ntot);
#endif
            if(E_a > 0.L && E_b > 0.L) 
            {
            fprintf(stderr,"T or eps out of bounds\n");
//            exit(0);
//            T=Tmax;
            T=Temp_in;
            iter_en=0;
/////////////////////////////////////////////////////////////////////////////
            for(int i=0;i<ncontr;i++)
             {
             int ii=i+1;
             float64_t Tdd=(double)T;
             at->partfun(ii,Tdd);
             uu1[i]=(long double)at->u1;
             uu2[i]=(long double)at->u2;
             uu3[i]=(long double)at->u3;
             }
            pe=inner_iteration(T,xnhtot,xntot);
            N_el=pe/(kbl*T);
            P_g=(N_el+ntot)*kbl*T; //update when computed
/////////////////////////////////////////////////////////////////////////////

            }
            if(E_a < 0.L && E_b < 0.L) 
            {
            fprintf(stderr,"T or eps out of bounds\n");
//            exit(0);
//            T=Tmin;
            T=Temp_in;
            iter_en=0;
/////////////////////////////////////////////////////////////////////////////
            for(int i=0;i<ncontr;i++)
             {
             int ii=i+1;
             float64_t Tdd=(double)T;
             at->partfun(ii,Tdd);
             uu1[i]=(long double)at->u1;
             uu2[i]=(long double)at->u2;
             uu3[i]=(long double)at->u3;
             }
            pe=inner_iteration(T,xnhtot,xntot);
            N_el=pe/(kbl*T);
            P_g=(N_el+ntot)*kbl*T; //update when computed
/////////////////////////////////////////////////////////////////////////////
            }
         

//            fprintf(stderr,"Tmin=%Le Tmax=%Le\n",Tmin,Tmax);
//            fprintf(stderr,"eps_min=%Le eps_max=%Le\n",E_min/rh,E_max/rh);
             E_current=E_min; ////Initialize
//            fprintf(stderr,"E_a0=%Le E_b0=%Le\n",E_a,E_b);
///////////////////////////////////////////////////////////////////
////////////////////////////////////T-search////////////////////////
          if(iter_en==1) 
          {
          T=0.5L*(Tmin+Tmax);
          T_search=T;//initialize

          rel_err=1.L+19;
//        fprintf(stderr,"crec=%Le\n",crec);
          iterTn=0;
////////INNER ITERATION STARTS///////////////////////////////////////
           do 
            {
            iterTn=iterTn+1;
            T=T>TK?T:TK;
            E_old=E_current;
            E_current=compute_energy(T,xnhtot,xntot);
            E_search=eps_in-(E_current)/rh; ////Difference in energy

            T_prev=T;
               if(E_search>0.L) 
                 {
                 if(E_a<0.L) 
                  {
                   ////replace(b,f(b)//////
                   T_b=T;
                   E_b=E_search;
                   T_search=0.5L*(T+T_a);
                   }
                 if(E_b<0.L) 
                   {
                   ////replace(a,f(a)/////
                   T_a=T;
                   E_a=E_search;
                   T_search=0.5L*(T+T_b);
                   }
                 }
                if(E_search<0.L) 
                 {
                 if(E_a>0.L) 
                  {
                  ////replace(b,f(b)//////
                  T_b=T;
                  E_b=E_search;
                  T_search=0.5L*(T+T_a);
                  }
                 if(E_b>0.L) 
                  {
                   ////replace(a,f(a)//////
                  T_a=T;
                  E_a=E_search;
                  T_search=0.5L*(T+T_b);
                  }
                 }
//            fprintf(stderr,"T=%Le T_search=%Le\n",T,T_search);
//            rel_err=fabs(eps_in-(E_current/rh))/eps_in;
//            rel_err=fabs(T_search-T_prev)/T_search;
            rel_err=fabs(T_search-T_prev);
            T=T_search; ////Update T
#ifdef print_ions
            fprintf(stderr,"T_a=%Le T_b=%Le T_search=%Le diff=%Le\n",T_a,T_b,T_search,T_search-T_prev);
            fprintf(stderr,"E_a=%Le E_b=%Le E_search=%Le\n",E_a,E_b,E_search);
            fprintf(stderr,"E_old=%Le E_current=%Le\n",E_old,E_current);
            fprintf(stderr,"eps_in=%Le e=%Le\n",eps_in,E_current/rh); 
            fprintf(stderr,"iterTn=%d rel_err=%Le\n",iterTn,rel_err);

            long double T_diff=(T_search-Temp_in)/(T_search+Temp_in);
            long double Nel_diff=(N_el-Nel_in)/(N_el+Nel_in);
            fprintf(stderr,"T_diff=%Le Nel_diff=%Le\n",T_diff,Nel_diff);  
#endif
#ifdef print_ions
         fprintf(stderr,"Hp_1=%Le Hn_1=%Le Hm_1=%Le H2p_1=%Le H2n_1=%Le\n",Hp_1,Hn_1,Hm_1,H2p_1,H2n_1);
         fprintf(stderr," \n");
#endif
            }while(rel_err > crec && iterTn < itermax);

            if(iterTn==itermax) 
            {
            fprintf(stderr,"Bisection method did not find a root: lte_switch=%d\n",lte_switch);
            fprintf(stderr,"T=%Le Hp_1=%Le Hn_1=%Le\n",T,Hp_1,Hn_1);
            fprintf(stderr,"H2p_1=%Le H2n_1=%Le\n",H2p_1,H2n_1);
             exit(0);
            }
//////////////////////////////////////////////////////////////
//            T=0.5L*(T_a+T_b);
//            T=0.5L*(T+T_prev);
//////////////////////////////////////////////////////////////////////
            pe=inner_iteration(T,xnhtot,xntot);
            N_el=pe/(kbl*T);
/////////////////////////////////////////////////////////////////////
            Temp=T;
            P_g=(N_el+ntot)*kbl*T; //update when computed
          }//iter_en condition
//////////////////output/////////////////////////////////////
            T_out=T; 
            Pg_out=P_g;
            Pe_out=pe;
            Nel_out=N_el;
            nHt=Hp_1+Hn_1; //atomic hydrogen
//////////////////output/////////////////////////////////////
#ifdef print_ions
         fprintf(stderr,"Hp_1=%Le Hn_1=%Le Hm_1=%Le H2p_1=%Le H2n_1=%Le\n",Hp_1,Hn_1,Hm_1,H2p_1,H2n_1);
         fprintf(stderr,"charge cons: Hp_1=%Le Hm_1=%Le H2p_1=%Le Nel=%Le g1*nh=%Le ch-n=%Le ch-p=%Le\n",
                 Hp_1,Hm_1,H2p_1,Nel_out,g1*phtot/(kbl*T),Hp_1+H2p_1+g1*phtot/(kbl*T),Hm_1+Nel_out);
         fprintf(stderr," \n");
        
         fprintf(stderr,"T=%Le Pg=%Le Pe=%Le Nel=%Le nHt=%Le\n",T_out,Pg_out,Pe_out,Nel_out,nHt);
         fprintf(stderr,"P_g=%Le Rho=%Le Eps=%Le\n",P_g,rh,eps_in);
         fprintf(stderr,"End iteration\n");
         fprintf(stderr," \n");
#endif
         delete [] uu1;
         delete [] uu2;
         delete [] uu3;
         delete [] chi1;
         delete [] chi2;
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
       
}

