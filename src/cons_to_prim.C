#include <mpi.h>
#include <math.h>
#include <limits>
#include <stdlib.h>
#include "grid.H"
#include "physics.H"
#include "precision.h"
#include "eos.H"
#include "run.H"

double eos_gamma = 1.65;
double eos_mu = 0.62;
double eos_mu_n = 1.3;
double mh = 1.67262158e-24;
double kb = 1.380622e-16; 

inline int invalid_eos(double a, double b){
  int status = 0;

  if( (a != a) or (a <= 0.) or (a == std::numeric_limits<double>::infinity()) )
    status = 1;

  if( (b != b) or (b <= 0.) or (b == std::numeric_limits<double>::infinity()) )
    status = 1;

  return status;
}

void ConsToPrim(GridData& Grid, const PhysicsData& Physics, const RunData& Run) {

  //double time,s_time;
  //static double t_time = 0.0 ,c_time = 0.0 , r_time = 0.0;
  //static int call_count = 0;
  //s_time = MPI_Wtime();
  
  int i,j,k,node,i_d,i_e;
  double ep,lr,dn,vv,xx;

  const double eps_max=exp(xeps[N_eps-1])-eps_off;

  int sz = Grid.lsize[0]+2*Grid.ghosts[0];
  const int v_nvar = Grid.v_nvar;
  /* ideal gas law for out of bound values */
  const double Rgas=8.314e7;

  const double c_temp= (eos_gamma-1.)*eos_mu/Rgas;
  const double c_pres= (eos_gamma-1.);  

  cState* U = Grid.U;

  double var[v_nvar][sz], coeff[4][sz], ptab[4][sz],
  netab[4][sz], rhoitab[4][sz],ambtab[4][sz],ttab[4][sz], lge[sz],eps[sz], 
  lgr[sz], pres[sz], temp[sz],nel[sz],amb[sz],rhoi[sz];

  for(k=Grid.lbeg[2]-Grid.ghosts[2];k<=Grid.lend[2]+Grid.ghosts[2];k++){
    for(j=Grid.lbeg[1]-Grid.ghosts[1];j<=Grid.lend[1]+Grid.ghosts[1];j++){
      //time = MPI_Wtime();
      for(i=0;i<sz;i++){
        node = Grid.node(i,j,k);
        var[0][i] = U[node].d;
        var[1][i] = U[node].M.x;
        var[2][i] = U[node].M.y;
        var[3][i] = U[node].M.z;
        var[4][i] = U[node].e;
        var[5][i] = U[node].B.x;
        var[6][i] = U[node].B.y;
        var[7][i] = U[node].B.z;
      }
      //r_time += MPI_Wtime()-time;

      //time = MPI_Wtime();
      for(i=0;i<sz;i++){
        dn        = 1.0/var[0][i];
        var[1][i] = var[1][i]*dn;
        var[2][i] = var[2][i]*dn;
        var[3][i] = var[3][i]*dn;

        vv  = (var[1][i]*var[1][i]+
        var[2][i]*var[2][i]+
        var[3][i]*var[3][i]);

        var[4][i] = var[4][i] - 0.5*vv*var[0][i];
  
        eps[i] = var[4][i]*dn +eps_off;

        lgr[i] = log(var[0][i]);
        lge[i] = log(eps[i]);
      }
      for(i=0;i<sz;i++){

        ep = lge[i];
        lr = lgr[i];

	i_d = (int) ( (lr-xlr[0])*inv_del_lr );
        i_d = max(0,i_d);
        i_d = min(N_lr-2,i_d);
        
        i_e = (int) ((ep-xeps[0])*inv_del_eps );
        i_e = max(0,i_e);
        i_e = min(N_eps-2,i_e);
  
        coeff[0][i] = (ep-xeps[i_e])   * (lr-xlr[i_d]);
        coeff[1][i] = (ep-xeps[i_e])   * (xlr[i_d+1]-lr);
        coeff[2][i] = (xeps[i_e+1]-ep) * (lr-xlr[i_d]);
        coeff[3][i] = (xeps[i_e+1]-ep) * (xlr[i_d+1]-lr);
          
        ptab[0][i]  = (double) p_eostab[i_e+1][i_d+1];
        ptab[1][i]  = (double) p_eostab[i_e+1][i_d];
        ptab[2][i]  = (double) p_eostab[i_e][i_d+1];
        ptab[3][i]  = (double) p_eostab[i_e][i_d];

        ttab[0][i]  = (double) T_eostab[i_e+1][i_d+1];
        ttab[1][i]  = (double) T_eostab[i_e+1][i_d];
        ttab[2][i]  = (double) T_eostab[i_e][i_d+1];
        ttab[3][i]  = (double) T_eostab[i_e][i_d];  

        netab[0][i]  = (double) ne_eostab[i_e+1][i_d+1];
        netab[1][i]  = (double) ne_eostab[i_e+1][i_d];  
        netab[2][i]  = (double) ne_eostab[i_e][i_d+1];  
        netab[3][i]  = (double) ne_eostab[i_e][i_d];    

        rhoitab[0][i]  = (double) rhoi_eostab[i_e+1][i_d+1];
        rhoitab[1][i]  = (double) rhoi_eostab[i_e+1][i_d];  
        rhoitab[2][i]  = (double) rhoi_eostab[i_e][i_d+1];  
        rhoitab[3][i]  = (double) rhoi_eostab[i_e][i_d];    

        ambtab[0][i]  = (double) amb_eostab[i_e+1][i_d+1];
        ambtab[1][i]  = (double) amb_eostab[i_e+1][i_d];  
        ambtab[2][i]  = (double) amb_eostab[i_e][i_d+1];  
        ambtab[3][i]  = (double) amb_eostab[i_e][i_d];    

      }
  
      for(i=0;i<sz;i++){
        double pp  = coeff[0][i] * ptab[0][i];
        pp += coeff[1][i] * ptab[1][i];
        pp += coeff[2][i] * ptab[2][i];
        pp += coeff[3][i] * ptab[3][i];
        
        double tt  = coeff[0][i] * ttab[0][i];
        tt += coeff[1][i] * ttab[1][i];
        tt += coeff[2][i] * ttab[2][i];
        tt += coeff[3][i] * ttab[3][i];

        double nn  = coeff[0][i] * netab[0][i];
        nn += coeff[1][i] * netab[1][i];
        nn += coeff[2][i] * netab[2][i];
        nn += coeff[3][i] * netab[3][i];
       
        double ii = coeff[0][i] * rhoitab[0][i]; 
        ii += coeff[1][i] * rhoitab[1][i]; 
        ii += coeff[2][i] * rhoitab[2][i]; 
        ii += coeff[3][i] * rhoitab[3][i];

        double aa  = coeff[0][i] * ambtab[0][i];
        aa += coeff[1][i] * ambtab[1][i];
        aa += coeff[2][i] * ambtab[2][i];
        aa += coeff[3][i] * ambtab[3][i];

        pres[i] = exp(pp);
        temp[i] = exp(tt);
        nel[i] = exp(nn);
        rhoi[i] = exp(ii);
        amb[i] = exp(aa);
      }
      
      for(i=0;i<sz;i++){
    
      // Out of upper table bounds
      if( (lgr[i] < xlr[0]) ) {
        pres[i]=c_pres*var[4][i];
        temp[i]=c_temp*eps[i];
        // Use PV = nrt for ne
        nel[i]= pres[i]/(temp[i]*kb)-var[0][i]/(eos_mu_n*mh);
      }
  
      if (eps[i] > 0.8*eps_max) {
        xx = max(0.0,5.0*(eps_max-eps[i])/eps_max);
        pres[i]=pres[i]*xx+(1.0-xx)*c_pres*var[4][i];
        temp[i]=temp[i]*xx+(1.0-xx)*c_temp*eps[i];
        // Use PV = nrt for ne
        nel[i]= pres[i]/(temp[i]*kb)-var[0][i]/(eos_mu_n*mh);
      }
  
      if(invalid_eos(pres[i],temp[i])){
        pres[i]=c_pres*var[4][i];
        temp[i]=c_temp*eps[i];
      }

      // if nel is inf use pv=nrt
      if(nel[i] == std::numeric_limits<double>::infinity()) nel[i]= pres[i]/(temp[i]*kb)-var[0][i]/(eos_mu_n*mh);

      // If electron number goes too small, or NaNs, or is still inf set to a minimum
      if ((nel[i] < 1.0) || (nel[i]!=nel[i]) || (nel[i] == std::numeric_limits<double>::infinity())) nel[i] = 1.0;

      }

      //time = MPI_Wtime();
      for(i=Grid.lbeg[0]-Grid.ghosts[0];i<=Grid.lend[0]+Grid.ghosts[0];i++){
        node = Grid.node(i,j,k);
        U[node].M.x=var[1][i];
        U[node].M.y=var[2][i];
        U[node].M.z=var[3][i];
        U[node].e  =var[4][i];
        
        Grid.pres[node] = pres[i];
        Grid.temp[node] = temp[i];
        Grid.ne[node] = nel[i];
        Grid.rhoi[node] = rhoi[i];
        Grid.amb[node] = amb[i];

      }
      //r_time += MPI_Wtime()-time;

    }
  }
}
