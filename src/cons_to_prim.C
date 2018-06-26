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

#pragma acc routine seq
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

  //double var[v_nvar][sz], 
  //double coeff[4][sz], 
  //double ptab[4][sz],
  //double netab[4][sz], 
  //double rhoitab[4][sz],
  //double ambtab[4][sz],
  //double ttab[4][sz], 
  //double lge[sz],eps[sz], 
  //lgr[sz], pres[sz], temp[sz],nel[sz],amb[sz],rhoi[sz];
  double *var = (double*) malloc(v_nvar*sz*sizeof(double));
  double *coeff = (double*) malloc(4*sz*sizeof(double));
  double *ptab = (double*) malloc(4*sz*sizeof(double));
  double *netab = (double*) malloc(4*sz*sizeof(double));
  double *rhoitab = (double*) malloc(4*sz*sizeof(double));
  double *ambtab = (double*) malloc(4*sz*sizeof(double));
  double *ttab = (double*) malloc(4*sz*sizeof(double));
  double *lge = (double*) malloc(sz*sizeof(double));
  double *eps = (double*) malloc(sz*sizeof(double));
  double *lgr = (double*) malloc(sz*sizeof(double));
  double *pres = (double*) malloc(sz*sizeof(double));
  double *temp = (double*) malloc(sz*sizeof(double));
  double *nel = (double*) malloc(sz*sizeof(double));
  double *amb = (double*) malloc(sz*sizeof(double));
  double *rhoi = (double*) malloc(sz*sizeof(double));

  int k1 = Grid.lbeg[2]-Grid.ghosts[2];
  int k2 = Grid.lend[2]+Grid.ghosts[2];
  int j1 = Grid.lbeg[1]-Grid.ghosts[1];
  int j2 = Grid.lend[1]+Grid.ghosts[1];
  int i1 = Grid.lbeg[0]-Grid.ghosts[0];
  int i2 = Grid.lend[0]+Grid.ghosts[0];
  int str1 = Grid.stride[0];
  int str2 = Grid.stride[1];
  int str3 = Grid.stride[2];

  eos_real *p_eostab_flat = p_eostab[0];
  eos_real *T_eostab_flat = T_eostab[0];
  eos_real *s_eostab_flat = s_eostab[0];
  eos_real *ne_eostab_flat = ne_eostab[0];
  eos_real *rhoi_eostab_flat = rhoi_eostab[0];
  eos_real *amb_eostab_flat = amb_eostab[0];

  //for(k=Grid.lbeg[2]-Grid.ghosts[2];k<=Grid.lend[2]+Grid.ghosts[2];k++){
  //  for(j=Grid.lbeg[1]-Grid.ghosts[1];j<=Grid.lend[1]+Grid.ghosts[1];j++){


#pragma acc enter data copyin(Grid)
#pragma acc enter data copyin(U[0:Grid.bufsize],Grid.pres[0:Grid.bufsize], \
  Grid.temp[0:Grid.bufsize],Grid.ne[0:Grid.bufsize],Grid.rhoi[0:Grid.bufsize], \
  Grid.amb[0:Grid.bufsize],xlr[0:N_lr],xeps[0:N_eps],p_eostab_flat[0:N_eps*N_lr], \
 T_eostab_flat[0:N_eps*N_lr],ne_eostab_flat[0:N_eps*N_lr],rhoi_eostab_flat[0:N_eps*N_lr], \
 amb_eostab_flat[0:N_eps*N_lr])

#pragma acc parallel default(present)
{
#pragma acc loop collapse(2) gang independent private(var[0:v_nvar*sz], \
  eps[0:sz],lgr[0:sz],lge[0:sz],coeff[0:4*sz],ptab[0:4*sz],ttab[0:4*sz], \
  netab[0:4*sz],rhoitab[0:4*sz],ambtab[0:4*sz],pres[0:sz],temp[0:sz], \
  nel[0:sz],rhoi[0:sz],amb[0:sz])
  for(k=k1; k<=k2; k++){
    for(j=j1; j<=j2; j++){
      //time = MPI_Wtime();
#pragma acc loop independent vector private(node)
      for(i=0;i<sz;i++){
        //node = Grid.node(i,j,k);
        node = i*str1+j*str2+k*str3;
        var[i*v_nvar+0] = U[node].d;
        var[i*v_nvar+1] = U[node].M.x;
        var[i*v_nvar+2] = U[node].M.y;
        var[i*v_nvar+3] = U[node].M.z;
        var[i*v_nvar+4] = U[node].e;
        var[i*v_nvar+5] = U[node].B.x;
        var[i*v_nvar+6] = U[node].B.y;
        var[i*v_nvar+7] = U[node].B.z;
      }
      //r_time += MPI_Wtime()-time;

      //time = MPI_Wtime();
#pragma acc loop independent vector private(dn,vv)
      for(i=0;i<sz;i++){
        dn        = 1.0/var[i*v_nvar+0];
        var[i*v_nvar+1] = var[i*v_nvar+1]*dn;
        var[i*v_nvar+2] = var[i*v_nvar+2]*dn;
        var[i*v_nvar+3] = var[i*v_nvar+3]*dn;

        vv  = (var[i*v_nvar+1]*var[i*v_nvar+1]+
        var[i*v_nvar+2]*var[i*v_nvar+2]+
        var[i*v_nvar+3]*var[i*v_nvar+3]);

        var[i*v_nvar+4] = var[i*v_nvar+4] - 0.5*vv*var[i*v_nvar+0];
  
        eps[i] = var[i*v_nvar+4]*dn +eps_off;

        lgr[i] = log(var[i*v_nvar+0]);
        lge[i] = log(eps[i]);
      }
#pragma acc loop independent vector private(ep,lr,i_d,i_e)
      for(i=0;i<sz;i++){

        ep = lge[i];
        lr = lgr[i];

	      i_d = (int) ( (lr-xlr[0])*inv_del_lr );
        i_d = max(0,i_d);
        i_d = min(N_lr-2,i_d);
        
        i_e = (int) ((ep-xeps[0])*inv_del_eps );
        i_e = max(0,i_e);
        i_e = min(N_eps-2,i_e);
  
        coeff[i*4+0] = (ep-xeps[i_e])   * (lr-xlr[i_d]);
        coeff[i*4+1] = (ep-xeps[i_e])   * (xlr[i_d+1]-lr);
        coeff[i*4+2] = (xeps[i_e+1]-ep) * (lr-xlr[i_d]);
        coeff[i*4+3] = (xeps[i_e+1]-ep) * (xlr[i_d+1]-lr);
          
        //ptab[i*4+0]  = (double) p_eostab[i_e+1][i_d+1];
        //ptab[i*4+1]  = (double) p_eostab[i_e+1][i_d];
        //ptab[i*4+2]  = (double) p_eostab[i_e][i_d+1];
        //ptab[i*4+3]  = (double) p_eostab[i_e][i_d];
        ptab[i*4+0]  = (double) p_eostab_flat[(i_e+1)*N_lr+(i_d+1)];
        ptab[i*4+1]  = (double) p_eostab_flat[(i_e+1)*N_lr+i_d];
        ptab[i*4+2]  = (double) p_eostab_flat[i_e*N_lr+(i_d+1)];
        ptab[i*4+3]  = (double) p_eostab_flat[i_e*N_lr+i_d];

        //ttab[i*4+0]  = (double) T_eostab[i_e+1][i_d+1];
        //ttab[i*4+1]  = (double) T_eostab[i_e+1][i_d];
        //ttab[i*4+2]  = (double) T_eostab[i_e][i_d+1];
        //ttab[i*4+3]  = (double) T_eostab[i_e][i_d];  
        ttab[i*4+0]  = (double) T_eostab_flat[(i_e+1)*N_lr+(i_d+1)];
        ttab[i*4+1]  = (double) T_eostab_flat[(i_e+1)*N_lr+i_d];
        ttab[i*4+2]  = (double) T_eostab_flat[i_e*N_lr+(i_d+1)];
        ttab[i*4+3]  = (double) T_eostab_flat[i_e*N_lr+i_d];

        //netab[i*4+0]  = (double) ne_eostab[i_e+1][i_d+1];
        //netab[i*4+1]  = (double) ne_eostab[i_e+1][i_d];  
        //netab[i*4+2]  = (double) ne_eostab[i_e][i_d+1];  
        //netab[i*4+3]  = (double) ne_eostab[i_e][i_d]; 
        netab[i*4+0]  = (double) ne_eostab_flat[(i_e+1)*N_lr+(i_d+1)];
        netab[i*4+1]  = (double) ne_eostab_flat[(i_e+1)*N_lr+i_d];
        netab[i*4+2]  = (double) ne_eostab_flat[i_e*N_lr+(i_d+1)];
        netab[i*4+3]  = (double) ne_eostab_flat[i_e*N_lr+i_d];   

        //rhoitab[i*4+0]  = (double) rhoi_eostab[i_e+1][i_d+1];
        //rhoitab[i*4+1]  = (double) rhoi_eostab[i_e+1][i_d];  
        //rhoitab[i*4+2]  = (double) rhoi_eostab[i_e][i_d+1];  
        //rhoitab[i*4+3]  = (double) rhoi_eostab[i_e][i_d];   
        rhoitab[i*4+0]  = (double) rhoi_eostab_flat[(i_e+1)*N_lr+(i_d+1)];
        rhoitab[i*4+1]  = (double) rhoi_eostab_flat[(i_e+1)*N_lr+i_d];
        rhoitab[i*4+2]  = (double) rhoi_eostab_flat[i_e*N_lr+(i_d+1)];
        rhoitab[i*4+3]  = (double) rhoi_eostab_flat[i_e*N_lr+i_d];   

        //ambtab[i*4+0]  = (double) amb_eostab[i_e+1][i_d+1];
        //ambtab[i*4+1]  = (double) amb_eostab[i_e+1][i_d];  
        //ambtab[i*4+2]  = (double) amb_eostab[i_e][i_d+1];  
        //ambtab[i*4+3]  = (double) amb_eostab[i_e][i_d];  
        ambtab[i*4+0]  = (double) amb_eostab_flat[(i_e+1)*N_lr+(i_d+1)];
        ambtab[i*4+1]  = (double) amb_eostab_flat[(i_e+1)*N_lr+i_d];
        ambtab[i*4+2]  = (double) amb_eostab_flat[i_e*N_lr+(i_d+1)];
        ambtab[i*4+3]  = (double) amb_eostab_flat[i_e*N_lr+i_d];  
      }
  
#pragma acc loop independent vector
      for(i=0;i<sz;i++){
        double pp  = coeff[i*4+0] * ptab[i*4+0];
        pp += coeff[i*4+1] * ptab[i*4+1];
        pp += coeff[i*4+2] * ptab[i*4+2];
        pp += coeff[i*4+3] * ptab[i*4+3];
        
        double tt  = coeff[i*4+0] * ttab[i*4+0];
        tt += coeff[i*4+1] * ttab[i*4+1];
        tt += coeff[i*4+2] * ttab[i*4+2];
        tt += coeff[i*4+3] * ttab[i*4+3];

        double nn  = coeff[i*4+0] * netab[i*4+0];
        nn += coeff[i*4+1] * netab[i*4+1];
        nn += coeff[i*4+2] * netab[i*4+2];
        nn += coeff[i*4+3] * netab[i*4+3];
       
        double ii = coeff[i*4+0] * rhoitab[i*4+0]; 
        ii += coeff[i*4+1] * rhoitab[i*4+1]; 
        ii += coeff[i*4+2] * rhoitab[i*4+2]; 
        ii += coeff[i*4+3] * rhoitab[i*4+3];

        double aa  = coeff[i*4+0] * ambtab[i*4+0];
        aa += coeff[i*4+1] * ambtab[i*4+1];
        aa += coeff[i*4+2] * ambtab[i*4+2];
        aa += coeff[i*4+3] * ambtab[i*4+3];

        pres[i] = exp(pp);
        temp[i] = exp(tt);
        nel[i] = exp(nn);
        rhoi[i] = exp(ii);
        amb[i] = exp(aa);
      }
      
#pragma acc loop independent private(xx)
      for(i=0;i<sz;i++){
    
        // Out of upper table bounds
        if( (lgr[i] < xlr[0]) ) {
          pres[i]=c_pres*var[i*v_nvar+4];
          temp[i]=c_temp*eps[i];
          // Use PV = nrt for ne
          nel[i]= pres[i]/(temp[i]*kb)-var[i*v_nvar+4]/(eos_mu_n*mh);
        }
    
        if (eps[i] > 0.8*eps_max) {
          xx = max(0.0,5.0*(eps_max-eps[i])/eps_max);
          pres[i]=pres[i]*xx+(1.0-xx)*c_pres*var[i*v_nvar+4];
          temp[i]=temp[i]*xx+(1.0-xx)*c_temp*eps[i];
          // Use PV = nrt for ne
          nel[i]= pres[i]/(temp[i]*kb)-var[i*v_nvar+4]/(eos_mu_n*mh);
        }
    
        if(invalid_eos(pres[i],temp[i])){
          pres[i]=c_pres*var[i*v_nvar+4];
          temp[i]=c_temp*eps[i];
        }

        // if nel is inf use pv=nrt
        if(nel[i] == std::numeric_limits<double>::infinity()) nel[i]= pres[i]/(temp[i]*kb)-var[i*v_nvar+4]/(eos_mu_n*mh);

        // If electron number goes too small, or NaNs, or is still inf set to a minimum
        if ((nel[i] < 1.0) || (nel[i]!=nel[i]) || (nel[i] == std::numeric_limits<double>::infinity())) nel[i] = 1.0;

      }

      //time = MPI_Wtime();
      //for(i=Grid.lbeg[0]-Grid.ghosts[0];i<=Grid.lend[0]+Grid.ghosts[0];i++){
#pragma acc loop independent vector private(node)
      for(i=i1; i<=i2; i++) {
        //node = Grid.node(i,j,k);
        node = i*str1+j*str2+k*str3;
        U[node].M.x=var[i*v_nvar+1];
        U[node].M.y=var[i*v_nvar+2];
        U[node].M.z=var[i*v_nvar+3];
        U[node].e  =var[i*v_nvar+4];
        
        Grid.pres[node] = pres[i];
        Grid.temp[node] = temp[i];
        Grid.ne[node] = nel[i];
        Grid.rhoi[node] = rhoi[i];
        Grid.amb[node] = amb[i];

      }
      //r_time += MPI_Wtime()-time;

    }
  }

} // end parallel

#pragma acc exit data copyout(U[0:Grid.bufsize],Grid.pres[0:Grid.bufsize], \
  Grid.temp[0:Grid.bufsize],Grid.ne[0:Grid.bufsize],Grid.rhoi[0:Grid.bufsize], \
  Grid.amb[0:Grid.bufsize])
#pragma acc exit data delete(Grid)

  free(var);
  free(coeff);
  free(ptab);
  free(netab);
  free(rhoitab);
  free(ambtab);
  free(ttab);
  free(lge);
  free(eps);
  free(lgr);
  free(pres);
  free(temp);
  free(nel);
  free(amb);
  free(rhoi);

}
