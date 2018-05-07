#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include "physics.H"
#include "grid.H"
#include "run.H"
#include "comm_split.H"
#include "src_int_tck.H"
#include "limit_va.H"
#include "exchange.H"
#include <stdio.h>

#define XZ_LOOP(G,i,k) \
  for((k)=(G).lbeg[2];(k)<=(G).lend[2];(k)++) \
  for((i)=(G).lbeg[0];(i)<=(G).lend[0];(i)++)

#define YZ_LOOP(G,j,k) \
  for((k)=(G).lbeg[2];(k)<=(G).lend[2];(k)++) \
  for((j)=(G).lbeg[1];(j)<=(G).lend[1];(j)++)

#define YZ_FULL_LOOP(G,j,k) \
  for((k)=(G).lbeg[2]-(G).ghosts[2];(k)<=(G).lend[2]+(G).ghosts[2];(k)++) \
  for((j)=(G).lbeg[1]-(G).ghosts[1];(j)<=(G).lend[1]+(G).ghosts[1];(j)++)

#define GLOOP(gs,i,j,d1,d2) \
  for((j)=(gs)[(d2)][0];(j)<=(gs)[(d2)][1];(j)++) \
  for((i)=(gs)[(d1)][0];(i)<=(gs)[(d1)][1];(i)++)

inline int imin(int a, int b) { return a < b ? a : b; }
inline int imax(int a, int b) { return a > b ? a : b; }

void get_damping(const RunData&, const GridData&, const PhysicsData&, double*);
/*****************************************************************************/
void Source_Integrate_Tcheck(const RunData& Run, GridData& Grid, 
			     const PhysicsData& Physics, const int stage) {
  
  int i,j,k,node,off0,v;

  const double wdt = Run.dt/double(maxstage+1-stage);

  const int vsize = Grid.vsize;
  const int i_beg = Grid.lbeg[0];
  const int i_end = Grid.lend[0];
  
  const int H_loss = imax(imin(Physics.rt_ext[i_ext_hlines],1),0);
  const int Mg_loss = imax(imin(Physics.rt_ext[i_ext_mglines],1),0);
  const int Ca_loss = imax(imin(Physics.rt_ext[i_ext_calines],1),0);

  int next[3];
  for(v=0;v<3;v++)
    next[v] = Grid.stride[v];
   
  const double eps_min = Physics.tchk[i_tchk_eps_min];
  const double rho_min = Physics.tchk[i_tchk_rho_min];

  const bool ambipolar = (Physics.params[i_param_ambipolar] > 0.0);
  const bool spitzer   = (Physics.params[i_param_spitzer] > 0.0);

  double vmax,vv,ekin,eps_max;
  double ee_low[vsize], ee_up[vsize], vfac[vsize], rho[vsize], Qrad[vsize], efac[vsize], ee[vsize], e_dmp[vsize];

  if(v_lim > 0.0)
    vmax = v_lim;
  else
    vmax = Physics.tchk[i_tchk_vmax];

  if(e_lim > 0.0)
    eps_max = e_lim;
  else
    eps_max = Physics.tchk[i_tchk_eps_max];

  const double g = Physics.g[0];

  double dn,fsr_x,fsr_y,fsr_z,va2,bb,bf,s,x2,x4;
  double rh,vx,vy,vz,bx,by,bz,rd;

  const double invdt  = 1.0/Run.dt; 
  const double invRdt = Physics.rt[i_rt_cfl]/Run.dt;
  
  double boris[4][vsize];

  int sz = Grid.lsize[0]+2*Grid.ghosts[0];
  static double* dmp =  (double*) malloc(sz*sizeof(double));

  if(stage == 1) get_damping(Run,Grid,Physics,dmp);
   
  YZ_LOOP(Grid,j,k){
    off0 = j*next[1]+k*next[2];

    if(Physics.rt_ext[i_ext_cor]==1){
      #pragma ivdep
      for(i=i_beg;i<=i_end;i++){
	node = off0+i;
	Qrad[i] = Grid.Qtot[node]+Grid.Qthin[node];
      }
    } else if(Physics.rt_ext[i_ext_cor]==2){
      #pragma ivdep
      for(i=i_beg;i<=i_end;i++){
	node = off0+i;
	Qrad[i] = Grid.Qtot[node]+Grid.Qthin[node]+Grid.QH[node]*H_loss
	  +Grid.QMg[node]*Mg_loss + Grid.QCa[node]*Ca_loss;
      }
    } else {
      #pragma ivdep
      for(i=i_beg;i<=i_end;i++){
	node = off0+i;
	Qrad[i] = Grid.Qtot[node];
      }
    }
    
    #pragma ivdep
    for(i=i_beg;i<=i_end;i++){
      node = off0+i;
 
      dn      = invRdt*Grid.U[node].e;
      efac[i] = dn/max(dn,fabs(Qrad[i]));
      
      // gravity + damping + RT
      dn       = g*Grid.U[node].d-dmp[i];
      e_dmp[i] = dn*Grid.U[node].M.x; 
      Grid.Res[node].M.x += dn;
      Grid.Res[node].e   += e_dmp[i] + Qrad[i]*efac[i];
    }

    if(Run.need_diagnostics){
      #pragma ivdep
      for(i=i_beg;i<=i_end;i++){
	node = off0+i;
	// diagnostics
	Grid.tvar4[node] = efac[i];
	Grid.tvar5[node] = e_dmp[i];
      }
    }

    // Boris correction
    #pragma ivdep
    for(i=i_beg;i<=i_end;i++){
      node = off0+i;

      rh = Grid.U[node].d;
      rd = Grid.Res[node].d;
      vx = Grid.U[node].M.x;
      vy = Grid.U[node].M.y;
      vz = Grid.U[node].M.z;
      bx = Grid.U[node].B.x;
      by = Grid.U[node].B.y;
      bz = Grid.U[node].B.z;
      
      fsr_x = Grid.Res[node].M.x; 
      fsr_y = Grid.Res[node].M.y;
      fsr_z = Grid.Res[node].M.z;

      bb  = bx*bx+by*by+bz*bz;
      va2 = bb/rh;

      // approximation of 1-1/sqrt(1+(va/c)^4)
      x2 = va2*inv_va2max;
      x4 = x2*x2;
      bf = x4/(1+x2+x4);

      // convert -div(rho*v:v) to -rho*(v*grad)v by adding v*div(rho*v)
      fsr_x -= vx*rd;
      fsr_y -= vy*rd;
      fsr_z -= vz*rd;
      
      dn = (fsr_x*bx+fsr_y*by+fsr_z*bz)/max(1e-100,bb);
      
      boris[0][i] = -bf*(fsr_x-dn*bx);
      boris[1][i] = -bf*(fsr_y-dn*by);
      boris[2][i] = -bf*(fsr_z-dn*bz);
      
      boris[3][i] = boris[0][i]*vx+boris[1][i]*vy+boris[2][i]*vz;
    }
	  
    #pragma ivdep
    for(i=i_beg;i<=i_end;i++){
      node = off0+i;	
      Grid.Res[node].M.x += boris[0][i];
      Grid.Res[node].M.y += boris[1][i];
      Grid.Res[node].M.z += boris[2][i];
      Grid.Res[node].e   += boris[3][i];
    }

    if(Run.need_diagnostics){
      #pragma ivdep
      for(i=i_beg;i<=i_end;i++){
	node = off0+i;	
	Grid.tvar6[node] =  boris[3][i];
      }
    }
    
    // time integration - MHD
    #pragma ivdep
    for(i=i_beg;i<=i_end;i++){
      node = off0+i;      
      Grid.U[node]   = Grid.U0[node] + wdt*Grid.Res[node];
    }
    memset(&Grid.Res[off0+i_beg],0.0,(i_end-i_beg+1)*sizeof(cState));

    // time integration - hyperbolic conduction   
    if(spitzer){
      #pragma ivdep
      for(i=i_beg;i<=i_end;i++){
	node = off0+i;
	Grid.sflx[node] = Grid.sflx0[node]+wdt*Grid.Rflx[node];
      }
      memset(&Grid.Rflx[off0+i_beg],0.0,(i_end-i_beg+1)*sizeof(double));
    }
       
    if(ambipolar){
      #pragma ivdep
      for(i=i_beg;i<=i_end;i++){
	node = off0+i;	
	Grid.v_amb[node] = Grid.v0_amb[node]+wdt*Grid.R_amb[node];
      }
      memset(&Grid.R_amb[off0+i_beg],0.0,(i_end-i_beg+1)*sizeof(Vector));
    }

    // Tcheck - catch extreme values before they cause problems
    #pragma ivdep
    for(i=i_beg;i<=i_end;i++){
      node = off0+i;
      
      dn = max(rho_min,Grid.U[node].d);

      vv = Grid.U[node].M.abs()/dn;
      s  = vmax/max(vmax,vv);

      ekin = 0.5*dn*vv*vv;
      
      ee_low[i] = eps_min*dn+ekin;
      ee_up[i]  = eps_max*dn+ekin*s*s;
      vfac[i]   = s;
      rho[i]    = dn;
    }

    #pragma ivdep
    for(i=i_beg;i<=i_end;i++){ 
      node = off0+i;

      ee[i] = Grid.U[node].e;
      
      Grid.U[node].d  = rho[i];
      Grid.U[node].M *= vfac[i];
      Grid.U[node].e  = min(max(ee[i],ee_low[i]),ee_up[i]);
    }

    if(Run.need_diagnostics){
      #pragma ivdep
      for(i=i_beg;i<=i_end;i++){
	node = off0+i;	   
	Grid.tvar7[node] += (Grid.U[node].e-ee[i])*invdt;
      }
    }

    if(ambipolar){
      #pragma ivdep
      for(i=i_beg;i<=i_end;i++){
	node = off0+i;
      
	vv = Grid.v_amb[node].abs();
	s  = Physics.params[i_param_ambvel_max]/max(Physics.params[i_param_ambvel_max],vv);
	Grid.v_amb[node] *= s;
      }
    }

  }
 
}
/*****************************************************************************/
void SaveCons(GridData& Grid,const PhysicsData& Physics) {

  const int i_beg    = Grid.lbeg[0];
  const int i_end    = Grid.lend[0];

  int off0,j,k;
  int bufsz = i_end-i_beg+1;

  int next[3];
  for(k=0;k<3;k++)
    next[k] = Grid.stride[k];
  
  YZ_LOOP(Grid,j,k){
    off0 = j*next[1]+k*next[2];
    memcpy(&Grid.U0[off0+i_beg],&Grid.U[off0+i_beg],bufsz*sizeof(cState));
    if(Physics.params[i_param_spitzer] > 0.0)
      memcpy(&Grid.sflx0[off0+i_beg],&Grid.sflx[off0+i_beg],bufsz*sizeof(double));
    if(Physics.params[i_param_ambipolar] > 0.0)
    memcpy(&Grid.v0_amb[off0+i_beg],&Grid.v_amb[off0+i_beg],bufsz*sizeof(Vector));
  }

} 
/*****************************************************************************/
void TCheck(const RunData&  Run, GridData& Grid, const PhysicsData& Physics) {

  static int ini_flag = 1;

  // called after Integrate(); U contains conservative Vars.

  register int i,j,k,v,off0,node;
 
  const int vsize = Grid.vsize;

  const double eps_min = Physics.tchk[i_tchk_eps_min];
  const double rho_min = Physics.tchk[i_tchk_rho_min];
  
  double dn,vmax,vv,ekin,s,eps_max;
  double ee_low[vsize], ee_up[vsize], vfac[vsize], rho[vsize], ee[vsize];

  const int i_beg = Grid.lbeg[0];
  const int i_end = Grid.lend[0];

  double invdt;
  if(Run.dt > 0.0)
    invdt = 1.0/Run.dt;
  else
    invdt = 0.0;

  int next[3];
  for(v=0;v<3;v++)
    next[v] = Grid.stride[v];

  if(v_lim > 0.0)
    vmax = v_lim;
  else
    vmax = Physics.tchk[i_tchk_vmax];

  if(e_lim > 0.0)
    eps_max = e_lim;
  else
    eps_max = Physics.tchk[i_tchk_eps_max];
  
  if(ini_flag){
    if (Run.rank == 0) {
      cout << "tcheck: rho_min   = "  << rho_min   << endl;
      cout << "tcheck: eps_min   = "  << eps_min   << endl;
      cout << "tcheck: eps_max   = "  << eps_max   << endl;    
      cout << "tcheck: vmax      = "  << vmax      << endl;
    }
  
    ini_flag = 0;
  }

  YZ_LOOP(Grid,j,k){
    off0 = j*next[1]+k*next[2];

    #pragma ivdep
    for(i=i_beg;i<=i_end;i++){
      node = off0+i;
      
      dn = max(rho_min,Grid.U[node].d);

      vv = Grid.U[node].M.abs()/dn;
      s  = vmax/max(vmax,vv);

      ekin = 0.5*dn*vv*vv;
      
      ee_low[i] = eps_min*dn+ekin;
      ee_up[i]  = eps_max*dn+ekin*s*s;
      vfac[i]   = s;
      rho[i]    = dn;
    }

    #pragma ivdep
    for(i=i_beg;i<=i_end;i++){
      node = off0+i;
 
      ee[i] = Grid.U[node].e;
      
      Grid.U[node].d  = rho[i];
      Grid.U[node].M *= vfac[i];
      Grid.U[node].e  = min(max(ee[i],ee_low[i]),ee_up[i]);
    }
    
    if(Run.need_diagnostics){
      #pragma ivdep
      for(i=i_beg;i<=i_end;i++){
	node = off0+i;	   
	Grid.tvar7[node] += (Grid.U[node].e-ee[i])*invdt;
      }
    }
      
  }

}
/*****************************************************************************/
void get_damping(const RunData&  Run, const GridData& Grid,
		 const PhysicsData& Physics, double* my_mean) {

  static int ini_flag = 1;

  int i,j,k,node;

  double hh,hmax,hphot,fmax,vmax;
  double sbuf[4],rbuf[4];

  int hsize = Grid.gsize[1]*Grid.gsize[2];
  int vsize = Grid.lsize[0]+2*Grid.ghosts[0];

  static double *lbuf    = (double*) malloc(vsize*sizeof(double));
  static double *rh_mean = (double*) malloc(vsize*sizeof(double));
  static double rho_max;

  const double tau_ref = 1.0/Physics.dmp[i_dmp_tau_ref];
  const double vel_ref = Physics.dmp[i_dmp_vel_ref];
  const double tau_min = 1.0/Physics.dmp[i_dmp_tau_min];

  //-------------------------------------------------------------------

  if (ini_flag == 1){
    hh = 0.0;
    for(i=0;i<vsize;i++) lbuf[i] = 0.0;
    LOCAL_LOOP(Grid,i,j,k){
      node     = Grid.node(i,j,k);
      hh       = max(hh,Grid.U[node].d);
      lbuf[i] += Grid.U[node].d/hsize;  
    }
    MPI_Allreduce(&hh,&rho_max,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    MPI_Allreduce(lbuf,rh_mean,vsize,MPI_DOUBLE,MPI_SUM,YZ_COMM);

    if ( Run.rank == 0 ) {
      cout << "damp: tau_ref = " << tau_ref         << endl;
      cout << "damp: vel_ref = " << vel_ref         << endl;
      cout << "damp: tau_min = " << tau_min         << endl;
      cout << "damp: mflx    = " << vel_ref*rho_max << endl;
    }

    ini_flag = 0;
  }

  //-----------------------------------------------------------------

  for(j=0;j<vsize;j++) lbuf[j] = 0.0;
  LOCAL_LOOP(Grid,i,j,k){
    node = Grid.node(i,j,k);
    lbuf[i] += Grid.U[node].M.x*Grid.U[node].d/hsize;  
  }
  
  MPI_Allreduce(lbuf,my_mean,vsize,MPI_DOUBLE,MPI_SUM,YZ_COMM);

  hmax  = 0.0; hphot = 0.0; fmax  = 0.0; vmax  = 0.0;
    
  for(j=Grid.lbeg[0];j <=Grid.lend[0];j++){

    hh = tau_ref*pow(my_mean[j]/(vel_ref*rho_max),4);
    hh = min(tau_min,hh);

    hmax=max(hmax,hh);
    if(rh_mean[j] < 1.e-06) hphot=max(hphot,hh);
    fmax=max(fmax,fabs(my_mean[j]));
    vmax=max(vmax,fabs(my_mean[j]/rh_mean[j]));

    my_mean[j] *= hh;
  }

  if(Run.verbose >1){
    sbuf[0] = fmax;
    sbuf[1] = vmax;
    sbuf[2] = hmax;
    sbuf[3] = hphot; 
    MPI_Reduce(sbuf,rbuf,4,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
    if (Run.rank == 0){
      cout << "damp: flx_max  = " << rbuf[0] << endl;
      cout << "damp: vel_max  = " << rbuf[1] << endl;
      cout << "damp: tau_max  = " << rbuf[2] << endl;
      cout << "damp: tau_phot = " << rbuf[3] << endl;
    }
  }
}
/*****************************************************************************/
