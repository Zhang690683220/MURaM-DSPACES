#include <mpi.h>
#include "physics.H"
#include "grid.H"
#include "run.H"
#include <math.h>
#include <stdlib.h>
#include "eos.H"
#include "limit_va.H"

/*
   - Special Corona version. Allows to set numerical pm in Corona. 
   - Resistive heating in top layer can be omitted.
   - Enhance viscosity for v>vmax_lim.
   - Output of Qres and Qvis for diagnostic reasons.
*/ 

double * hfb;
double * hft;
double * qft;

#define OUTER_LOOP(G,i,j,d1,d2) \
  for((j)=(G)[(d2)][0];(j)<=(G)[(d2)][1];(j)++) \
  for((i)=(G)[(d1)][0];(i)<=(G)[(d1)][1];(i)++)

/********************************************************************/
void TVDlimit(const RunData&  Run, GridData& Grid, 
	      const PhysicsData& Physics, const double dt_tvd){

  //double time,s_time;
  //static double t_time = 0.0 ,c_time = 0.0 , r_time = 0.0;
  //static int call_count = 0;
  //s_time = MPI_Wtime();

  const double rho_min = Physics.tchk[i_tchk_rho_min];
  const double eps_min = Physics.tchk[i_tchk_eps_min];

  const double* tvd_h    = Physics.tvd_h;
  const double* tvd_cs   = Physics.tvd_cs;  
  const double rho_lev   = pow(Physics.tvd[i_tvd_rholev],2);
  const int    rho_log   = (int) Physics.tvd[i_tvd_rholog];
  const double q_rho_max = log(Physics.tvd[i_tvd_qrho]);

  const double visc_coeff_bot = min(1.0,Physics.tvd_visc_bnd[0]);
  const double visc_slope_bot = max(0.0,Physics.tvd_visc_bnd[0]-1.0);
  const double visc_coeff_top = min(1.0,Physics.tvd_visc_bnd[1]);
  const double visc_slope_top = max(0.0,Physics.tvd_visc_bnd[1]-1.0);

  const double eta_coeff_bot = min(1.0,Physics.tvd_eta_bnd[0]);
  const double eta_slope_bot = max(0.0,Physics.tvd_eta_bnd[0]-1.0);
  const double eta_coeff_top = min(1.0,Physics.tvd_eta_bnd[1]);
  const double eta_slope_top = max(0.0,Physics.tvd_eta_bnd[1]-1.0);

  const double B_par_diff = Physics.tvd[i_tvd_Bpar];
  const double vhyp       = Physics.tvd[i_tvd_vhyp]; // Additional Hyperdiffusion in vertical direction
  const double Qdiff_bnd  = Physics.tvd[i_tvd_Qdiff_bnd]; // Value of Qdiff at upper boundary
  const double tvd_pm_v   = Physics.tvd[i_tvd_pm_v]; //
  const double tvd_pm_B   = Physics.tvd[i_tvd_pm_B]; //
  const double vmax_lim   = Physics.tvd[i_tvd_vmax_lim]; // (<1) relative to vlim, (>1) absolut, (0) disable
  const double CME_thresh = Physics.tvd[i_tvd_CME_thresh];

  const int nvar = Physics.NVAR;
  const int vsize = Grid.vsize;
  const int v_nvar = Grid.v_nvar;
 

  double h_bnd;
  if(Physics.tvd_h_bnd[0] < 1)
    h_bnd = Physics.tvd_h_bnd[0];
  else
    h_bnd = Physics.tvd_h_bnd[0]/double(Grid.gsize[0]);

  const double h_bnd_bot = h_bnd;

  if(Physics.tvd_h_bnd[1] < 1)
    h_bnd = Physics.tvd_h_bnd[1];
  else
    h_bnd = Physics.tvd_h_bnd[1]/double(Grid.gsize[0]);

  const double h_bnd_top = h_bnd;

  static int tvd_ini_flag = 1;

  const bool need_diagnostics = Run.need_diagnostics;
  const bool ambipolar        = (Physics.params[i_param_ambipolar] > 0.0);

  register int i,j,k,node,d,d1,d2,d3,ivar,offset,str,istart;
  int sz,gh;

  int k1, k2, j1, j2; // NEW

  cState* U  = (cState*) Grid.U;
    
  int stride[3];
  int dim_order[3];
 
  int bounds[3][2];

  dim_order[0] = 0;
  dim_order[1] = 1;
  dim_order[2] = 2;
  
  stride[0] = Grid.stride[0];
  stride[1] = Grid.stride[1];
  stride[2] = Grid.stride[2];

  const int loop_order[3][3] = {{ 0, 1, 2 },{ 1, 0, 2 },{ 2, 0, 1 }};

  double dn,vv,bb,cs2,va2,lf,rf,hh,cf,cmax,sl,sr,sl_lim,vsqr_diff,lower,upper,
    slm,x2,x4,s,cfast,CME_mode;

  double idx[3],tvd_coeff[4];

  double idt_full = 1.0/Run.dt;
  double dt_fac = dt_tvd/Run.dt;

  idx[0] = 1./Grid.dx[0];
  if (Grid.NDIM >= 2) idx[1] = 1./Grid.dx[1];
  if (Grid.NDIM == 3) idx[2] = 1.0/Grid.dx[2];

  tvd_coeff[0] = Physics.tvd_coeff[0];
  tvd_coeff[1] = Physics.tvd_coeff[1];
  tvd_coeff[2] = Physics.tvd_coeff[2];
  tvd_coeff[3] = Physics.tvd_coeff[3];

  int need_tvd_coeff;
  if( tvd_coeff[0]*tvd_coeff[1]*tvd_coeff[2]*tvd_coeff[3] != 1.0 )
    need_tvd_coeff = 1;
  else
    need_tvd_coeff = 0;

  int is_gend = Grid.is_gend[0]; // NEW
  int is_gbeg = Grid.is_gbeg[0]; // NEW

  if (tvd_ini_flag == 1){
    hfb = new double [vsize];
    qft = new double [vsize];
    hft = new double [vsize];

     // parabolic profile, zero for h>2*h_bnd
    for(i=Grid.lbeg[0]-Grid.ghosts[0];i<=Grid.lend[0]+Grid.ghosts[0];i++) {
      hh      = Grid.coord(i,0)/Grid.gxmax[0];
      vv = max(0.0,1.0-0.5*(1.0-hh)/h_bnd_top);
      hft[i]  = min(1.0,vv*vv);  
      vv = max(0.0,1.0-0.5*hh/h_bnd_bot);     
      hfb[i]  = min(1.0,vv*vv);

      if(hh < 1-2*h_bnd_top)
	      qft[i] = 1.0;
      else
	      qft[i] = Qdiff_bnd;
    }

    if (Run.rank == 0) {
      cout << " *** SPLIT TVD Coronal VERSION *** "           << endl;    
      if(rho_log){
	cout << " *** use log(rho) and log(eps) *** "    << endl;
      }
      cout << "tvd: apply additional 4th order hyperdiff to 0,2,4" << endl;
      cout << "tvd: special CME mode for v > CME_thresh*c_fast   " << endl;
      cout << "tvd: disable mass diffusion correction in corona  " << endl;
      cout << "tvd: change diffusivity based on lfac and rho     " << endl; 
      cout << "tvd: h              = "  << tvd_h[0] << ' ' 
	   << tvd_h[1] << ' ' << tvd_h[2] << ' ' << tvd_h[3] << endl; 
      cout << "tvd: cs             = "  << tvd_cs[0]  << ' ' 
	   << tvd_cs[1] << ' ' << tvd_cs[2] << ' ' << tvd_cs[3] << endl;
      cout << "tvd: tvd_coeff      = "
	   << tvd_coeff[0] << ' ' << tvd_coeff[1] << ' ' 
	   << tvd_coeff[2] << ' ' << tvd_coeff[3] << endl;
      cout << "tvd: rho_lev        = "  << sqrt(rho_lev)  << endl; 
      cout << "tvd: rho_log        = "  << rho_log        << endl;
      cout << "tvd: q_rho_max      = "  << exp(q_rho_max) << endl;
      cout << "tvd: vhyp           = "  << vhyp         << endl;
      cout << "tvd: B_par_diff     = "  << B_par_diff   << endl;
      cout << "tvd: Qdiff_bnd      = "  << Qdiff_bnd    << endl;
      cout << "tvd: tvd_pm_v       = "  << tvd_pm_v     << endl;
      cout << "tvd: tvd_pm_B       = "  << tvd_pm_B     << endl;
      cout << "tvd: vmax_lim       = "  << vmax_lim     << endl;
      cout << "tvd: CME_thresh     = "  << CME_thresh   << endl;
      cout << "tvd: h_bnd          = " <<  h_bnd_bot << ' ' << h_bnd_top << endl;  
      cout << "tvd: visc_coeff_bnd = " <<  visc_coeff_bot << ' '
	   << visc_coeff_top  << endl;
      cout << "tvd: visc_slope_bnd = " <<  visc_slope_bot << ' '
	   << visc_slope_top << endl;
      cout << "tvd: eta_coeff_bnd  = " <<  eta_coeff_bot << ' '
	   << eta_coeff_top << endl;
      cout << "tvd: eta_slope_bnd  = " <<  eta_slope_bot << ' ' 
	   << eta_slope_top << endl;     
    }
      
    tvd_ini_flag =0;
  }

  if(vmax_lim < 1.0)
    vsqr_diff = pow(vmax_lim*v_lim,2.0);
  else
    vsqr_diff = pow(vmax_lim,2.0);

  sz = 1;
  for(d=0;d<3;d++){
    bounds[d][0] = Grid.lbeg[d]-Grid.ghosts[d];
    bounds[d][1] = Grid.lend[d]+Grid.ghosts[d];

    sz *= Grid.lsize[d]+2*Grid.ghosts[d];
  }
  
  //double var[v_nvar][vsize], slp[v_nvar][vsize], res[v_nvar][vsize], flx[v_nvar][vsize], 
  //  uif[v_nvar][vsize], boris[4][vsize], vv_amb[vsize]; 
  //double qrho[vsize], hft1[vsize], hfb1[vsize], cm[vsize],
  //  tvd_fac[vsize], bsqr[vsize], vsqr[vsize], BC2[vsize],pres[vsize],
  //  hyperdiff[vsize], qres[vsize], qvis[vsize], qft1[vsize],
  //  hyp_e[vsize],hyp_v[vsize],hyp_B[vsize],CZ_fac[vsize];

  double *var       = (double*) malloc(v_nvar*vsize*sizeof(double));
  double *slp       = (double*) malloc(v_nvar*vsize*sizeof(double));
  double *res       = (double*) malloc(v_nvar*vsize*sizeof(double));
  double *flx       = (double*) malloc(v_nvar*vsize*sizeof(double));
  double *uif       = (double*) malloc(v_nvar*vsize*sizeof(double));
  double *boris     = (double*) malloc(   4  *vsize*sizeof(double));
  double *vv_amb    = (double*) malloc(       vsize*sizeof(double));
  double *qrho      = (double*) malloc(       vsize*sizeof(double));
  double *hft1      = (double*) malloc(       vsize*sizeof(double));
  double *hfb1      = (double*) malloc(       vsize*sizeof(double));
  double *cm        = (double*) malloc(       vsize*sizeof(double));
  double *tvd_fac   = (double*) malloc(       vsize*sizeof(double));
  double *bsqr      = (double*) malloc(       vsize*sizeof(double));
  double *vsqr      = (double*) malloc(       vsize*sizeof(double));
  double *BC2       = (double*) malloc(       vsize*sizeof(double));
  double *pres      = (double*) malloc(       vsize*sizeof(double));
  double *hyperdiff = (double*) malloc(       vsize*sizeof(double));
  double *qres      = (double*) malloc(       vsize*sizeof(double));
  double *qvis      = (double*) malloc(       vsize*sizeof(double));
  double *qft1      = (double*) malloc(       vsize*sizeof(double));
  double *hyp_e     = (double*) malloc(       vsize*sizeof(double));
  double *hyp_v     = (double*) malloc(       vsize*sizeof(double));
  double *hyp_B     = (double*) malloc(       vsize*sizeof(double));
  double *CZ_fac    = (double*) malloc(       vsize*sizeof(double));

  /* y direction first to be consistent with vertical boundary */
  for (d=0;d<Grid.NDIM;d++){
    d1=loop_order[dim_order[d]][0];
    d2=loop_order[dim_order[d]][1];
    d3=loop_order[dim_order[d]][2];

    cmax = 0.975*Grid.dx[d1]/dt_tvd;

    sz = bounds[d1][1]-bounds[d1][0]+1;
    gh = Grid.ghosts[d1];

    str=stride[d1];

    k1 = bounds[d3][0];
    k2 = bounds[d3][1];
    j1 = bounds[d2][0];
    j2 = bounds[d2][1];


#pragma acc enter data copyin(Grid)
#pragma acc enter data copyin(stride[0:3], U[0:Grid.bufsize], Grid.Tau[0:Grid.bufsize], \
  Grid.v_amb[0:Grid.bufsize], hfb[0:vsize], qft[0:vsize], hft[0:vsize], tvd_coeff[0:4], \
  Grid.tvar8[0:Grid.bufsize], Grid.Qres[0:Grid.bufsize], Grid.Qvis[0:Grid.bufsize], \
  Grid.tvar6[0:Grid.bufsize], Grid.tvar7[0:Grid.bufsize], idx[0:3], tvd_cs[0:4], tvd_h[0:4], \
  Grid.pres[0:Grid.bufsize])
#pragma acc parallel default(present)
{
#pragma acc loop collapse(2) gang independent private(offset, var[0:v_nvar*vsize], slp[0:v_nvar*vsize], \
  pres[0:vsize], CZ_fac[0:vsize], vv_amb[0:vsize], hfb1[0:vsize], hft1[0:vsize], qft1[0:vsize], \
  vsqr[0:vsize], bsqr[0:vsize], cm[0:vsize], hyp_e[0:vsize], hyp_v[0:vsize], hyp_B[0:vsize], \
  BC2[0:vsize], hyperdiff[0:vsize], qrho[0:vsize], flx[0:v_nvar*vsize], uif[0:v_nvar*vsize], \
  tvd_fac[0:vsize], res[0:v_nvar*vsize], boris[0:4*vsize], qres[0:vsize], qvis[0:vsize])
    //OUTER_LOOP(bounds,j,k,d2,d3){
    for(k = k1; k <= k2; k++) {
      for(j = j1; j <= j2; j++) {

        offset = j*stride[d2]+k*stride[d3];
        //time = MPI_Wtime();
        #pragma ivdep
#pragma acc loop vector independent private(node)
        for(i=0;i<sz;i++){
	        node = offset+i*str;
	        var[i*v_nvar+0] = max(rho_min,U[node].d);
	        var[i*v_nvar+1] = U[node].M.x;
	        var[i*v_nvar+2] = U[node].M.y;
	        var[i*v_nvar+3] = U[node].M.z;
	        var[i*v_nvar+4] = U[node].e;
	        var[i*v_nvar+5] = U[node].B.x;
	        var[i*v_nvar+6] = U[node].B.y;	
	        var[i*v_nvar+7] = U[node].B.z;
	
	        pres[i] = Grid.pres[node];
	        CZ_fac[i] = (double) (Grid.Tau[node] > 1.0e-5);

	        vv_amb[i] = 0.0;
        }

        if(ambipolar){
	        #pragma ivdep
#pragma acc loop vector independent private(node)
	        for(i=0;i<sz;i++){
	          node = offset+i*str;
	          vv_amb[i] = sqrt(Grid.v_amb[node].x*Grid.v_amb[node].x
                           + Grid.v_amb[node].y*Grid.v_amb[node].y
                           + Grid.v_amb[node].z*Grid.v_amb[node].z);
	        }
        }
	
        if(d1 == 0){
#pragma acc loop vector independent
	        for(i=0;i<sz;i++){
	          hfb1[i] = hfb[i];
	          hft1[i] = hft[i];
	          qft1[i] = qft[i];
	        }
        } else if (d2 == 0){
#pragma acc loop vector independent
	        for(i=0;i<sz;i++){
	          hfb1[i] = hfb[j];
	          hft1[i] = hft[j];
	          qft1[i] = qft[j];
	        }
        } else if(d3 == 0){
#pragma acc loop vector independent
	        for(i=0;i<sz;i++){
	          hfb1[i] = hfb[k];
	          hft1[i] = hft[k];
	          qft1[i] = qft[k];
	        }
        }
        //r_time += MPI_Wtime()-time;

        //time = MPI_Wtime();
#pragma acc loop vector independent private(dn, vv, bb)
        for(i=0;i<sz;i++){
	        dn = 1.0/var[i*v_nvar+0];
          var[i*v_nvar+1] = var[i*v_nvar+1]*dn;
	        var[i*v_nvar+2] = var[i*v_nvar+2]*dn;
	        var[i*v_nvar+3] = var[i*v_nvar+3]*dn;
	
	        vv = (var[i*v_nvar+1]*var[i*v_nvar+1]+
		            var[i*v_nvar+2]*var[i*v_nvar+2]+
		            var[i*v_nvar+3]*var[i*v_nvar+3]);
	             
	        bb = (var[i*v_nvar+5]*var[i*v_nvar+5]+
		            var[i*v_nvar+6]*var[i*v_nvar+6]+
		            var[i*v_nvar+7]*var[i*v_nvar+7]);

	        var[i*v_nvar+4] = var[i*v_nvar+4]*dn-0.5*vv;

	        vsqr[i] = vv;
	        bsqr[i] = bb;
        }
       
#pragma acc loop vector independent private(dn, vv, va2, x2, x4, s, lf, cs2, cfast, CME_mode, \
  rf, hh, cf, cs2)
        for(i=0;i<sz;i++){
	        dn   = 1.0/var[i*v_nvar+0];
	        vv   = sqrt(vsqr[i]);
	        va2  = bsqr[i]*dn;
	
	        x2 = va2*inv_va2max;
	        x4 = x2*x2;
	        s  = 1.0+x2;
	        lf = s/(s+x4);

	        cs2  = 5.0/3.0*pres[i]*dn;

	        cfast = sqrt(va2+cs2);

	        CME_mode = (1.0-CZ_fac[i])*( (double) (vv >= CME_thresh*cfast) );
	
	        dn = var[i*v_nvar+0]*var[i*v_nvar+0];
	        rf = min(lf,dn/(dn+rho_lev));
	
	        hh   = 1.0-hfb1[i]-hft1[i];
	        cf   = hfb1[i]*tvd_cs[0]+(rf*tvd_cs[1]+(1.0-rf)
                 *tvd_cs[2])*hh+hft1[i]*tvd_cs[3];

	        cf   = max(cf,CME_mode);
	
	        cs2 *= cf*cf;
	
	        cm[i]   = vv_amb[i]+vv+sqrt(max(cs2,lf*(cs2+va2)));
	
	        hyp_e[i]  = hfb1[i]*tvd_h[0]+(rf*tvd_h[1]+(1.0-rf)
                      *tvd_h[2])*hh+hft1[i]*tvd_h[3];
	        hyp_v[i]  = hfb1[i]*tvd_h[0]+(rf*tvd_h[1]+(1.0-rf)
                      *tvd_h[2]*tvd_pm_v)*hh+hft1[i]*tvd_h[3]*tvd_pm_v;
	        hyp_B[i]  = hfb1[i]*tvd_h[0]+(rf*tvd_h[1]+(1.0-rf)
                      *tvd_h[2]*tvd_pm_B)*hh+hft1[i]*tvd_h[3]*tvd_pm_B;

	        hyp_e[i] *= (1.0-CME_mode);
	        hyp_v[i] *= (1.0-CME_mode);
	
	        BC2[i]  = 1.0-lf;
        }
      
        if( (d1 == 0) && (vhyp > 0.0) ){
#pragma acc loop vector independent
	        for(i=0;i<sz;i++){
	          hyperdiff[i] = vhyp*sqrt(vsqr[i]);	  
	          hyperdiff[i] = 0.5*idx[d1]*max(0.0,min(hyperdiff[i],cmax-cm[i]));
	        }
        }

        if(rho_log){
#pragma acc loop vector independent
	        for(i=0;i<sz;i++){ 
	          var[i*v_nvar+0] = log(var[i*v_nvar+0]);
	          var[i*v_nvar+4] = log(var[i*v_nvar+4]);
	        }
#pragma acc loop vector independent
	        for(i=0;i<sz-1;i++)
	          qrho[i] = fabs(var[i*v_nvar+0]-var[(i+1)*v_nvar+0]);
        } else {
#pragma acc loop vector independent
	        for(i=0;i<sz-1;i++)
	          qrho[i] = fabs(log(var[i*v_nvar+0])-log(var[(i+1)*v_nvar+0])); 
        }
	
        // reconstruction slopes
#pragma acc loop collapse(2) vector independent private(sl, sr, slm, lower, upper)
        for(ivar=0;ivar<nvar;ivar++){
	        for(i=1;i<sz-1;i++){
	          sl    = var[i*v_nvar+ivar]-var[(i-1)*v_nvar+ivar];
	          sr    = var[(i+1)*v_nvar+ivar]-var[i*v_nvar+ivar];
	          slm   = 0.25*(sl+sr);
	          lower = min(min(sl,sr),slm);
	          upper = max(max(sl,sr),slm);
	          slp[i*v_nvar+ivar]=max(0.0,lower)+min(0.0,upper);
	        }
        }

        // more diffusivity at vertical boundaries
        if( (visc_slope_bot > 0.0) or ( visc_slope_top > 0.0) ){
#pragma acc loop vector independent private(cf)
	        for(i=1;i<sz-1;i++){
	          cf =(1.0-visc_slope_bot*hfb1[i])*(1.0-visc_slope_top*hft1[i]);
	          slp[i*v_nvar+1] *= cf; slp[i*v_nvar+2] *= cf; slp[i*v_nvar+3] *= cf;	
	        }
        }

        if( (eta_slope_bot > 0.0) or ( eta_slope_top > 0.0) ){
#pragma acc loop vector independent private(cf)
	        for(i=1;i<sz-1;i++){
	          cf =(1.0-eta_slope_bot*hfb1[i])*(1.0-eta_slope_top*hft1[i]);
	          slp[i*v_nvar+5] *= cf; slp[i*v_nvar+6] *= cf; slp[i*v_nvar+7] *= cf;	
	        }
        }

        // more diffusivity around large density jumps
#pragma acc loop vector independent
        for(i=1;i<sz-1;i++){
	        if(qrho[i] > q_rho_max)
	          slp[i*v_nvar+0] = 0.0;
        }

        // more diffusivity for very fast flows
        if(vmax_lim > 0.0){
#pragma acc loop vector independent
	        for(i=1;i<sz-1;i++){
	          if(vsqr[i] > vsqr_diff){
	            slp[i*v_nvar+1]  = 0.0;
	            slp[i*v_nvar+2]  = 0.0;
	            slp[i*v_nvar+3]  = 0.0;
	          }
	        }
        }

        // This loop is kinda weird and seems unnecessary.
        // I'm pretty sure nvar will always be size 8, and I'll test it.
        /*for(ivar=0;ivar<nvar;ivar++){
	        if( (ivar >=1) && (ivar <= 3) ){
	         for(i=1;i<sz-2;i++){
	            sl     = var[(i+1)*v_nvar+ivar]-var[i*v_nvar+ivar];
	            sl_lim = var[(i+1)*v_nvar+ivar]-var[i*v_nvar+ivar]-(slp[i*v_nvar+ivar] +slp[(i+1)*v_nvar+ivar]);
	            dn     = hyp_v[i];
	            rf     = fabs(sl_lim)/max(1e-100,fabs(sl));
	            cf     = max(0.0,min(sl_lim,sl))+min(0.0,max(sl_lim,sl));
	            flx[i*v_nvar+ivar] = max(0.0,1.0+dn*(rf-1.0))*cf;
	            uif[i*v_nvar+ivar] = 0.5*(var[i*v_nvar+ivar]+slp[i*v_nvar+ivar]+
				        var[(i+1)*v_nvar+ivar]-slp[(i+1)*v_nvar+ivar]);
	          }
	        } else if( (ivar >=5) && (ivar <= 7) ){
	          for(i=1;i<sz-2;i++){
	            sl     = var[(i+1)*v_nvar+ivar]-var[i*v_nvar+ivar];
	            sl_lim = var[(i+1)*v_nvar+ivar]-var[i*v_nvar+ivar]-(slp[i*v_nvar+ivar]+slp[(i+1)*v_nvar+ivar]);
	            dn     = hyp_B[i];
	            rf     = fabs(sl_lim)/max(1e-100,fabs(sl));
	            cf     = max(0.0,min(sl_lim,sl))+min(0.0,max(sl_lim,sl));
	            flx[i*v_nvar+ivar] = max(0.0,1.0+dn*(rf-1.0))*cf;
	            uif[i*v_nvar+ivar] = 0.5*(var[i*v_nvar+ivar]+slp[i*v_nvar+ivar]+
				        var[(i+1)*v_nvar+ivar]-slp[(i*1)*v_nvar+ivar]);
	          }
	        } else {
	          for(i=1;i<sz-2;i++){
	            sl     = var[(i+1)*v_nvar+ivar]-var[i*v_nvar+ivar];
	            sl_lim = var[(i+1)*v_nvar+ivar]-var[i*v_nvar+ivar]-(slp[i*v_nvar+ivar]+slp[(i+1)*v_nvar+ivar]);
	           dn     = hyp_e[i];
	            rf     = fabs(sl_lim)/max(1e-100,fabs(sl));
	            cf     = max(0.0,min(sl_lim,sl))+min(0.0,max(sl_lim,sl));
	            flx[i*v_nvar+ivar] = max(0.0,1.0+dn*(rf-1.0))*cf;
	            uif[i*v_nvar+ivar] = 0.5*(var[i*v_nvar+ivar]+slp[i*v_nvar+ivar]+
				        var[(i+1)*v_nvar+ivar]-slp[(i+1)*v_nvar+ivar]);
	          }
	       }
        }*/

// Experimental change -- probably okay but does need to be double checked
#pragma acc loop vector independent private(sl, sl_lim, dn, rf, cf)
        for(i=1;i<sz-2;i++){
          // ivar = 0
      	  sl     = var[(i+1)*v_nvar+0]-var[i*v_nvar+0];
          sl_lim = var[(i+1)*v_nvar+0]-var[i*v_nvar+0]-(slp[i*v_nvar+0] +slp[(i+1)*v_nvar+0]);
          dn     = hyp_e[i];
          rf     = fabs(sl_lim)/max(1e-100,fabs(sl));
          cf     = max(0.0,min(sl_lim,sl))+min(0.0,max(sl_lim,sl));
          flx[i*v_nvar+0] = max(0.0,1.0+dn*(rf-1.0))*cf;
          uif[i*v_nvar+0] = 0.5*(var[i*v_nvar+0]+slp[i*v_nvar+0]+
		        var[(i+1)*v_nvar+0]-slp[(i+1)*v_nvar+0]);
          // ivar = 1
      	  sl     = var[(i+1)*v_nvar+1]-var[i*v_nvar+1];
          sl_lim = var[(i+1)*v_nvar+1]-var[i*v_nvar+1]-(slp[i*v_nvar+1] +slp[(i+1)*v_nvar+1]);
          dn     = hyp_v[i];
          rf     = fabs(sl_lim)/max(1e-100,fabs(sl));
          cf     = max(0.0,min(sl_lim,sl))+min(0.0,max(sl_lim,sl));
          flx[i*v_nvar+1] = max(0.0,1.0+dn*(rf-1.0))*cf;
          uif[i*v_nvar+1] = 0.5*(var[i*v_nvar+1]+slp[i*v_nvar+1]+
		        var[(i+1)*v_nvar+1]-slp[(i+1)*v_nvar+1]);
          // ivar = 2
      	  sl     = var[(i+1)*v_nvar+2]-var[i*v_nvar+2];
          sl_lim = var[(i+1)*v_nvar+2]-var[i*v_nvar+2]-(slp[i*v_nvar+2] +slp[(i+1)*v_nvar+2]);
          dn     = hyp_v[i];
          rf     = fabs(sl_lim)/max(1e-100,fabs(sl));
          cf     = max(0.0,min(sl_lim,sl))+min(0.0,max(sl_lim,sl));
          flx[i*v_nvar+2] = max(0.0,1.0+dn*(rf-1.0))*cf;
          uif[i*v_nvar+2] = 0.5*(var[i*v_nvar+2]+slp[i*v_nvar+2]+
		        var[(i+1)*v_nvar+2]-slp[(i+1)*v_nvar+2]);
          // ivar = 3
      	  sl     = var[(i+1)*v_nvar+3]-var[i*v_nvar+3];
          sl_lim = var[(i+1)*v_nvar+3]-var[i*v_nvar+3]-(slp[i*v_nvar+3] +slp[(i+1)*v_nvar+3]);
          dn     = hyp_v[i];
          rf     = fabs(sl_lim)/max(1e-100,fabs(sl));
          cf     = max(0.0,min(sl_lim,sl))+min(0.0,max(sl_lim,sl));
          flx[i*v_nvar+3] = max(0.0,1.0+dn*(rf-1.0))*cf;
          uif[i*v_nvar+3] = 0.5*(var[i*v_nvar+3]+slp[i*v_nvar+3]+
		        var[(i+1)*v_nvar+3]-slp[(i+1)*v_nvar+3]);
          // ivar = 4
      	  sl     = var[(i+1)*v_nvar+4]-var[i*v_nvar+4];
          sl_lim = var[(i+1)*v_nvar+4]-var[i*v_nvar+4]-(slp[i*v_nvar+4] +slp[(i+1)*v_nvar+4]);
          dn     = hyp_e[i];
          rf     = fabs(sl_lim)/max(1e-100,fabs(sl));
          cf     = max(0.0,min(sl_lim,sl))+min(0.0,max(sl_lim,sl));
          flx[i*v_nvar+4] = max(0.0,1.0+dn*(rf-1.0))*cf;
          uif[i*v_nvar+4] = 0.5*(var[i*v_nvar+4]+slp[i*v_nvar+4]+
		        var[(i+1)*v_nvar+4]-slp[(i+1)*v_nvar+4]);
          // ivar = 5
      	  sl     = var[(i+1)*v_nvar+5]-var[i*v_nvar+5];
          sl_lim = var[(i+1)*v_nvar+5]-var[i*v_nvar+5]-(slp[i*v_nvar+5] +slp[(i+1)*v_nvar+5]);
          dn     = hyp_B[i];
          rf     = fabs(sl_lim)/max(1e-100,fabs(sl));
          cf     = max(0.0,min(sl_lim,sl))+min(0.0,max(sl_lim,sl));
          flx[i*v_nvar+5] = max(0.0,1.0+dn*(rf-1.0))*cf;
          uif[i*v_nvar+5] = 0.5*(var[i*v_nvar+5]+slp[i*v_nvar+5]+
		        var[(i+1)*v_nvar+5]-slp[(i+1)*v_nvar+5]);
          // ivar = 6
      	  sl     = var[(i+1)*v_nvar+6]-var[i*v_nvar+6];
          sl_lim = var[(i+1)*v_nvar+6]-var[i*v_nvar+6]-(slp[i*v_nvar+6] +slp[(i+1)*v_nvar+6]);
          dn     = hyp_B[i];
          rf     = fabs(sl_lim)/max(1e-100,fabs(sl));
          cf     = max(0.0,min(sl_lim,sl))+min(0.0,max(sl_lim,sl));
          flx[i*v_nvar+6] = max(0.0,1.0+dn*(rf-1.0))*cf;
          uif[i*v_nvar+6] = 0.5*(var[i*v_nvar+6]+slp[i*v_nvar+6]+
		        var[(i+1)*v_nvar+6]-slp[(i+1)*v_nvar+6]);
          // ivar = 7
      	  sl     = var[(i+1)*v_nvar+7]-var[i*v_nvar+7];
          sl_lim = var[(i+1)*v_nvar+7]-var[i*v_nvar+7]-(slp[i*v_nvar+7] +slp[(i+1)*v_nvar+7]);
          dn     = hyp_B[i];
          rf     = fabs(sl_lim)/max(1e-100,fabs(sl));
          cf     = max(0.0,min(sl_lim,sl))+min(0.0,max(sl_lim,sl));
          flx[i*v_nvar+7] = max(0.0,1.0+dn*(rf-1.0))*cf;
          uif[i*v_nvar+7] = 0.5*(var[i*v_nvar+7]+slp[i*v_nvar+7]+
		        var[(i+1)*v_nvar+7]-slp[(i+1)*v_nvar+7]);
        }



        // diffusion coefficient
#pragma acc loop vector independent
        for(i=0;i<sz-1;i++)
	        tvd_fac[i] = 0.5*idx[d1]*min(cmax,max(cm[i],cm[i+1]));
      
#pragma acc loop collapse(2) vector independent
        for(ivar=0;ivar<nvar;ivar++)
	        for(i=1;i<sz-2;i++)
	          flx[i*v_nvar+ivar] *= tvd_fac[i];
      
        // avoid terms that produce large divB error -> dBx/dx, dBy/dy, dBz/dz
#pragma acc loop vector independent
        for(i=1;i<sz-2;i++)
	        flx[i*v_nvar+5+d1] *= B_par_diff;

        if(need_tvd_coeff){
#pragma acc loop vector independent private(cf)
	        for(i=1;i<sz-2;i++){
	          flx[i*v_nvar+0] *= tvd_coeff[0];
	          
	          cf = max(tvd_coeff[1],visc_coeff_bot*hfb1[i]+visc_coeff_top*hft1[i]);
	          flx[i*v_nvar+1] *= cf;
	          flx[i*v_nvar+2] *= cf;
	          flx[i*v_nvar+3] *= cf;
	          
	          flx[i*v_nvar+4] *= tvd_coeff[2];
	          
	          cf = max(tvd_coeff[3],eta_coeff_bot*hfb1[i]+eta_coeff_top*hft1[i]);
	          flx[i*v_nvar+5] *= cf;
	          flx[i*v_nvar+6] *= cf;
	          flx[i*v_nvar+7] *= cf; 
	        }
        } 

        // additional 4th order hyperdiffusion for HD varibales in the
        // y-direction to damp spurious oscillations
        if( (d1 == 0) && (vhyp > 0.0) ){
	        istart = 1;
	        if(Grid.is_gbeg[0]) istart += gh;
#pragma acc loop vector independent private(dn, ivar)
	        for(i=istart;i<sz-2;i++){
	          dn = max(hyperdiff[i],hyperdiff[i+1]);
	          ivar = 0;
	          flx[i*v_nvar+ivar] -= dn*(0.125*(var[(i+2)*v_nvar+ivar]-var[(i-1)*v_nvar+ivar]) -
			              0.375*(var[(i+1)*v_nvar+ivar]-var[i*v_nvar+ivar]));
	          ivar = 1;
	          flx[i*v_nvar+ivar] -= dn*(0.125*(var[(i+2)*v_nvar+ivar]-var[(i-1)*v_nvar+ivar]) -
			              0.375*(var[(i+1)*v_nvar+ivar]-var[i*v_nvar+ivar]));
	          ivar = 4;
	          flx[i*v_nvar+ivar] -= dn*(0.125*(var[(i+2)*v_nvar+ivar]-var[(i-1)*v_nvar+ivar]) -
			              0.375*(var[(i+1)*v_nvar+ivar]-var[i*v_nvar+ivar]));
	        }
        }
      
        if(rho_log){
#pragma acc loop vector independent
	        for(i=1;i<sz-2;i++){   
	          uif[i*v_nvar+0]  = exp(uif[i*v_nvar+0]);
	          flx[i*v_nvar+0] *= uif[i*v_nvar+0];

	          uif[i*v_nvar+4]  = exp(uif[i*v_nvar+4]);
	          flx[i*v_nvar+4] *= uif[i*v_nvar+4];
	        }
        }

#pragma acc loop collapse(2) vector independent
        for(ivar=1;ivar<=4;ivar++){
	        for(i=1;i<sz-2;i++){
	          flx[i*v_nvar+ivar] *= uif[i*v_nvar+0];
	        }
        }

        // correction of momentum and energy flux for mass diffusion
#pragma acc loop collapse(2) vector independent
        for(ivar=1;ivar<=4;ivar++){
	        for(i=1;i<sz-2;i++){
	          flx[i*v_nvar+ivar] += uif[i*v_nvar+ivar]*flx[i*v_nvar+0]*min(CZ_fac[i],CZ_fac[i+1]);
	        }
        } 
	  
        // viscous energy flux
#pragma acc loop vector independent
        for(i=1;i<sz-2;i++)   
	        flx[i*v_nvar+4] += qft1[i]*(uif[i*v_nvar+1]*flx[i*v_nvar+1]+uif[i*v_nvar+2]*flx[i*v_nvar+2]+uif[i*v_nvar+3]*flx[i*v_nvar+3]);

        // no diffusion for vertical field at boundaries (conserve vertical magnetic flux)
        if( d1 == 0) {
	        //if( Grid.is_gend[0] ) flx[(sz-gh-1)*v_nvar+5] = 0.0;
	        //if( Grid.is_gbeg[0] ) flx[(gh-1)*v_nvar+5] = 0.0;
	        if(is_gend) flx[(sz-gh-1)*v_nvar+5] = 0.0;
	        if(is_gbeg) flx[(gh-1)*v_nvar+5] = 0.0;
        } 
      
#pragma acc loop collapse(2) vector independent
        for(ivar=0;ivar<nvar;ivar++){
	        for(i=gh;i<sz-gh;i++){
	          res[i*v_nvar+ivar] = dt_tvd*(flx[i*v_nvar+ivar]-flx[(i-1)*v_nvar+ivar]);
	        }
        }

        // Remove viscous heating at top boundary
#pragma acc loop vector independent
        for(i=gh;i<sz-gh;i++)
	        res[i*v_nvar+4] += (1.0-qft1[i])*(res[i*v_nvar+1]*var[i*v_nvar+1] + res[i*v_nvar+2]*var[i*v_nvar+2] + res[i*v_nvar+3]*var[i*v_nvar+3]);

        if(need_diagnostics){
          #pragma ivdep
#pragma acc loop vector independent private(node)
	        for(i=gh;i<sz-gh;i++){
	          node = offset+i*str;
	          Grid.tvar8[node] += res[i*v_nvar+4]*idt_full;
	        }
        }

        // projection of fvisc for consistency with SR treatment
#pragma acc loop vector independent private(dn)
        for(i=gh;i<sz-gh;i++){
	        dn = (res[i*v_nvar+1]*var[i*v_nvar+5] + res[i*v_nvar+2]*var[i*v_nvar+6] + res[i*v_nvar+3]*var[i*v_nvar+7])/max(1e-100,bsqr[i]);

	        boris[i*4+0] = -BC2[i]*(res[i*v_nvar+1]-dn*var[i*v_nvar+5]);
	        boris[i*4+1] = -BC2[i]*(res[i*v_nvar+2]-dn*var[i*v_nvar+6]);
	        boris[i*4+2] = -BC2[i]*(res[i*v_nvar+3]-dn*var[i*v_nvar+7]);
	        boris[i*4+3] =  boris[i*4+0]*var[i*v_nvar+1]+boris[i*4+1]*var[i*v_nvar+2]+boris[i*4+2]*var[i*v_nvar+3];	
        }

#pragma acc loop vector independent
        for(i=gh;i<sz-gh;i++){
	        res[i*v_nvar+1] += boris[i*4+0];
	        res[i*v_nvar+2] += boris[i*4+1];
	        res[i*v_nvar+3] += boris[i*4+2];
	        res[i*v_nvar+4] += boris[i*4+3];
        }
      
        // resistive heating 
#pragma acc loop vector independent
        for(i=1;i<sz-2;i++){
	        qrho[i] = 
	          (var[(i+1)*v_nvar+5]-var[i*v_nvar+5])*flx[i*v_nvar+5] +
	          (var[(i+1)*v_nvar+6]-var[i*v_nvar+6])*flx[i*v_nvar+6] +
	          (var[(i+1)*v_nvar+7]-var[i*v_nvar+7])*flx[i*v_nvar+7];
        }
      
#pragma acc loop vector independent
        for(i=gh;i<sz-gh;i++){
	        qres[i] = 0.5*(qrho[i]+qrho[i-1])*qft1[i];
	        res[i*v_nvar+4] += dt_tvd*qres[i];
        }
        
        // viscous heating for diagnostics
        if(need_diagnostics){
#pragma acc loop vector independent
	        for(i=1;i<sz-2;i++){
	          qrho[i] = 
	            (var[(i+1)*v_nvar+1]-var[i*v_nvar+1])*flx[i*v_nvar+1] +
	            (var[(i+1)*v_nvar+2]-var[i*v_nvar+2])*flx[i*v_nvar+2] +
	            (var[(i+1)*v_nvar+3]-var[i*v_nvar+3])*flx[i*v_nvar+3];
	        }

#pragma acc loop vector independent
	        for(i=gh;i<sz-gh;i++){
	          qvis[i] = 0.5*(qrho[i]+qrho[i-1])*qft1[i];
	        }
        }
 
        //c_time += MPI_Wtime()-time;

        //time = MPI_Wtime();
        #pragma ivdep
#pragma acc loop vector independent private(node)
        for(i=gh;i<sz-gh;i++){
	        node = offset+i*str;
	        U[node].d   += res[i*v_nvar+0];
	        U[node].M.x += res[i*v_nvar+1];
	        U[node].M.y += res[i*v_nvar+2];
	        U[node].M.z += res[i*v_nvar+3];
	        U[node].e   += res[i*v_nvar+4];
	        U[node].B.x += res[i*v_nvar+5];
	        U[node].B.y += res[i*v_nvar+6];
	        U[node].B.z += res[i*v_nvar+7];
	        U[node].d    = max(rho_min,U[node].d);
        }

        if (need_diagnostics){ 
          #pragma ivdep
#pragma acc loop vector independent private(node)
	        for(i=gh;i<sz-gh;i++){
	          node = offset+i*str;
	          Grid.Qres[node] += dt_fac*qres[i];
	          Grid.Qvis[node] += dt_fac*qvis[i];
	
	          Grid.tvar6[node] += boris[i*4+3]*idt_full;
	        }
        }

        #pragma ivdep
#pragma acc loop vector independent private(node)
        for(i=gh;i<sz-gh;i++){
	        node = offset+i*str;
	        bsqr[i] = U[node].e;
	        U[node].e = max(U[node].e,eps_min*U[node].d+0.5*U[node].M.sqr()/U[node].d);
        }

        if (need_diagnostics){ 
          #pragma ivdep
#pragma acc loop vector independent private(node)
	        for(i=gh;i<sz-gh;i++){
	          node = offset+i*str;     
	          Grid.tvar7[node] += (U[node].e-bsqr[i])*idt_full;
	        }
        }
        //r_time += MPI_Wtime()-time;
      } // end J loop
    } // end K loop

} // END PARALLEL

#pragma acc exit data copyout(stride[0:3], U[0:Grid.bufsize], Grid.Tau[0:Grid.bufsize], \
  Grid.v_amb[0:Grid.bufsize], hfb[0:vsize], qft[0:vsize], hft[0:vsize], tvd_coeff[0:4], \
  Grid.tvar8[0:Grid.bufsize], Grid.Qres[0:Grid.bufsize], Grid.Qvis[0:Grid.bufsize], \
  Grid.tvar6[0:Grid.bufsize], Grid.tvar7[0:Grid.bufsize], idx[0:3], tvd_cs[0:4], tvd_h[0:4], \
  Grid.pres[0:Grid.bufsize])
#pragma acc exit data delete(Grid)


    bounds[d1][0] = Grid.lbeg[d1];
    bounds[d1][1] = Grid.lend[d1];
  } // end D loop

  //t_time += MPI_Wtime()-s_time;
  //call_count += 1;

  // if(Run.rank == 0)
  // cout << "TVD: " << t_time/call_count << ' ' 
  // << c_time/call_count << ' ' 
  // << r_time/call_count << endl; 

  free(var);
  free(slp);
  free(res);
  free(flx);
  free(uif);
  free(boris);
  free(vv_amb);
  free(qrho);
  free(hft1);
  free(hfb1);
  free(cm);
  free(tvd_fac);
  free(bsqr);
  free(vsqr);
  free(BC2);
  free(pres);
  free(hyperdiff);
  free(qres);
  free(qvis);
  free(qft1);
  free(hyp_e);
  free(hyp_v);
  free(hyp_B);
  free(CZ_fac);


}


