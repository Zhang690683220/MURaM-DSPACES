#include <mpi.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "physics.H"
#include "grid.H"
#include "run.H"
#include "limit_va.H"

#define OUTER_LOOP(G,i,j,d1,d2) \
  for((j)=(G)[(d2)][0];(j)<=(G)[(d2)][1];(j)++) \
  for((i)=(G)[(d1)][0];(i)<=(G)[(d1)][1];(i)++)

// Hall coefficient is c/(4*pi*e*n_e) , our table is ne so this is the rest.
const double hall_const = 2.99792458e10/(16.0*atan(1)*1.60217733E-19);

//-------------------------------------------------------------




/**************************************************************
                             N O T E S                         
 **************************************************************
 * Jun 12
 * From profiles when using 1 MPI rank, MHD_Residual was the
 * most time-consuming at 30% runtime.
 * Array sign[3] is constant size, may not need to copy in.
 * Need to copyin R array
 ***************************************************************/

double MHD_Residual(const RunData&  Run, GridData& Grid, 
		    const PhysicsData& Physics) {

  //double time,s_time;
  //static double t_time = 0.0 ,c_time = 0.0 , r_time = 0.0;
  //static int call_count = 0;
  //s_time = MPI_Wtime();

  static int call_count = 1;

  const double CFL_hyp  = 1.8;
  const double avr_max  = 0.25;
  const double eps_hyp  = 0.02; // -> hyperdiffusion used on heatflux, blows up during flare without it, otherwise not needed 

  const double k0_spt   = Physics.params[i_param_spitzer];

  const double Rgas = 8.314e7;
  const double X_H  = 0.7;

  // 1/6 of free streaming limit 0.25*ne*K*T*sqrt(K*T/m_e)
  const double sat_flx = (1.0+X_H)*pow(Rgas,1.5)*sqrt(1836.)/8.0;

  // corresponds to 1.5 rho cs^1.5
  //const double sat_flx  = 4.98277e12;

  cState* U  = (cState*) Grid.U;
  cState* R  = (cState*) Grid.Res;

  register int i,j,k,d,node,d1,d2,d3,ivar,offset,str;

  const double eta = Physics.params[i_param_eta];

  const bool ambipolar        = (Physics.params[i_param_ambipolar] > 0.0);
  const bool spitzer          = (Physics.params[i_param_spitzer] > 0.0);
  const bool needs_curlB      = (eta > 0.0);
  const bool need_diagnostics = Run.need_diagnostics;
  
  const int nvar   = Physics.NVAR;
  const int vsize  = Grid.vsize;
  const int v_nvar = Grid.v_nvar;
  
  double dn,lf,cs2,va2,x2,x4,s,cmax,dx,dy,dz,dxmin,dxmax,dt_res,dt_cfl,p1,p2,p3,vel,mag,ekin,pmag,mflx;
  double dt_amb,amb_diff,c_lim,cmax_hyp,amb_diff_max,amb_vel_max,nu_hyp,nu_hyp_max,rho_i,nu_in;
  double kap_spt,invBB,tt32,f_lim,F_sat,F_spt,hyp_diff;

  double amb_diff_max_acc, amb_vel_max_acc; // NEW

  double w1,w2,ww0,ww1,ww2;
  double dB2,dB3;

  int stride[3];
  int sz,gh;
 
  stride[0] = Grid.stride[0];
  stride[1] = Grid.stride[1];
  stride[2] = Grid.stride[2];

  const int loop_order[3][3] = {{ 0, 1, 2 },{ 1, 0, 2 },{ 2, 0, 1 }};
  const double del[3][3]     = {{ 1.,0.,0.},{ 0.,1.,0.},{ 0.,0.,1.}};

  const double sign[3] = {-1.,1.,-1.};

  int bounds[3][2];

  for(d=0;d<3;d++){
    bounds[d][0] = Grid.lbeg[d];
    bounds[d][1] = Grid.lend[d];
  }

  dxmin = Grid.dx[0];
  dxmax = Grid.dx[0]; 
  for(d=1;d<Grid.NDIM;d++){
    dxmin  = min(dxmin,Grid.dx[d]);
    dxmax  = max(dxmax,Grid.dx[d]);
  }

  //double var[v_nvar][vsize]; // flx[v_nvar][vsize]; //, res[v_nvar][vsize]; // curlBxB[3][vsize];
  //double lf1[vsize], msx[vsize], msy[vsize], msz[vsize]; // lrt[3][vsize], pres[vsize],
    //sflx[vsize],fcond[vsize]; //,curlB[3][vsize]; //, bb[vsize], vv[vsize], vv_amb[vsize];
    //v_amb[3][vsize]; //, amb_fac[vsize]; //, D_n[vsize], temp[vsize], BgradT[vsize], hyp_spt[vsize];

  double *curlBxB = (double*) malloc(3*vsize*sizeof(double));
  double *res     = (double*) malloc(v_nvar*vsize*sizeof(double));
  double *flx     = (double*) malloc(v_nvar*vsize*sizeof(double));
  double *var     = (double*) malloc(v_nvar*vsize*sizeof(double));
  double *D_n     = (double*) malloc(vsize*sizeof(double));
  double *temp    = (double*) malloc(vsize*sizeof(double));
  double *BgradT  = (double*) malloc(vsize*sizeof(double));
  double *hyp_spt = (double*) malloc(vsize*sizeof(double));
  double *amb_fac = (double*) malloc(vsize*sizeof(double));
  double *v_amb   = (double*) malloc(3*vsize*sizeof(double));
  double *vv_amb  = (double*) malloc(vsize*sizeof(double));
  double *vv      = (double*) malloc(vsize*sizeof(double));
  double *bb      = (double*) malloc(vsize*sizeof(double));
  double *curlB   = (double*) malloc(3*vsize*sizeof(double));
  double *fcond   = (double*) malloc(vsize*sizeof(double));
  double *sflx    = (double*) malloc(vsize*sizeof(double));
  double *pres    = (double*) malloc(vsize*sizeof(double));
  double *lrt     = (double*) malloc(3*vsize*sizeof(double));
  double *msx     = (double*) malloc(vsize*sizeof(double));
  double *msy     = (double*) malloc(vsize*sizeof(double));
  double *msz     = (double*) malloc(vsize*sizeof(double));
  double *lf1     = (double*) malloc(vsize*sizeof(double));
    
  // Estimate max characteristic velocity from previous time step (not defined at restart, use va_max instead) 
  if (Run.dt > 0.0) {
    cmax_hyp   = min(Run.CFL,CFL_hyp)*dxmin/Run.dt;
    nu_hyp_max = avr_max/Run.dt;
    hyp_diff   = eps_hyp/Run.dt;
  } else {
    cmax_hyp   = 1.0/sqrt(inv_va2max);
    nu_hyp_max = avr_max*cmax_hyp/(min(Run.CFL,CFL_hyp)*dxmin);
    hyp_diff   = eps_hyp*cmax_hyp/(min(Run.CFL,CFL_hyp)*dxmin);
  }
  
  if(needs_curlB)
    memset(Grid.curlB,0.0,Grid.bufsize*sizeof(Vector));
  if(spitzer)
    memset(Grid.BgradT,0.0,Grid.bufsize*sizeof(double));
  if(ambipolar)
    memset(Grid.curlBxB,0.0,Grid.bufsize*sizeof(Vector));
  
  cmax = 0.0;
  amb_diff_max = 0.0;
  amb_diff_max_acc = 0.0; // NEW
  amb_vel_max  = 0.0;
  amb_vel_max_acc = 0.0; // NEW
//cout << "********************************************************************" << endl;
//cout << "*                         LOOP DIMENSIONS                          *" << endl;
//cout << "********************************************************************" << endl;
  for (d=0;d<Grid.NDIM;d++){
    d1=loop_order[d][0];
    d2=loop_order[d][1];
    d3=loop_order[d][2];

    dx = Grid.dx[d1];

    w1  = 8./(12.*dx);
    w2  =-1./(12.*dx);

    ww0 = -30./(12.*dx*dx);
    ww1 =  16./(12.*dx*dx);
    ww2 =  -1./(12.*dx*dx);

    p1 = del[d1][0];
    p2 = del[d1][1];
    p3 = del[d1][2];

    sz = bounds[d1][1]-bounds[d1][0]+1+2*Grid.ghosts[d1];
    gh = Grid.ghosts[d1];

    str = stride[d1];

    // Replacing call to OUTER_LOOP with a regular double for loop
    //OUTER_LOOP(bounds,j,k,d2,d3){
//cout << "Grid " << d << endl;
//cout << "for(k=" << bounds[d3][0] << "; k<=" << bounds[d3][1] << "; k++)" << endl;
//cout << "for(j=" << bounds[d2][0] << "; j<=" << bounds[d2][1] << "; j++)" << endl;
//cout << "sz=" << sz << endl;
//cout << "gh=" << gh << endl;

    amb_diff_max = 0.0; // NEW
    amb_vel_max = 0.0; // NEW
#pragma acc parallel
{
#pragma acc loop gang collapse(2) private(offset) reduction(max:amb_diff_max, max:amb_vel_max)
    for(k=bounds[d3][0]; k<=bounds[d3][1]; k++) { // NEW
      for(j=bounds[d2][0]; j<=bounds[d2][1]; j++) { // NEW
        offset = j*stride[d2]+k*stride[d3];
        //time = MPI_Wtime();
        #pragma ivdep
#pragma acc loop vector private(node)
        for(i=0;i<sz;i++){
	        node = offset+i*str;
	        var[i*v_nvar+0] = U[node].d;
	        var[i*v_nvar+1] = U[node].M.x;
	        var[i*v_nvar+2] = U[node].M.y;
	        var[i*v_nvar+3] = U[node].M.z;
	        var[i*v_nvar+4] = U[node].e;
	        var[i*v_nvar+5] = U[node].B.x;
	        var[i*v_nvar+6] = U[node].B.y;	
	        var[i*v_nvar+7] = U[node].B.z;

	        pres[i]   = Grid.pres[node];
	        temp[i]   = Grid.temp[node];
	
	        sflx[i]   = 0.0;
	        vv_amb[i] = 0.0;
        }

        if(spitzer){
          #pragma ivdep
#pragma acc loop vector private(node)
	        for(i=0;i<sz;i++){
	          node = offset+i*str;
	          sflx[i] = Grid.sflx[node];
          }
        }

        if(ambipolar){
          #pragma ivdep
#pragma acc loop vector private(node, rho_i, nu_in)
	        for(i=0;i<sz;i++){
	          node = offset+i*str;

	          rho_i = Grid.rhoi[node];
	          nu_in = Grid.amb[node];

	          D_n[i] = max(0.0,1.0-rho_i/var[i*v_nvar+0]);
	          
	          amb_fac[i] = 1.0/(rho_i*nu_in);
	          
	          // Switch off for T>1e4 as a saveguard for now
	          amb_fac[i] *= max(0.0,min(1.0,10.0-1e-3*temp[i]));
	          
	          amb_fac[i] = min(amb_fac[i],Physics.params[i_param_ambfac_max]);

	          // test profile
	          //amb_fac[i] = 1e12*max(0.0,min(1.0,5.0-1e-3*temp[i]));
	          
	          v_amb[i*3+0] = Grid.v_amb[node].x*D_n[i];
	          v_amb[i*3+1] = Grid.v_amb[node].y*D_n[i];
	          v_amb[i*3+2] = Grid.v_amb[node].z*D_n[i];

	          vv_amb[i] = sqrt(v_amb[i*3+0]*v_amb[i*3+0]
                            +v_amb[i*3+1]*v_amb[i*3+1]
                            +v_amb[i*3+2]*v_amb[i*3+2]);
	        }
        } // end if(ambipolar)
      
        //r_time += MPI_Wtime()-time;

        //time = MPI_Wtime();
#pragma acc loop vector private(vel, mag, dn, va2, ekin, pmag, x2, x4, s, lf, cs2, cmax, mflx)
        for(i=0;i<sz;i++){ 
          res[i*v_nvar+0] = 0.0;
          res[i*v_nvar+1] = 0.0;
          res[i*v_nvar+2] = 0.0;
          res[i*v_nvar+3] = 0.0;
          res[i*v_nvar+4] = 0.0;
          res[i*v_nvar+5] = 0.0;
          res[i*v_nvar+6] = 0.0;
          res[i*v_nvar+7] = 0.0;
              
          vel = var[i*v_nvar+1+d1];
          mag = var[i*v_nvar+5+d1];

          dn   = 1.0/var[i*v_nvar+0];

          vv[i] = var[i*v_nvar+1]*var[i*v_nvar+1]+var[i*v_nvar+2]*var[i*v_nvar+2]
                + var[i*v_nvar+3]*var[i*v_nvar+3];
          bb[i] = var[i*v_nvar+5]*var[i*v_nvar+5]+var[i*v_nvar+6]*var[i*v_nvar+6]
                + var[i*v_nvar+7]*var[i*v_nvar+7];

          va2  = bb[i]*dn;
          ekin = 0.5*var[i*v_nvar+0]*vv[i];
          pmag = 0.5*bb[i];

          // approximation of 1/sqrt(1+(va/c)^4)
          x2 = va2*inv_va2max;
          x4 = x2*x2;
          s  = 1.0+x2;
          lf = s/(s+x4);

          cs2  = 5.0/3.0*pres[i]*dn;

          vv[i] = sqrt(vv[i]);
          cmax = max(cmax,vv[i]+vv_amb[i]+sqrt(max(cs2,lf*(cs2+va2))));

          lf1[i] = lf;

          mflx = vel*var[i*v_nvar+0];

          // Mass, momentum and energy fluxes without Maxwell stress
          flx[i*v_nvar+0] = mflx;
          flx[i*v_nvar+1] = mflx*var[i*v_nvar+1] + p1*pres[i];
          flx[i*v_nvar+2] = mflx*var[i*v_nvar+2] + p2*pres[i];
          flx[i*v_nvar+3] = mflx*var[i*v_nvar+3] + p3*pres[i];
          flx[i*v_nvar+4] = vel*(var[i*v_nvar+4]+ekin+pres[i]);

          //Induction  equation
          flx[i*v_nvar+5] = vel*var[i*v_nvar+5] - var[i*v_nvar+1]*mag;
          flx[i*v_nvar+6] = vel*var[i*v_nvar+6] - var[i*v_nvar+2]*mag;
          flx[i*v_nvar+7] = vel*var[i*v_nvar+7] - var[i*v_nvar+3]*mag;

          // Maxwell stress
          msx[i] = mag*var[i*v_nvar+5] - p1*pmag;
          msy[i] = mag*var[i*v_nvar+6] - p2*pmag;
          msz[i] = mag*var[i*v_nvar+7] - p3*pmag;
        }

        // Ambipolar induction
        if(ambipolar){
          for(i=0;i<sz;i++){
            flx[i*v_nvar+5] += v_amb[i*3+d1]*var[i*v_nvar+5] - v_amb[i*3+0]*var[i*v_nvar+5+d1];
            flx[i*v_nvar+6] += v_amb[i*3+d1]*var[i*v_nvar+6] - v_amb[i*3+1]*var[i*v_nvar+5+d1];
            flx[i*v_nvar+7] += v_amb[i*3+d1]*var[i*v_nvar+7] - v_amb[i*3+2]*var[i*v_nvar+5+d1];
          }
        }

        // Heat conduction
        if(spitzer){
#pragma acc loop vector
          for(i=0;i<sz;i++)
            flx[i*v_nvar+4] += sflx[i]*var[i*v_nvar+5+d1];
#pragma acc loop vector
          for(i=gh;i<sz-gh;i++){
            BgradT[i]  = var[i*v_nvar+5+d1]*(w1*(temp[i+1]-temp[i-1])+w2*(temp[i+2]-temp[i-2]));
            hyp_spt[i] = (sflx[i+2]+sflx[i-2])-4.*(sflx[i+1]+sflx[i-1])+6.*sflx[i];
          }
        }
	     
#pragma acc loop vector collapse(2)
        for (ivar=0;ivar<nvar;ivar++)
	        for(i=gh;i<sz-gh;i++)
	          res[i*v_nvar+ivar] -= 
	            w1*(flx[(i+1)*v_nvar+ivar]-flx[(i-1)*v_nvar+ivar])+
	            w2*(flx[(i+2)*v_nvar+ivar]-flx[(i-2)*v_nvar+ivar]);

        if (need_diagnostics){
	        if(spitzer){
#pragma acc loop vector
	          for(i=0;i<sz;i++)
	            fcond[i] = sflx[i]*var[i*v_nvar+5+d1];
#pragma acc loop vector private(node, dn)
	          for(i=gh;i<sz-gh;i++){
	            node = offset+i*str;
	            dn =-(w1*(fcond[i+1]-fcond[i-1])+w2*(fcond[i+2]-fcond[i-2]));
	            Grid.tvar1[node] += res[i*v_nvar+4]-dn;
	            Grid.tvar2[node] += dn;
	          }
	        } else {
#pragma acc loop vector private(node)
	          for(i=gh;i<sz-gh;i++){
	            node = offset+i*str;
	            Grid.tvar1[node] += res[i*v_nvar+4];
	          }
	        }
        } // end if(need_diagnostics)
        
        // Special treatment of Lorentz force for Boris Correction
#pragma acc loop vector private(dB2, dB3)
        for(i=gh;i<sz-gh;i++){
	        dB2 = w1*(var[(i+1)*v_nvar+5+d2]-var[(i-1)*v_nvar+5+d2])+w2*(var[(i+2)*v_nvar+5+d2]-var[(i-2)*v_nvar+5+d2]);
	        dB3 = w1*(var[(i+1)*v_nvar+5+d3]-var[(i-1)*v_nvar+5+d3])+w2*(var[(i+2)*v_nvar+5+d3]-var[(i-2)*v_nvar+5+d3]);

	        // Lorentz force
	        curlBxB[i*3+d1] =-var[i*v_nvar+5+d2]*dB2-var[i*v_nvar+5+d3]*dB3;
	        curlBxB[i*3+d2] = var[i*v_nvar+5+d1]*dB2;
	        curlBxB[i*3+d3] = var[i*v_nvar+5+d1]*dB3;

	        // curlB
	        curlB[i*3+d1] = 0;
	        curlB[i*3+d2] = sign[d1]*dB3;
	        curlB[i*3+d3]=-sign[d1]*dB2;
        }
        
#pragma acc loop vector private(lf)
        for(i=gh;i<sz-gh;i++){
	        lrt[i*3+0] = w1*(msx[i+1]-msx[i-1])+w2*(msx[i+2]-msx[i-2]);
	        lrt[i*3+1] = w1*(msy[i+1]-msy[i-1])+w2*(msy[i+2]-msy[i-2]);
	        lrt[i*3+2] = w1*(msz[i+1]-msz[i-1])+w2*(msz[i+2]-msz[i-2]);

	        // Transition from stress to jxB based on lfac
	        lf=lf1[i];
	        lrt[i*3+0] = lf*lrt[i*3+0] + (1.0-lf)*curlBxB[i*3+0];
	        lrt[i*3+1] = lf*lrt[i*3+1] + (1.0-lf)*curlBxB[i*3+1];
	        lrt[i*3+2] = lf*lrt[i*3+2] + (1.0-lf)*curlBxB[i*3+2];
	
	        res[i*v_nvar+1] += lrt[i*3+0];
	        res[i*v_nvar+2] += lrt[i*3+1];
	        res[i*v_nvar+3] += lrt[i*3+2];
	        res[i*v_nvar+4] += var[i*v_nvar+1]*lrt[i*3+0]
                          +  var[i*v_nvar+2]*lrt[i*3+1]
                          +  var[i*v_nvar+3]*lrt[i*3+2];
        }

        // Ambipolar heating + hyperdiffusion
        if(ambipolar){
#pragma acc loop vector
	        for(i=gh;i<sz-gh;i++){
	          res[i*v_nvar+4] += v_amb[i*3+0]*lrt[i*3+0]
                            +  v_amb[i*3+1]*lrt[i*3+1]
                            +  v_amb[i*3+2]*lrt[i*3+2]; 
	        }
        } 

        if (need_diagnostics){
#pragma acc loop vector private(node)
          for(i=gh;i<sz-gh;i++){
            node = offset+i*str;
            Grid.tvar3[node] += var[i*v_nvar+1]*lrt[i*3+0]
                             +  var[i*v_nvar+2]*lrt[i*3+1]
                             +  var[i*v_nvar+3]*lrt[i*3+2];
          }
          if(ambipolar)
#pragma acc loop vector private(node)
            for(i=gh;i<sz-gh;i++){
              node = offset+i*str;
              Grid.Qamb[node]  += v_amb[i*3+0]*lrt[i*3+0]
                               +  v_amb[i*3+1]*lrt[i*3+1]
                               +  v_amb[i*3+2]*lrt[i*3+2];  
            }
        } // end if(diagnostics)

        if(eta>0.0){
#pragma acc loop vector colapse(2)
          for(ivar=5;ivar<=7;ivar++){
            for(i=2;i<sz-2;i++){	    
              res[i*v_nvar+ivar] += eta*(ww0*var[i*v_nvar+ivar]+
                 ww1*(var[(i+1)*v_nvar+ivar]+var[(i-1)*v_nvar+ivar])+
                 ww2*(var[(i+2)*v_nvar+ivar]+var[(i-2)*v_nvar+ivar])
                 );	
            }
          }	      
        }

        //c_time += MPI_Wtime()-time; 

        //time = MPI_Wtime();
        #pragma ivdep
#pragma acc loop vector private(node)
        for(i=gh;i<sz-gh;i++){
          node = offset+i*str;
          R[node].d   += res[i*v_nvar+0];
          R[node].M.x += res[i*v_nvar+1];
          R[node].M.y += res[i*v_nvar+2];
          R[node].M.z += res[i*v_nvar+3];
          R[node].e   += res[i*v_nvar+4];
          R[node].B.x += res[i*v_nvar+5];
          R[node].B.y += res[i*v_nvar+6];
          R[node].B.z += res[i*v_nvar+7];
        }

        // use this to temporarily store stuff while building up jxB
        if(ambipolar){
          #pragma ivdep
#pragma acc loop vector private(node)
          for(i=gh;i<sz-gh;i++){
            node = offset+i*str;
            Grid.curlBxB[node].x += lrt[i*3+0];
            Grid.curlBxB[node].y += lrt[i*3+1];
            Grid.curlBxB[node].z += lrt[i*3+2];
          }
        }

        if(spitzer){
          #pragma ivdep
#pragma acc loop vector private(node)
          for(i=gh;i<sz-gh;i++){
            node = offset+i*str;
            Grid.BgradT[node] += BgradT[i];
            Grid.Rflx[node]   -= hyp_diff*hyp_spt[i];
          }
        }
        
        if(needs_curlB){
          #pragma ivdep
#pragma acc loop vector private(node)
          for(i=gh;i<sz-gh;i++){
            node = offset+i*str;
            Grid.curlB[node].x += curlB[i*3+0];
            Grid.curlB[node].y += curlB[i*3+1];
            Grid.curlB[node].z += curlB[i*3+2];
          }
        }

        // Last sweep, now j, jxB, BgradT are fully computed
        if(d == Grid.NDIM-1){
          // Ohmic heating
          if(eta > 0){
                  #pragma ivdep
#pragma acc loop vector private(node)
            for(i=gh;i<sz-gh;i++){
              node = offset+i*str;
              R[node].e += eta*(Grid.curlB[node].x*Grid.curlB[node].x+
		                Grid.curlB[node].y*Grid.curlB[node].y+
		                Grid.curlB[node].z*Grid.curlB[node].z);	   
            }
          }

        	// now curlBxB is fully computed, compute additional terms in Grid.R_amb
          if(ambipolar){
            #pragma ivdep
#pragma acc loop vector private(node, amb_diff, c_lim, nu_hyp) reduction(max:amb_diff_max, max:amb_vel_max)
            for(i=gh;i<sz-gh;i++){
              node = offset+i*str;

              amb_diff     = max(1e-20,D_n[i]*amb_fac[i]*bb[i]);	    
              amb_diff_max = max(amb_diff_max,amb_diff);

              amb_vel_max = max(amb_vel_max,vv_amb[i]);
              
              c_lim  = cmax_hyp - vv[i];
              nu_hyp = min(c_lim*c_lim/amb_diff,nu_hyp_max);
              
              Grid.R_amb[node] += nu_hyp*(amb_fac[i]*Grid.curlBxB[node]-Grid.v_amb[node]);
            }
          }

          if(spitzer){
            #pragma ivdep
#pragma acc loop vector private(node, c_lim, tt32, kap_spt, F_sat, invBB, F_spt, f_lim, nu_hyp)
            for(i=gh;i<sz-gh;i++){
              node = offset+i*str;
            
              c_lim    = cmax_hyp - vv[i];
              tt32     = temp[i]*sqrt(temp[i]);
              kap_spt  = k0_spt*temp[i]*tt32;
              F_sat    = sat_flx*var[i*v_nvar+0]*tt32;  
              invBB    = 1.0/max(1e-100,bb[i]);
              F_spt    =-kap_spt*Grid.BgradT[node];
              f_lim    = F_sat/(F_sat+fabs(F_spt*sqrt(invBB)));
              F_spt   *= invBB*f_lim;
              kap_spt *= f_lim;
              
              nu_hyp = min(c_lim*c_lim*var[i*v_nvar+4]/(temp[i]*kap_spt),nu_hyp_max);

              Grid.Rflx[node] += nu_hyp*(F_spt-Grid.sflx[node]);
            }
          }
        }
        //r_time += MPI_Wtime()-time;
      } // end for j
    } // end for k

} // END ACC PARALLEL REGION

// NEED TO MAKE SURE OPENACC IS PROPERLY GETTING AMB_DIFF_MAX TO HOST
  amb_diff_max_acc = max(amb_diff_max_acc, amb_diff_max); // NEW
  amb_vel_max_acc = max(amb_vel_max_acc, amb_vel_max); // NEW

  }// end for d

  amb_diff_max = amb_diff_max_acc; // NEW
  amb_vel_max = amb_vel_max_acc; // NEW
//cout << "********************************************************************" << endl;
//cout << "*                      END LOOP DIMENSIONS                         *" << endl;
//cout << "********************************************************************" << endl;

  dt_cfl = dxmin/cmax;
  cout << dt_cfl << endl;

  if(ambipolar && !(call_count%4)){
    dt_amb = Physics.params[i_param_ambipolar]*dxmin*dxmin/(amb_diff_max)/Run.CFL;
    dt_cfl = min(dt_cfl,dt_amb);

    if(Run.verbose > 3){
      double sbuf[2],rbuf[2];
      sbuf[0] = amb_diff_max;
      sbuf[1] = amb_vel_max;
      MPI_Reduce(sbuf,rbuf,2,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
      if( Run.rank == 0 ){
        cout << " Amb_diff_max = " << rbuf[0]      << " cm^2/s"         << endl;
        cout << " Amb_vel_max  = " << rbuf[1]*1e-5 << " km/s"           << endl;
        cout << " Amb_speedup  = " << rbuf[0]*Run.dt/(dxmin*dxmin) << endl;;
      }
    }
  }
  
  if (eta > 0.0) {
    if(Grid.NDIM == 1){
      dx=Grid.dx[0];
      dt_res = 0.5*dx*dx / eta;
    } else if(Grid.NDIM == 2){
      dx=Grid.dx[0];
      dy=Grid.dx[1];
      dt_res = 0.5/(1./(dx*dx)+1./(dy*dy)) / eta;
    } else {
      dx=Grid.dx[0];
      dy=Grid.dx[1];
      dz=Grid.dx[2]; 
      dt_res = 0.5/(1./(dx*dx)+1./(dy*dy)+1./(dz*dz)) / eta;
    }
    dt_res /= Run.CFL;
    
    dt_cfl = min(dt_cfl,dt_res);
  }

  //t_time += MPI_Wtime()-s_time;
  //call_count += 1;

  //if(Run.rank == 0)
  //cout << "MHD: " << t_time/call_count << ' ' 
  // << c_time/call_count << ' ' 
  // << r_time/call_count << endl;

  free(curlBxB);
  free(res);
  free(flx);
  free(var);
  free(D_n);
  free(temp);
  free(BgradT);
  free(hyp_spt);
  free(amb_fac);
  free(v_amb);
  free(vv_amb);
  free(vv);
  free(bb);
  free(curlB);
  free(fcond);
  free(sflx);
  free(pres);
  free(lrt);
  free(msx);
  free(msy);
  free(msz);
  free(lf1);

  call_count += 1;

  return dt_cfl;
}







