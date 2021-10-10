#include <mpi.h>
#include <fstream>
#include <stdlib.h>
#include "analysis.H"
#include "physics.H"
#include "grid.H"
#include "run.H"
#include "rt/rt.h"
#include "comm_split.H"
#include "limit_va.H"
#include "io.H"

using namespace std;

extern est_total_slice_iters;
extern struct log *io_file_log, *io_dspaces_log;

float *analyzevp_buf = NULL;
int analyzevp_nvar;
dspaces_put_req_t* analyzevp_dspaces_put_req_list = NULL;

//======================================================================
void AnalyzeSolution_VP(const RunData& Run,const GridData& Grid,
		        const PhysicsData& Physics, RTS * rts) {

  static int ini_flag = 1;
  static int analyzevp_ref_count = 0;

  static MPI_Datatype x_subarray;

  

  register int i,j,k,node,ind,ioff,v;

  int next[3];
  double w1[3],w2[3];
  
  for(v=0;v<Grid.NDIM;v++){
    next[v] = Grid.stride[v];
    w1[v]   = 8./(12.*Grid.dx[v]);
    w2[v]   =-1./(12.*Grid.dx[v]);
  }

  double clk, file_time, dspaces_time, dspaces_wait_time;
	file_time = 0.0;
	// dspaces_put_req_t* dspaces_put_req_list;
	if(Run.use_dspaces_io) {
		dspaces_time = 0.0;
    dspaces_wait_time = 0.0;
    // dspaces_put_req_list = NULL;
	}

  char filename[128];

  MPI_File fhandle_mpi;
  int offset;

  static const cState* U = Grid.U;

  const double b_unit = sqrt(8.0*asin(1.0));

  double ihsz=1.0/((double (Grid.gsize[1]))*(double (Grid.gsize[2])));

  double dn,eps,vx,vy,vz,bx,by,bz,vsqr,idn,rfac;

  static int nvar=50;

  if (ini_flag) {
    int array_of_sizes[1];
    int array_of_subsizes[1];
    int array_of_starts[1];

    array_of_sizes[0]=Grid.gsize[0];
    array_of_subsizes[0]=Grid.lsize[0];
    array_of_starts[0]=Grid.beg[0]-Grid.gbeg[0];

    MPI_Type_create_subarray(1,array_of_sizes,array_of_subsizes,
			     array_of_starts,MPI_ORDER_FORTRAN,
			     MPI_FLOAT,&x_subarray);
    MPI_Type_commit(&x_subarray);

    int gsize[1];
    gsize[0] = Grid.gsize[0];
    io_file_log->analyze_vp = (struct log_entry*) malloc(sizeof(struct log_entry));
    log_entry_init(io_file_log->analyze_vp, "ANALYZE_VP", est_total_slice_iters, 1, gsize, nvar);
		if(Run.use_dspaces_io) {
      io_dspaces_log->analyze_vp = (struct log_entry*) malloc(sizeof(struct log_entry));
      log_entry_init(io_dspaces_log->analyze_vp, "ANALYZE_VP", est_total_slice_iters, 1, gsize, nvar);
      
      analyzevp_nvar = nvar;
      // dspaces_iput() is only called in ranks whose yz_rank == iroot
			// so the dspaces_put_req_list is only malloced there
      if(yz_rank == 0) {
        analyzevp_dspaces_put_req_list = (dspaces_put_req_t*) malloc(nvar*sizeof(dspaces_put_req_t));
      }
      analyzevp_buf = (float*) malloc(nvar*Grid.lsize[0]*sizeof(float));
    }

    ini_flag = 0;
  }

  double  *loc, *glo;
  float *iobuf;
  
  int bufsz = nvar*Grid.lsize[0];
 
  loc   = new double [bufsz];
  glo   = new double [bufsz];
  iobuf = new float[bufsz];

  for(j=0;j<bufsz;j++){ 
    loc[j] = 0.0;
  }
  
  ioff = Grid.lsize[0];

  LOCAL_LOOP(Grid,i,j,k) {
    node = Grid.node(i,j,k);

    ind = i-Grid.ghosts[0];

    dn   = U[node].d;
    idn  = 1.0/dn;

    vx   = U[node].M.x;
    vy   = U[node].M.y;   
    vz   = U[node].M.z;   
    bx   = U[node].B.x;
    by   = U[node].B.y;   
    bz   = U[node].B.z;

    vsqr = vx*vx+vy*vy+vz*vz;

    eps  = U[node].e*idn;

    // Output only stuff in corona
    if(Physics.rt_ext[i_ext_cor] == 1)
      rfac = Grid.Qthin[node] != 0.0;
    else if (Physics.rt_ext[i_ext_cor] == 2)
      rfac = (Grid.temp[node] > Physics.rt[i_rt_tr_tem]) && (Grid.pres[node] < Physics.rt[i_rt_tr_pre]);
    else
      rfac = 0.0;

    // mean state
    loc[ind+0*ioff]  += dn;
    loc[ind+1*ioff]  += eps;
    loc[ind+2*ioff]  += Grid.pres[node];
    loc[ind+3*ioff]  += Grid.temp[node];
    loc[ind+4*ioff]  += Grid.Tau[node];

    // v, B
    loc[ind+5*ioff]  += vx*vx;
    loc[ind+6*ioff]  += vy*vy+vz*vz;
    loc[ind+7*ioff]  += bx*bx;  
    loc[ind+8*ioff]  += (by*by+bz*bz);
    loc[ind+9*ioff]  += bx;
    loc[ind+10*ioff] += fabs(bx);
    loc[ind+11*ioff] += sqrt(by*by+bz*bz);
    loc[ind+12*ioff] += Grid.divB[node]*Grid.divB[node];
    
    // vertical MHD + conductive fluxes
    loc[ind+13*ioff] += vx*dn;
    loc[ind+14*ioff] += fabs(vx*dn);   
    loc[ind+15*ioff] += vx*dn*0.5*vsqr;
    loc[ind+16*ioff] += vx*(eps*dn+Grid.pres[node]);
    loc[ind+17*ioff] += eps+Grid.pres[node]/dn;
    loc[ind+18*ioff] += vx*(by*by+bz*bz)-bx*(vy*by+vz*bz);
    loc[ind+19*ioff] -= bx*(vy*by+vz*bz);

    if(Physics.params[i_param_spitzer] > 0.0)
      loc[ind+20*ioff] += Grid.sflx[node]*bx;

    double tvar1 = 0.0;
    double tvar2 = 0.0;
    double tvar3 = 0.0;
    double tvar4 = 1.0;
    double tvar5 = 0.0;
    double tvar6 = 0.0;
    double tvar7 = 0.0;
    double tvar8 = 0.0;
    double Qres  = 0.0;
    double Qvis  = 0.0;
    double Qamb  = 0.0;
    
    if (Run.diagnostics){
      tvar1 = Grid.tvar1[node];
      tvar2 = Grid.tvar2[node];
      tvar3 = Grid.tvar3[node];
      tvar4 = Grid.tvar4[node];
      tvar5 = Grid.tvar5[node];
      tvar6 = Grid.tvar6[node];
      tvar7 = Grid.tvar7[node];
      tvar8 = Grid.tvar8[node];
      Qres  = Grid.Qres[node];
      Qvis  = Grid.Qvis[node];
      if(Physics.params[i_param_ambipolar] > 0.0)
	Qamb  = Grid.Qamb[node];
    }

    double Qthin = 0.0;
    double QChr  = 0.0;
    double QCa   = 0.0;
    double QMg   = 0.0;
    double QH    = 0.0;

    if(Physics.rt_ext[i_ext_cor] >= 1)
      Qthin = Grid.Qthin[node];
    if(Physics.rt_ext[i_ext_cor] == 2){
      QChr  = Grid.QChr[node];
      QCa   = Grid.QCa[node];
      QMg   = Grid.QMg[node];
      QH    = Grid.QH[node];
    }

    // Radiation quantities
    loc[ind+21*ioff] += Grid.Qtot[node]*tvar4; // tvar4 -> efac
    loc[ind+22*ioff] += Qthin*tvar4;// tvar4 -> efac
    loc[ind+23*ioff] += QChr*tvar4;// tvar4 -> efac
    loc[ind+24*ioff] += QCa*tvar4;// tvar4 -> efac
    loc[ind+25*ioff] += QMg*tvar4;// tvar4 -> efac
    loc[ind+26*ioff] += QH*tvar4;// tvar4 -> efac

    // terms in plasma energy equation
    loc[ind+27*ioff] += dn*(eps+0.5*vsqr);
    loc[ind+28*ioff] += tvar1; // div(F_adv)
    loc[ind+29*ioff] += tvar2; // div(F_cond)
    loc[ind+30*ioff] += tvar3; // wlrt
    loc[ind+31*ioff] += tvar5; // m*g + damping
    loc[ind+32*ioff] += tvar6; // Boris
    loc[ind+33*ioff] += tvar7; // Tcheck
    loc[ind+34*ioff] += tvar8; // numerical diff
    loc[ind+35*ioff] += Qres;
    loc[ind+36*ioff] += Qvis;
    loc[ind+37*ioff] += Qamb;
    
    // terms in plasma energy equation in Corona only
    loc[ind+38*ioff] += rfac; // mask defining corona volume
    loc[ind+39*ioff] += rfac*dn*(eps+0.5*vsqr);
    loc[ind+40*ioff] += rfac*Grid.Qtot[node]*tvar4; // tvar4 -> efac
    loc[ind+41*ioff] += Qthin*tvar4;
    loc[ind+42*ioff] += rfac*tvar1; // div(F_adv)
    loc[ind+43*ioff] += rfac*tvar3; // wlrt
    loc[ind+44*ioff] += rfac*tvar5; // m*g + damping
    loc[ind+45*ioff] += rfac*tvar6; // Boris
    loc[ind+46*ioff] += rfac*tvar7; // Tcheck
    loc[ind+47*ioff] += rfac*tvar8; // numerical diff
    loc[ind+48*ioff] += rfac*Qres;
    loc[ind+49*ioff] += rfac*Qvis;  
  }

  MPI_Reduce(loc,glo,bufsz,MPI_DOUBLE,MPI_SUM,0,YZ_COMM);
  
  if(yz_rank==0){ // MPI_Reduce results are only meaningful on rank 0!

    // dspaces_put_req_list = (dspaces_put_req_t*) malloc(nvar*sizeof(dspaces_put_req_t));
    
    for(ind=0;ind<bufsz;ind++){
      glo[ind] *=  ihsz;
    }

    // B
    for(ind=0;ind<Grid.lsize[0];ind++){
      glo[ind+ 7*ioff] *= b_unit*b_unit;
      glo[ind+ 8*ioff] *= b_unit*b_unit;
      glo[ind+ 9*ioff] *= b_unit;
      glo[ind+10*ioff] *= b_unit;
      glo[ind+11*ioff] *= b_unit;
      glo[ind+12*ioff] *= b_unit*b_unit;
    }

    if(Run.use_dspaces_io && analyzevp_ref_count > 0) {
      clk = MPI_Wtime();
      for(int i=0; i<nvar; i++) {
        dspaces_check_put(ds_client, analyzevp_dspaces_put_req_list[i], 1);
      }
      double dspaces_check_time = MPI_Wtime() - clk;
      if(dspaces_check_time > nvar*dspaces_check_overhead) {
        dspaces_wait_time += MPI_Wtime() - clk - nvar*dspaces_check_overhead;
      }
    }
    if(Run.use_dspaces_io) {
      for(ind=0;ind<bufsz;ind++){ 
        iobuf[ind] = (float) glo[ind];
        analyzevp_buf[ind] = (float) glo[ind];
      }
    } else {
      for(ind=0;ind<bufsz;ind++){ 
        iobuf[ind] = (float) glo[ind];
      }
    }

    sprintf(filename,"%s%s.%06d",Run.path_2D,"hmean1D",Run.globiter);
    clk = MPI_Wtime();
    MPI_File_open(XCOL_COMM,filename,MPI_MODE_CREATE | MPI_MODE_WRONLY,
		  MPI_INFO_NULL,&fhandle_mpi);

    if( xcol_rank == 0 ){
      float header[4];            
      header[0] = (float) nvar;
      header[1] = (float) Grid.gsize[0];
      header[2] = (float) 1.0;      
      header[3] = (float) Run.time;
      MPI_File_write(fhandle_mpi,header,4,MPI_FLOAT,MPI_STATUS_IGNORE);
    }
    
    offset = 4*sizeof(float);
      
    MPI_File_set_view(fhandle_mpi,offset,MPI_FLOAT,x_subarray,(char*) "native",
		      MPI_INFO_NULL);
    MPI_File_write_all(fhandle_mpi,iobuf,nvar*Grid.lsize[0],MPI_FLOAT,
		       MPI_STATUS_IGNORE);
      
    MPI_File_close(&fhandle_mpi);

    file_time += MPI_Wtime() - clk;

    if(Run.use_dspaces_io) {
      // write 50 vars separately
      char ds_var_name[128];
      uint64_t lb[1], ub[1];
      lb[0] = Grid.beg[0]-Grid.gbeg[0];
      ub[0] = lb[0] + Grid.lsize[0] - 1;
      for(int v=0; v<nvar; v++) {
        sprintf(ds_var_name, "%s%s_%d", Run.path_2D, "hmean1D", v);
        clk = MPI_Wtime();
        analyzevp_dspaces_put_req_list[v] = dspaces_iput(ds_client, ds_var_name, Run.globiter,
                                                         sizeof(float), 1, lb, ub,
                                                         &analyzevp_buf[v*Grid.lsize[0]], 0, 0);
        dspaces_time += MPI_Wtime() - clk;
      }
      char header_name[128];
      FILE * hfhandle = NULL;
      if( xcol_rank == 0 ){
        sprintf(header_name, "%s%s.header", Run.path_2D, "hmean1D");
        hfhandle=fopen(header_name, "w");
        float header[4];            
        header[0] = (float) nvar;
        header[1] = (float) Grid.gsize[0];
        header[2] = (float) 1.0;      
        header[3] = (float) Run.time;
        fwrite(header, sizeof(float), 4, hfhandle);
        fclose(hfhandle);
      }
    }

    

  }
  
  delete[] loc;
  delete[] glo;

  // if(yz_rank == 0 && Run.use_dspaces_io) {
  //   double dspaces_overlap_time = MPI_Wtime() - clk;
  //   clk = MPI_Wtime();
  //   for(int i=0; i<nvar; i++) {
  //     dspaces_check_put(ds_client, dspaces_put_req_list[i], 1);
  //   }
  //   dspaces_wait_time += MPI_Wtime() - clk;
  //   free(dspaces_put_req_list);
  // }
  delete[] iobuf;

  if(yz_rank == 0 && xcol_rank == 0) {
    io_file_log->analyze_vp->iter[io_file_log->analyze_vp->count] = Run.globiter;
		io_file_log->analyze_vp->api_time[io_file_log->analyze_vp->count] = file_time;
		io_file_log->analyze_vp->time[io_file_log->analyze_vp->count] = file_time;
    io_file_log->analyze_vp->count++ ;
		if(Run.use_dspaces_io) {
			io_dspaces_log->analyze_vp->iter[io_dspaces_log->analyze_vp->count] = Run.globiter;
      io_dspaces_log->analyze_vp->api_time[io_dspaces_log->analyze_vp->count] = dspaces_time;
      if(io_dspaces_log->analyze_vp->count > 0) {
        io_dspaces_log->analyze_vp->wait_time[io_dspaces_log->analyze_vp->count-1] = dspaces_wait_time;
        io_dspaces_log->analyze_vp->time[io_dspaces_log->analyze_vp->count-1] = dspaces_wait_time
                                  + io_dspaces_log->analyze_vp->api_time[io_dspaces_log->analyze_vp->count-1];
      }
      io_dspaces_log->analyze_vp->count++ ;
		}
    if(Run.verbose > 0) {
      std::cout << "File Output (ANALYZE_VP) in " << file_time << " seconds" << std::endl;
      if(Run.use_dspaces_io) {
        std::cout << "DataSpaces API Call (ANALYZE_VP) in " << dspaces_time
                  << " seconds" << std::endl;
        // std::cout << "DataSpaces Wait (ANALYZE_VP) in " << dspaces_wait_time
        //           << " seconds" << std::endl;
        // std::cout << "DataSpaces Output (ANALYZE_VP) in " << dspaces_time
        //           << " seconds" << std::endl;
      }
    }
  }
  analyzevp_ref_count++;
}
