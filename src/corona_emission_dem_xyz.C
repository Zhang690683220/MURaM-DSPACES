#include <mpi.h>
#include "physics.H"
#include "grid.H"
#include "run.H"
#include <math.h>
#include <stdlib.h>
#include <fstream>
#include <stdio.h>
#include "io.H"
#include "comm_split.H"

using namespace std;

#define OUTER_LOOP(G,i,j,d1,d2) \
  for((j)=(G)[(d2)][0];(j)<=(G)[(d2)][1];(j)++) \
  for((i)=(G)[(d1)][0];(i)<=(G)[(d1)][1];(i)++)

//extern void slice_write(const GridData&,const int,float*,int,int,const int,
//			const int,FILE*);

extern est_total_slice_iters;
extern struct log *io_file_log, *io_dspaces_log;

extern double slice_write_rebin(const GridData&,const int,float*,const int,const int,
                              const int,const int,const int,const int,FILE*);

extern dspaces_put_req_t* slice_write_rebin_dspaces(const GridData& Grid,
           											const int iroot, float* vloc,
		       										const int nloc,const int nvar,const int n0,
		       										const int n1,const int sm_x,const int sm_y,
           											char* filename, const int iter, const int ndim);

float *coronaxy_buf = NULL;
float *coronaxz_buf = NULL;
float *coronayz_buf = NULL;
int corona_nout;
int corona_nslvar;
dspaces_put_req_t** coronaxy_dspaces_put_req_list = NULL;
dspaces_put_req_t** coronaxz_dspaces_put_req_list = NULL;
dspaces_put_req_t** coronayz_dspaces_put_req_list = NULL;

inline int imin(int a, int b) { return a < b ? a : b; }
inline int imax(int a, int b) { return a > b ? a : b; }
//======================================================================
void corona_emission_dem_xyz(const RunData&  Run, const GridData& Grid, 
		               const PhysicsData& Physics) {

  static int ini_flag = 1;
	static int corona_ref_count = 0;
  
  const int iroot = 0;
  
  const double lgTmin = 4.5;
  const double dellgT = 0.1;
  
  const int nout = 4;
  
  const int rebin[3] = {3,1,1};
  
  register int i, j, k, node1, node2, ind, v, v1, v2, ind1, d, d1, d2, d3;
  
  int str, offset;
  
  int stride[3];
  
  stride[0] = Grid.stride[0];
  stride[1] = Grid.stride[1];
  stride[2] = Grid.stride[2];
  
  int bounds[3][2];
  
  for(d=0;d<3;d++){
    bounds[d][0] = Grid.lbeg[d];
    bounds[d][1] = Grid.lend[d];
  }
  
  const int loop_order[3][3] = {{ 1, 2, 0 },{ 0, 1, 2 },{ 2, 1, 0 }};
  
  double* los_sum_loc;
  double* los_sum;
  float*  io_buf;

  char filename[128];

  double r14a,r14b,r14,t6a,t6b,va,vb,vv,rfac,dx,p_a,p_b,s,t1,t2,tmin,tmax,lgTmax;

  double lgT0,del0;

  tmax = 0;
  LOCAL_LOOP(Grid,i,j,k){
    tmax = max(tmax,Grid.temp[Grid.node(i,j,k)]);
  }
  
  MPI_Allreduce(&tmax,&lgTmax,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  
  lgTmax = log10(lgTmax);

  if (lgTmax < lgTmin)
    lgTmax = lgTmin;

  lgT0=lgTmin;
  del0=1.0/dellgT;
  
  int nslvar = (int) ((lgTmax-lgTmin)/dellgT + 1);
   
  double* tlev = new double[nslvar+1];
  
  for(i=0;i<nslvar+1;i++)
    tlev[i] = lgTmin+double(i)*dellgT;

  if(ini_flag){
    if(Run.rank==0){
      cout << "DEM: Use log(T), log(rho)" << endl;
    }

		int gsize[2]={0};
		io_file_log->corona = (struct log_entry*) malloc(sizeof(struct log_entry));
    log_entry_init(io_file_log->corona, "CORONA", est_total_slice_iters, 2, gsize, nout*nslvar);
		io_file_log->corona->corona_gsize[0][0] = Grid.gsize[2];
		io_file_log->corona->corona_gsize[0][1] = Grid.gsize[0];
		io_file_log->corona->corona_gsize[1][0] = Grid.gsize[1];
		io_file_log->corona->corona_gsize[1][1] = Grid.gsize[2];
		io_file_log->corona->corona_gsize[2][0] = Grid.gsize[1];
		io_file_log->corona->corona_gsize[2][1] = Grid.gsize[0];
		if(Run.use_dspaces_io) {
      io_dspaces_log->corona = (struct log_entry*) malloc(sizeof(struct log_entry));
      log_entry_init(io_dspaces_log->corona, "CORONA", est_total_slice_iters,2, gsize, nout*nslvar);
			io_dspaces_log->corona->corona_gsize[0][0] = Grid.gsize[2];
			io_dspaces_log->corona->corona_gsize[0][1] = Grid.gsize[0];
			io_dspaces_log->corona->corona_gsize[1][0] = Grid.gsize[1];
			io_dspaces_log->corona->corona_gsize[1][1] = Grid.gsize[2];
			io_dspaces_log->corona->corona_gsize[2][0] = Grid.gsize[1];
			io_dspaces_log->corona->corona_gsize[2][1] = Grid.gsize[0];

			corona_nout = nout;
			corona_nslvar = nslvar;
			// dspaces_iput() is only called in ranks whose xcol_rank == iroot, ycol_rank == iroot,
			// zcol_rank == iroot, so the according dspaces_put_req_list is only malloced there
			if(zcol_rank == iroot) {
				coronaxy_dspaces_put_req_list = (dspaces_put_req_t**) malloc(nout*sizeof(dspaces_put_req_t*));
			}
			if(ycol_rank == iroot) {
				coronaxz_dspaces_put_req_list = (dspaces_put_req_t**) malloc(nout*sizeof(dspaces_put_req_t*));
			}
			if(xcol_rank == iroot) {
				coronayz_dspaces_put_req_list = (dspaces_put_req_t**) malloc(nout*sizeof(dspaces_put_req_t*));
			}
    	coronaxy_buf = (float*) malloc(nout*nslvar*Grid.lsize[1]*Grid.lsize[0]*sizeof(float));
			coronaxz_buf = (float*) malloc(nout*nslvar*Grid.lsize[2]*Grid.lsize[0]*sizeof(float));
			coronayz_buf = (float*) malloc(nout*nslvar*Grid.lsize[1]*Grid.lsize[2]*sizeof(float));
    }

    ini_flag = 0;
  }

  if(Run.rank==0)
    cout << "DEM: lgTmin = " << lgTmin <<  " lgTmax = " << lgTmax << ' ' << nslvar  << endl;

  FILE* fhandle=NULL;

	double clk, file_time, dspaces_time, dspaces_wait_time;
	int put_count;
	enum rank_group {X_COL, Y_COL, Z_COL};
	enum rank_group rank_history;
	file_time = 0.0;
	// dspaces_put_req_t* dspaces_put_req_list;
	if(Run.use_dspaces_io) {
	  dspaces_time = 0.0;
	  dspaces_wait_time = 0.0;
      // dspaces_put_req_list = NULL;
	  put_count = 0;
	}

  for(d=0;d<Grid.NDIM;d++){
    d1=loop_order[d][0];
    d2=loop_order[d][1];
    d3=loop_order[d][2];

    int localsize  = Grid.lsize[d2]*Grid.lsize[d3];
 
    los_sum_loc = new double[nout*nslvar*localsize];
    los_sum     = new double[nout*nslvar*localsize];
    io_buf      = new float[nout*nslvar*localsize];
     
    for(v=0;v<nout*nslvar*localsize;v++){
      los_sum_loc[v] = 0.0;
      los_sum[v]     = 0.0;
      io_buf[v]      = 0.0;
    }
    
    str = stride[d1];
    dx  = Grid.dx[d1];
    OUTER_LOOP(bounds,j,k,d2,d3){
      offset = j*stride[d2]+k*stride[d3];
      ind = j-bounds[d2][0] +(k-bounds[d3][0])*Grid.lsize[d2]; 
      for(i=bounds[d1][0]-1;i<=bounds[d1][1];i++){
        node1 = offset+i*str;
	node2 = offset+(i+1)*str;

	t6a  = log10(Grid.temp[node1]);
	t6b  = log10(Grid.temp[node2]);
	r14a = log(Grid.U[node1].d*1e14);
	r14b = log(Grid.U[node2].d*1e14);
	va   = Grid.U[node1].M[d1];
	vb   = Grid.U[node2].M[d1];
	
	ind1 = ind;

	tmin = min(t6a,t6b);
	tmax = max(t6a,t6b);

	v1 = (int) ( (tmin-lgT0)*del0 );
	v2 = (int) ( (tmax-lgT0)*del0 );

	v1 = imax(0,v1);
	v2 = imin(nslvar-1,v2);

	ind1 += v1*localsize;

	for(v=v1;v<=v2;v++){
  
	  t1=max(tmin,tlev[v]);
	  t2=min(tmax,tlev[v+1]);
		 
	  if(t2 > t1){
	    rfac = (t2-t1)/(tmax-tmin);
	    s    = 0.5*(t1+t2);
	    p_a  = (t6b-s)/(t6b-t6a);
	    p_b  = (s-t6a)/(t6b-t6a);
	  } else if (t2 == t1) {
	    rfac = 1.0;
	    p_a  = 0.5;
	    p_b  = 0.5;
	  } else {
	    rfac = 0.0;
	    p_a  = 0.5;
	    p_b  = 0.5;
	  }

	  r14  = exp(p_a*r14a+p_b*r14b);
	  r14  = r14*r14;
	  vv   = p_a*va+p_b*vb;

	  if(r14 > 1e10) rfac = 0.0;
	    
	  los_sum_loc[ind1]                    += rfac*dx;
	  los_sum_loc[ind1+1*nslvar*localsize] += rfac*dx*r14;
	  los_sum_loc[ind1+2*nslvar*localsize] += rfac*dx*r14*vv;
	  los_sum_loc[ind1+3*nslvar*localsize] += rfac*dx*r14*vv*vv;
	  
	  ind1 += localsize;
	}
	
      }
    }

	// // check dspaces_iput() except for the first iter
	// if(Run.use_dspaces_io && d > 0) {
	//   double dspaces_overlap_time = MPI_Wtime() - clk;
    //   clk = MPI_Wtime();
	//   for(int i=0; i<nslvar; i++) {
	// 	dspaces_check_put(ds_client, dspaces_put_req_list[i], 1);
	//   }
	//   double dspaces_wait_time = MPI_Wtime() - clk;
	//   if(dspaces_wait_time > 1e-6) {
	// 	dspaces_wait_time += dspaces_wait_time + dspaces_overlap_time;
	//   }
	//   free(dspaces_put_req_list);
	// }

	
    if(d1 == 0){
      MPI_Reduce(los_sum_loc,los_sum,nout*nslvar*localsize,MPI_DOUBLE,MPI_SUM,iroot,
		 XCOL_COMM);

      for(i=0;i<nslvar*localsize;i++){
	rfac = los_sum[i+nslvar*localsize]; 
	if( rfac > 0.0 ){
	  los_sum[i+2*nslvar*localsize] /= rfac; 
	  los_sum[i+3*nslvar*localsize] /= rfac;
	}
      }

      for(i=0;i<nslvar*localsize;i++){
	rfac = los_sum[i+3*nslvar*localsize]-los_sum[i+2*nslvar*localsize]*los_sum[i+2*nslvar*localsize];
	los_sum[i+3*nslvar*localsize] = sqrt(max(rfac,0.0));
      }

	  // check dspaces_iput() except for the first iter
	  // if(Run.use_dspaces_io && dspaces_put_req_list != NULL) {
		// 	double dspaces_overlap_time = MPI_Wtime() - clk;
    // 	clk = MPI_Wtime();
	  // 	for(int i=0; i<nslvar; i++) {
		//   	dspaces_check_put(ds_client, dspaces_put_req_list[i], 1);
	  // 	}
	  // 	double dspaces_wait_time = MPI_Wtime() - clk;
	  // 	if(dspaces_wait_time > 1e-6) {
		//   	dspaces_wait_time += dspaces_wait_time + dspaces_overlap_time;
	  // 	}
	  // 	free(dspaces_put_req_list);
	  // }

		// if(corona_ref_count > 0) {
		// 	for(int j=0; j<nout; j++) {
		// 		for(int i=0; i<nslvar; i++) {
		// 			dspaces_check_put(ds_client, coronayz_dspaces_put_req_list[j][i], 1);
		// 		}
		// 	}
		// }
	
	  // update io_buf values
      // for(v=0;v<nout*nslvar*localsize;v++) {
			// 	io_buf[v] = (float) los_sum[v];
			// 	coronayz_buf[v] = (float) los_sum[v];
			// }
      if (xcol_rank == iroot){
	for (v=0;v<nout;v++){

		if(Run.use_dspaces_io && corona_ref_count > 0) {
			clk = MPI_Wtime();
			for(int i=0; i<nslvar; i++) {
				// slice_write_rebin() sometimes has all involved ranks write
				// sometimes gather the data to the root rank
				if(coronayz_dspaces_put_req_list[v][i] != NULL) {
					dspaces_check_put(ds_client, coronayz_dspaces_put_req_list[v][i], 1);
				}
			}
			double dspaces_check_time = MPI_Wtime() - clk;
			if(dspaces_check_time > nslvar*1e-6) {
				dspaces_wait_time += MPI_Wtime() - clk;
			}
			free(coronayz_dspaces_put_req_list[v]);
		}

		if(Run.use_dspaces_io) {
			// move iobuf update inside the nout loop
			for(int i=0; i<nslvar*localsize; i++) {
				io_buf[i+v*nslvar*localsize] = (float) los_sum[i+v*nslvar*localsize];
				coronayz_buf[i+v*nslvar*localsize] = (float) los_sum[i+v*nslvar*localsize];
			}
		} else {
			for(int i=0; i<nslvar*localsize; i++) {
				io_buf[i+v*nslvar*localsize] = (float) los_sum[i+v*nslvar*localsize];
			}
		}
	  if(yz_rank == 0) {
	    if(v == 0)
	      sprintf(filename,"%s%s.%06d",Run.path_2D,"corona_emission_adj_fil_x",Run.globiter);
	    if(v == 1)
	      sprintf(filename,"%s%s.%06d",Run.path_2D,"corona_emission_adj_dem_x",Run.globiter);
	    if(v == 2)
	      sprintf(filename,"%s%s.%06d",Run.path_2D,"corona_emission_adj_vlos_x",Run.globiter);
	    if(v == 3)
	      sprintf(filename,"%s%s.%06d",Run.path_2D,"corona_emission_adj_vrms_x",Run.globiter);
	    
	    fhandle=fopen(filename,"w");

	    float header[6];            
	    header[0] = (float) nslvar;
	    header[1] = (float) Grid.gsize[d2]/rebin[d2];
	    header[2] = (float) Grid.gsize[d3]/rebin[d3];
	    header[3] = (float) Run.time;
	    header[4] = (float) lgTmin;
	    header[5] = (float) dellgT;
	    fwrite(header,sizeof(float),6,fhandle);
	  }
	
	  //slice_write(Grid,0,&(io_buf[v*nslvar*localsize]),localsize,nslvar,d2,d3,fhandle);
		clk = MPI_Wtime();
	  slice_write_rebin(Grid,0,&(io_buf[v*nslvar*localsize]),localsize,nslvar,d2,d3,rebin[d2],rebin[d3],fhandle);
	  file_time += MPI_Wtime() - clk;

	  if(yz_rank == 0)
	    fclose(fhandle);

		if(Run.use_dspaces_io) {
			// var name passed to dspaces
			if(v == 0)
	      sprintf(filename,"%s%s",Run.path_2D,"corona_emission_adj_fil_x");
	    if(v == 1)
	      sprintf(filename,"%s%s",Run.path_2D,"corona_emission_adj_dem_x");
	    if(v == 2)
	      sprintf(filename,"%s%s",Run.path_2D,"corona_emission_adj_vlos_x");
	    if(v == 3)
	      sprintf(filename,"%s%s",Run.path_2D,"corona_emission_adj_vrms_x");

			// header file to be written by rank 0
			char header_filename[128];
			FILE * hfhandle = NULL;
			if(yz_rank == 0) {
	      sprintf(header_filename,"%s%s",filename,".header");

				hfhandle=fopen(header_filename, "w");
				float header[6];            
	    	header[0] = (float) nslvar;
	    	header[1] = (float) Grid.gsize[d2]/rebin[d2];
	    	header[2] = (float) Grid.gsize[d3]/rebin[d3];
	    	header[3] = (float) Run.time;
	    	header[4] = (float) lgTmin;
	    	header[5] = (float) dellgT;
	    	fwrite(header,sizeof(float),6,hfhandle);
				fclose(hfhandle);
			}
			clk = MPI_Wtime();
			coronayz_dspaces_put_req_list[v] = slice_write_rebin_dspaces(Grid, 0,
																																	 &coronayz_buf[v*nslvar*localsize],
															 																		 localsize, nslvar, d2, d3, rebin[d2],
															 																		 rebin[d3], filename, Run.globiter, 2);
			dspaces_time += MPI_Wtime() - clk;
		}
	}
	if(Run.use_dspaces_io) {
	  put_count++;
	  rank_history = X_COL;
	  clk = MPI_Wtime();
	}
      }
    }
    
    if(d1 == 1){    
      MPI_Reduce(los_sum_loc,los_sum,nout*nslvar*localsize,MPI_DOUBLE,MPI_SUM,iroot,
		 YCOL_COMM);
      
      for(i=0;i<nslvar*localsize;i++){
	rfac = los_sum[i+nslvar*localsize]; 
	if( rfac > 0.0 ){
	  los_sum[i+2*nslvar*localsize] /= rfac; 
	  los_sum[i+3*nslvar*localsize] /= rfac;
	}
      }

      for(i=0;i<nslvar*localsize;i++){
	rfac = los_sum[i+3*nslvar*localsize]-los_sum[i+2*nslvar*localsize]*los_sum[i+2*nslvar*localsize];
	los_sum[i+3*nslvar*localsize] = sqrt(max(rfac,0.0));
      }

	  // check dspaces_iput() except for the first iter
	  // if(Run.use_dspaces_io && dspaces_put_req_list != NULL) {
		// 	double dspaces_overlap_time = MPI_Wtime() - clk;
    // 	clk = MPI_Wtime();
	  // 	for(int i=0; i<nslvar; i++) {
		//   	dspaces_check_put(ds_client, dspaces_put_req_list[i], 1);
	  // 	}
	  // 	double dspaces_wait_time = MPI_Wtime() - clk;
	  // 	if(dspaces_wait_time > 1e-6) {
		//   	dspaces_wait_time += dspaces_wait_time + dspaces_overlap_time;
	  // 	}
	  // 	free(dspaces_put_req_list);
	  // }
	  
	  // update io_buf values
      // for(v=0;v<nout*nslvar*localsize;v++) {
			// 	io_buf[v] = (float) los_sum[v];
			// 	coronaxz_buf = (float) los_sum[v];
			// }
      
      if (ycol_rank == iroot){
	for (v=0;v<nout;v++){
		if(Run.use_dspaces_io && corona_ref_count > 0) {
			clk = MPI_Wtime();
			for(int i=0; i<nslvar; i++) {
				// slice_write_rebin() sometimes has all involved ranks write
				// sometimes gather the data to the root rank
				if(coronaxz_dspaces_put_req_list[v][i] != NULL) {
					dspaces_check_put(ds_client, coronaxz_dspaces_put_req_list[v][i], 1);
				}
			}
			double dspaces_check_time = MPI_Wtime() - clk;
			if(dspaces_check_time > nslvar*1e-6) {
				dspaces_wait_time += MPI_Wtime() - clk;
			}
			free(coronaxz_dspaces_put_req_list[v]);
		}
		if(Run.use_dspaces_io) {
			// move iobuf update inside nout loop
			for(int i=0; i<nslvar*localsize; i++) {
				io_buf[i+v*nslvar*localsize] = (float) los_sum[i+v*nslvar*localsize];
				coronaxz_buf[i+v*nslvar*localsize] = (float) los_sum[i+v*nslvar*localsize];
			}
		} else {
			for(int i=0; i<nslvar*localsize; i++) {
				io_buf[i+v*nslvar*localsize] = (float) los_sum[i+v*nslvar*localsize];
			}
		}
	  if(xz_rank == 0) {
	    if(v == 0)
	      sprintf(filename,"%s%s.%06d",Run.path_2D,"corona_emission_adj_fil_y",Run.globiter);
	    if(v == 1)
	      sprintf(filename,"%s%s.%06d",Run.path_2D,"corona_emission_adj_dem_y",Run.globiter);
	    if(v == 2)
	      sprintf(filename,"%s%s.%06d",Run.path_2D,"corona_emission_adj_vlos_y",Run.globiter);
	    if(v == 3)
	      sprintf(filename,"%s%s.%06d",Run.path_2D,"corona_emission_adj_vrms_y",Run.globiter);
	    
	    fhandle=fopen(filename,"w");
	    
	    float header[6];            
	    header[0] = (float) nslvar;
	    header[1] = (float) Grid.gsize[d2]/rebin[d2];
	    header[2] = (float) Grid.gsize[d3]/rebin[d3];
	    header[3] = (float) Run.time;
	    header[4] = (float) lgTmin;
	    header[5] = (float) dellgT;
	    fwrite(header,sizeof(float),6,fhandle);
	  }
	
	  //slice_write(Grid,0,&(io_buf[v*nslvar*localsize]),localsize,nslvar,d2,d3,fhandle);
		clk = MPI_Wtime();
	  slice_write_rebin(Grid,0,&(io_buf[v*nslvar*localsize]),localsize,nslvar,d2,d3,rebin[d2],rebin[d3],fhandle);
	  file_time += MPI_Wtime() - clk;

	  if(xz_rank == 0)
	    fclose(fhandle);

		if(Run.use_dspaces_io) {
			// var name passed to dspaces
			if(v == 0)
	      sprintf(filename,"%s%s",Run.path_2D,"corona_emission_adj_fil_y");
	    if(v == 1)
	      sprintf(filename,"%s%s",Run.path_2D,"corona_emission_adj_dem_y");
	    if(v == 2)
	      sprintf(filename,"%s%s",Run.path_2D,"corona_emission_adj_vlos_y");
	    if(v == 3)
	      sprintf(filename,"%s%s",Run.path_2D,"corona_emission_adj_vrms_y");

			// header file to be written by rank 0
			char header_filename[128];
			FILE * hfhandle = NULL;
			if(xz_rank == 0) {
				sprintf(header_filename,"%s%s",filename,".header");

				hfhandle=fopen(header_filename, "w");
				float header[6];            
	    	header[0] = (float) nslvar;
	    	header[1] = (float) Grid.gsize[d2]/rebin[d2];
	    	header[2] = (float) Grid.gsize[d3]/rebin[d3];
	    	header[3] = (float) Run.time;
	    	header[4] = (float) lgTmin;
	    	header[5] = (float) dellgT;
	    	fwrite(header,sizeof(float),6,hfhandle);
				fclose(hfhandle);
			}
			clk = MPI_Wtime();
			coronaxz_dspaces_put_req_list[v] = slice_write_rebin_dspaces(Grid, 0,
																																	 &coronaxz_buf[v*nslvar*localsize], 
															 																		 localsize, nslvar, d2, d3, rebin[d2],
															 																		 rebin[d3], filename, Run.globiter, 2);
			dspaces_time += MPI_Wtime() - clk;
		}
	}
	if(Run.use_dspaces_io) {
	  put_count++;
	  rank_history = Y_COL;
	  clk = MPI_Wtime();
	}
      }
    }
      
    if(d1 == 2){    
      MPI_Reduce(los_sum_loc,los_sum,nout*nslvar*localsize,MPI_DOUBLE,MPI_SUM,iroot,
		 ZCOL_COMM);
    
      for(i=0;i<nslvar*localsize;i++){
	rfac = los_sum[i+nslvar*localsize]; 
	if( rfac > 0.0 ){
	  los_sum[i+2*nslvar*localsize] /= rfac; 
	  los_sum[i+3*nslvar*localsize] /= rfac;
	}
      }

      for(i=0;i<nslvar*localsize;i++){
	rfac = los_sum[i+3*nslvar*localsize]-los_sum[i+2*nslvar*localsize]*los_sum[i+2*nslvar*localsize];
	los_sum[i+3*nslvar*localsize] = sqrt(max(rfac,0.0));
      }

	  // check dspaces_iput() except for the first iter
	  // if(Run.use_dspaces_io && dspaces_put_req_list != NULL) {
		// 	double dspaces_overlap_time = MPI_Wtime() - clk;
    // 	clk = MPI_Wtime();
	  // 	for(int i=0; i<nslvar; i++) {
		//   	dspaces_check_put(ds_client, dspaces_put_req_list[i], 1);
	  // 	}
	  // 	double dspaces_wait_time = MPI_Wtime() - clk;
	  // 	if(dspaces_wait_time > 1e-6) {
		//   	dspaces_wait_time += dspaces_wait_time + dspaces_overlap_time;
	  // 	}
	  // 	free(dspaces_put_req_list);	
	  // }
	  // update io_buf values
      for(v=0;v<nout*nslvar*localsize;v++) {
				io_buf[v] = (float) los_sum[v];
				coronaxy_buf[v] = (float) los_sum[v];
			}

      if (zcol_rank == iroot){
	for (v=0;v<nout;v++){
		if(Run.use_dspaces_io && corona_ref_count > 0) {
			clk = MPI_Wtime();
			for(int i=0; i<nslvar; i++) {
				// slice_write_rebin() sometimes has all involved ranks write
				// sometimes gather the data to the root rank
				if(coronaxy_dspaces_put_req_list[v][i]) {
					dspaces_check_put(ds_client, coronaxy_dspaces_put_req_list[v][i], 1);
				}
			}
			double dspaces_check_time = MPI_Wtime() - clk;
			if(dspaces_check_time > nslvar*1e-6) {
				dspaces_wait_time += MPI_Wtime() - clk;
			}
			free(coronaxy_dspaces_put_req_list[v]);
		}
		if(Run.use_dspaces_io) {
			// move iobuf update inside nout loop
			for(int i=0; i<nslvar*localsize; i++) {
				io_buf[i+v*nslvar*localsize] = (float) los_sum[i+v*nslvar*localsize];
				coronaxy_buf[i+v*nslvar*localsize] = (float) los_sum[i+v*nslvar*localsize];
			}
		} else {
			for(int i=0; i<nslvar*localsize; i++) {
				io_buf[i+v*nslvar*localsize] = (float) los_sum[i+v*nslvar*localsize];
			}
		}

	  if(xy_rank == 0) {
	    if(v == 0)
	      sprintf(filename,"%s%s.%06d",Run.path_2D,"corona_emission_adj_fil_z",Run.globiter);
	    if(v == 1)
	      sprintf(filename,"%s%s.%06d",Run.path_2D,"corona_emission_adj_dem_z",Run.globiter);
	    if(v == 2)
	      sprintf(filename,"%s%s.%06d",Run.path_2D,"corona_emission_adj_vlos_z",Run.globiter);
	    if(v == 3)
	      sprintf(filename,"%s%s.%06d",Run.path_2D,"corona_emission_adj_vrms_z",Run.globiter);
	
	    fhandle=fopen(filename,"w");

	    float header[6];            
	    header[0] = (float) nslvar;
	    header[1] = (float) Grid.gsize[d2]/rebin[d2];
	    header[2] = (float) Grid.gsize[d3]/rebin[d3];
	    header[3] = (float) Run.time;
	    header[4] = (float) lgTmin;
	    header[5] = (float) dellgT;
	    fwrite(header,sizeof(float),6,fhandle);
	  }
	
	  //slice_write(Grid,0,&(io_buf[v*nslvar*localsize]),localsize,nslvar,d2,d3,fhandle);
		clk = MPI_Wtime();
	  slice_write_rebin(Grid,0,&(io_buf[v*nslvar*localsize]),localsize,nslvar,d2,d3,rebin[d2],rebin[d3],fhandle);
	  file_time += MPI_Wtime() - clk;

	  if(xy_rank == 0)
	    fclose(fhandle);

		if(Run.use_dspaces_io) {
			// var name passed to dspaces
			if(v == 0)
	      sprintf(filename,"%s%s",Run.path_2D,"corona_emission_adj_fil_z");
	    if(v == 1)
	      sprintf(filename,"%s%s",Run.path_2D,"corona_emission_adj_dem_z");
	    if(v == 2)
	      sprintf(filename,"%s%s",Run.path_2D,"corona_emission_adj_vlos_z");
	    if(v == 3)
	      sprintf(filename,"%s%s",Run.path_2D,"corona_emission_adj_vrms_z");

			// header file to be written by rank 0
			char header_filename[128];
			FILE * hfhandle = NULL;
			if(xy_rank == 0) {
				sprintf(header_filename,"%s%s",filename,".header");

				hfhandle=fopen(header_filename, "w");
				float header[6];            
	    	header[0] = (float) nslvar;
	    	header[1] = (float) Grid.gsize[d2]/rebin[d2];
	    	header[2] = (float) Grid.gsize[d3]/rebin[d3];
	    	header[3] = (float) Run.time;
	    	header[4] = (float) lgTmin;
	    	header[5] = (float) dellgT;
	    	fwrite(header,sizeof(float),6,hfhandle);
				fclose(hfhandle);
			}
			clk = MPI_Wtime();
			coronaxy_dspaces_put_req_list[v] = slice_write_rebin_dspaces(Grid, 0,
																																	 &coronaxy_buf[v*nslvar*localsize],
															 																		 localsize, nslvar, d2, d3, rebin[d2],
															 																		 rebin[d3], filename, Run.globiter, 2);
			dspaces_time += MPI_Wtime() - clk;
		}
	}
	if(Run.use_dspaces_io) {
	  put_count++;
	  rank_history = Z_COL;
	  clk = MPI_Wtime();
	}
      }
	
    }
    
    delete[] los_sum_loc;
    delete[] los_sum;
    delete[] io_buf;
  }

  // if(Run.use_dspaces_io && dspaces_put_req_list != NULL) {
  //   for(int i=0; i<nslvar; i++) {
  //     dspaces_check_put(ds_client, dspaces_put_req_list[i], 1);
  //   }
  //   dspaces_wait_time += MPI_Wtime() - clk;
  //   free(dspaces_put_req_list);
  // }
  

  delete[] tlev;

	if(Run.rank == 0) {
		io_file_log->corona->iter[io_file_log->corona->count] = Run.globiter;
		io_file_log->corona->api_time[io_file_log->corona->count] = file_time;
		io_file_log->corona->time[io_file_log->corona->count] = file_time;
		io_file_log->corona->count++ ;
		if(Run.use_dspaces_io) {
			io_dspaces_log->corona->iter[io_dspaces_log->corona->count] = Run.globiter;
      io_dspaces_log->corona->api_time[io_dspaces_log->corona->count] = dspaces_time;
      if(io_dspaces_log->corona->count > 0) {
        io_dspaces_log->corona->wait_time[io_dspaces_log->corona->count-1] = dspaces_wait_time;
        io_dspaces_log->corona->time[io_dspaces_log->corona->count-1] = dspaces_wait_time
                                    + io_dspaces_log->corona->api_time[io_dspaces_log->corona->count-1];
      }
			io_dspaces_log->corona->count++ ;
		}
		if(Run.verbose > 0) {
	  	std::cout << "File Output (Corona_XYZ) in " << file_time << " seconds" << std::endl;
      if(Run.use_dspaces_io) {
				std::cout << "DataSpaces API Call (Corona_XYZ) in " << dspaces_time
									<< " seconds" << std::endl;
    		// std::cout << "DataSpaces Wait (Corona_XYZ) in " << dspaces_wait_time
				// 					<< " seconds" << std::endl;
    		// std::cout << "DataSpaces Output (Corona_XYZ) in " << dspaces_time
				// 					<< " seconds" << std::endl;
			}
		}
	}
	corona_ref_count++;
}

