#include <mpi.h>
#include <math.h>
#include <stdlib.h>
#include <fstream>
#include <stdio.h>
#include "physics.H"
#include "grid.H"
#include "run.H"
#include "comm_split.H"
#include "rt/rt.h"
#include "io.H"

using namespace std;

typedef double realtype;

extern est_total_slice_iters;
extern struct log *io_file_log, *io_dspaces_log;

extern void slice_write(const GridData&,const int,float*,int,int,const int,
			const int,FILE*);

extern dspaces_put_req_t* slice_write_dspaces(const GridData& Grid, const int iroot,
                                float* vloc, int nloc, int nvar,int n0,
                                int n1, char* filename, const int iter,
                                const int ndim);

//======================================================================
void tau_slice(const RunData&  Run, const GridData& Grid, 
	      const PhysicsData& Physics,RTS *rts) {

  static int ini_flag = 1;

  const int iroot = 0;

  register int i, j, k, ind, nsl, v, ind1, node1, node2;

  int ibeg = Grid.lbeg[0], iend = Grid.lend[0];
  int jbeg = Grid.lbeg[1], jend = Grid.lend[1];
  int kbeg = Grid.lbeg[2], kend = Grid.lend[2];

  int localsize  = Grid.lsize[1]*Grid.lsize[2];
   
  float* iobuf;
  float* iosum;

  char filename[128];

  double q1,q2;

  static int nslice; 
  static double* tau_lev;
  static int nslvar;

  FILE* fhandle=NULL;
	double clk, file_time, dspaces_time, dspaces_wait_time;
	file_time = 0.0;
	dspaces_put_req_t* dspaces_put_req_list;
	if(Run.use_dspaces_io) {
		dspaces_time = 0.0;
    dspaces_wait_time = 0.0;
    dspaces_put_req_list = NULL;
	}

  //MPI_File fhandle_mpi;
  //int offset;
 
  if (ini_flag ==1){
    nslice=Physics.slice[i_sl_tau];
    tau_lev = (realtype*) malloc(nslice*sizeof(realtype));
    for (i=0;i<nslice;i++){
      tau_lev[i] = Physics.tau_lev[i];
    }
    if (Run.rank == 0) {
      cout << "tau_slice: " << nslice << endl;
      for (i=0;i<nslice;i++){
	cout << "tau_slice: " << tau_lev[i]<< endl;
      } 
    }   

    nslvar = 0;
    for (v=0;v<14;v++){
      if (Physics.tau_var[v] == 1) nslvar+=1;
    }

		io_file_log->tau = (struct log_entry*) malloc(sizeof(struct log_entry));
    log_entry_init(io_file_log->tau, "TAU", est_total_slice_iters);
		if(Run.use_dspaces_io) {
      io_dspaces_log->tau = (struct log_entry*) malloc(sizeof(struct log_entry));
      log_entry_init(io_dspaces_log->tau, "TAU", est_total_slice_iters);
    }

    ini_flag = 0;
  }

  iobuf = (float*) malloc(nslvar*localsize*sizeof(float));
  iosum = (float*) malloc(nslvar*localsize*sizeof(float));

  for (nsl = 0; nsl<nslice; nsl++){
 
    for(v=0;v<nslvar*localsize;v++){
      iobuf[v] = 0.0;
      iosum[v] = 0.0;
    }

    for (k=kbeg; k<=kend; k++){
      for (j=jbeg; j<=jend; j++){
	ind  = j-jbeg + (k-kbeg)*Grid.lsize[1];
        for (i=ibeg; i<=iend; i++){
	  node1 = Grid.node(i,j,k);
	  node2 = Grid.node(i+1,j,k);
	  
	  if( (Grid.Tau[node1] >= tau_lev[nsl]) && (tau_lev[nsl] > Grid.Tau[node2]) ){
	    q1 = (tau_lev[nsl]-Grid.Tau[node2])/(Grid.Tau[node1]-Grid.Tau[node2]);
	    //q1 = (log(tau_lev[nsl])-log(Grid.Tau[node2]))/(log(Grid.Tau[node1])-log(Grid.Tau[node2]));
	    q2 = 1.0-q1;
	    
	    ind1 = ind;
	    
	    if (Physics.tau_var[0] == 1){
	      iobuf[ind1] = (float) (Grid.U[node1].d*q1+Grid.U[node2].d*q2);
	      ind1 += localsize;
	    }
	    if (Physics.tau_var[1] == 1){
	      iobuf[ind1] = (float) (Grid.U[node1].M.x*q1+Grid.U[node2].M.x*q2);
	      ind1 += localsize;
	    }
	    if (Physics.tau_var[2] == 1){
	      iobuf[ind1] = (float) (Grid.U[node1].M.y*q1+Grid.U[node2].M.y*q2); 
	      ind1 += localsize;
	    }
	    if (Physics.tau_var[3] == 1){
	      iobuf[ind1] = (float) (Grid.U[node1].M.z*q1+Grid.U[node2].M.z*q2);
	      ind1 += localsize;
	    }
	    if (Physics.tau_var[4] == 1){
	      iobuf[ind1] = (float) (Grid.U[node1].e/Grid.U[node1].d*q1+Grid.U[node2].e/Grid.U[node2].d*q2);
	      ind1 += localsize;
	    }
	    if (Physics.tau_var[5] == 1){
	      iobuf[ind1] = (float) (Grid.U[node1].B.x*q1+Grid.U[node2].B.x*q2);
	      ind1 += localsize;
	    }
	    if (Physics.tau_var[6] == 1){
	      iobuf[ind1] = (float) (Grid.U[node1].B.y*q1+Grid.U[node2].B.y*q2);  
	      ind1 += localsize;
	    }
	    if (Physics.tau_var[7] == 1){
	      iobuf[ind1] = (float) (Grid.U[node1].B.z*q1+Grid.U[node2].B.z*q2);
	      ind1 += localsize;
	    }
	    if (Physics.tau_var[8] == 1){
	      iobuf[ind1] = (float) (sqrt(Grid.U[node1].M.sqr())*q1+sqrt(Grid.U[node2].M.sqr())*q2);
	      ind1 += localsize;
	    }
	    if (Physics.tau_var[9] == 1){
	      iobuf[ind1] = (float) (sqrt(Grid.U[node1].B.sqr())*q1+sqrt(Grid.U[node2].B.sqr())*q2);
	      ind1 += localsize;
	    }
	    if (Physics.tau_var[10] == 1){
	      iobuf[ind1] = (float) (Grid.temp[node1]*q1+Grid.temp[node2]*q2);
	      ind1 += localsize;
	    }
	    if (Physics.tau_var[11] == 1){
	      iobuf[ind1] = (float) (Grid.pres[node1]*q1+Grid.pres[node2]*q2);
	      ind1 += localsize;
	    }
	    if (Physics.tau_var[12] == 1){
	      iobuf[ind1] = (float) rts->Iout(j,k);
	      ind1 += localsize;
	    }
	    if (Physics.tau_var[13] == 1){
	      iobuf[ind1] = float(Grid.coord(i,0)*q1+Grid.coord(i+1,0)*q2)/float(Grid.gxmax[0]);
	    }
	  }
	}
      }
    }

	// check dspaces_iput() except for the first iter
    if(xcol_rank == iroot && Run.use_dspaces_io && nsl > 0) {
	  double dspaces_overlap_time = MPI_Wtime() - clk;
      clk = MPI_Wtime();
      for(int i=0; i<nslvar; i++) {
        dspaces_check_put(ds_client, dspaces_put_req_list[i], 1);
      }
	  double dspaces_wait_time = MPI_Wtime() - clk;
	  if(dspaces_wait_time > 1e-6) {
		dspaces_wait_time += dspaces_wait_time + dspaces_overlap_time;
	  }
      
      free(dspaces_put_req_list);
    }

	// update iosum values
    MPI_Reduce(iobuf,iosum,nslvar*localsize,MPI_FLOAT,MPI_SUM,iroot,
		  XCOL_COMM);

    if (xcol_rank == iroot){

      /*
      sprintf(filename,"%s_%.3f.%06d","tau_slice",tau_lev[nsl],
	      Run.globiter);
      MPI_File_open(YZ_COMM,filename,MPI_MODE_CREATE | MPI_MODE_WRONLY,
		    MPI_INFO_NULL,&fhandle_mpi);
      
      if(yz_rank == 0) {
	float header[4];            
	header[0] = (float) nslvar;
	header[1] = (float) Grid.gsize[1];
	header[2] = (float) Grid.gsize[2];
	header[3] = (float) Run.time;
	MPI_File_write(fhandle_mpi,header,4,MPI_FLOAT,MPI_STATUS_IGNORE);
      }
      
      offset = 4*sizeof(float);
      
      MPI_File_set_view(fhandle_mpi,offset,MPI_FLOAT,yz_subarray,"native",
			MPI_INFO_NULL);
      MPI_File_write_all(fhandle_mpi,&(iosum[0]),nslvar*localsize,
			 MPI_FLOAT,MPI_STATUS_IGNORE);
      
      MPI_File_close(&fhandle_mpi);   
      */

      if(Physics.slice[i_sl_collect] == 0) {	
	if(yz_rank == 0) {
	  if(tau_lev[nsl] >= 1e-3)
	    sprintf(filename,"%s%s_%.3f.%06d",Run.path_2D,"tau_slice",
		    tau_lev[nsl],Run.globiter);
	  else
	    sprintf(filename,"%s%s_%.6f.%06d",Run.path_2D,"tau_slice",
		    tau_lev[nsl],Run.globiter);
	  fhandle=fopen(filename,"w");
	  
	  float header[4];            
	  header[0] = (float) nslvar;
	  header[1] = (float) Grid.gsize[1];
	  header[2] = (float) Grid.gsize[2];
	  header[3] = (float) Run.time;
	  fwrite(header,sizeof(float),4,fhandle);
	}
	
	clk = MPI_Wtime();
	slice_write(Grid,0,&(iosum[0]),localsize,nslvar,1,2,fhandle);
	file_time += MPI_Wtime() - clk;
	
	if(yz_rank == 0)
	  fclose(fhandle);
	
      } else {
	if(yz_rank == 0){
	  sprintf(filename,"%s_%.3f.dat","tau_slice",tau_lev[nsl]);
	  fhandle=fopen(filename,"a");
	}
	
	clk = MPI_Wtime();
	slice_write(Grid,0,&(iosum[0]),localsize,nslvar,1,2,fhandle);
	file_time += MPI_Wtime() - clk;
	
	if(yz_rank == 0){ 
	  fclose(fhandle);
	  
	  sprintf(filename,"%s_%.3f.log","tau_slice",tau_lev[nsl]);
	  
	  fstream fptr;
	  int newfile = 0;
	  fptr.open(filename,ios::in);
	  if (!fptr) newfile = 1;
	  fptr.close();
	  
	  fptr.open(filename,ios::out|ios::app);
	  fptr.precision(10);
	  if (newfile) {      
	    fptr <<  nslvar << ' ' <<  Grid.gsize[1] << ' ' 
		 << Grid.gsize[2] << endl;
	    fptr << Physics.tau_var[0]  << ' ' 
		 << Physics.tau_var[1]  << ' ' 
		 << Physics.tau_var[2]  << ' ' 
		 << Physics.tau_var[3]  << ' ' 
		 << Physics.tau_var[4]  << ' ' 
		 << Physics.tau_var[5]  << ' ' 
		 << Physics.tau_var[6]  << ' ' 
		 << Physics.tau_var[7]  << ' ' 
		 << Physics.tau_var[8]  << ' ' 
		 << Physics.tau_var[9]  << ' ' 
		 << Physics.tau_var[10] << ' ' 
		 << Physics.tau_var[11] << ' '
		 << Physics.tau_var[12] << ' '        
		 << Physics.tau_var[13] << endl;
	  }
	  fptr << Run.globiter << ' ' << Run.time << endl;
	  fptr.close();
	}
      }

		if(Run.use_dspaces_io) {
        char ds_var_name[128];
				if(tau_lev[nsl] >= 1e-3)
        	sprintf(ds_var_name, "%s%s_%.3f", Run.path_2D,"tau_slice",tau_lev[nsl]);
				else
					sprintf(ds_var_name, "%s%s_%.6f", Run.path_2D,"tau_slice",tau_lev[nsl]);
        clk = MPI_Wtime();
        dspaces_put_req_list = slice_write_dspaces(Grid, 0, &(iosum[0]), localsize, nslvar, 1, 2,
													ds_var_name, Run.globiter, 2);
        dspaces_time += MPI_Wtime() - clk;
        char header_filename[128];
        if(yz_rank == 0) {
          sprintf(header_filename, "%s.header", ds_var_name);
          fstream fptr;
          int newfile = 0;
          fptr.open(filename,ios::in);
          if (!fptr) newfile = 1;
          fptr.close();
      
          fptr.open(filename,ios::out|ios::app);
          fptr.precision(10);
          if (newfile) {      
	    			fptr <<  nslvar << ' ' <<  Grid.gsize[1] << ' ' 
		 				<< Grid.gsize[2] << endl;
	    			fptr << Physics.tau_var[0]  << ' ' 
		 				<< Physics.tau_var[1]  << ' ' 
		 				<< Physics.tau_var[2]  << ' ' 
		 				<< Physics.tau_var[3]  << ' ' 
		 				<< Physics.tau_var[4]  << ' ' 
		 				<< Physics.tau_var[5]  << ' ' 
		 				<< Physics.tau_var[6]  << ' ' 
		 				<< Physics.tau_var[7]  << ' ' 
		 				<< Physics.tau_var[8]  << ' ' 
		 				<< Physics.tau_var[9]  << ' ' 
		 				<< Physics.tau_var[10] << ' '
		 				<< Physics.tau_var[11] << endl;
	  			}
          fptr << Run.globiter << ' ' << Run.time << endl;
          fptr.close(); 
        }
		clk = MPI_Wtime();
      }

    }
  }

  free(iobuf);
  if(xcol_rank == iroot && Run.use_dspaces_io) {
    for(int i=0; i<nslvar; i++) {
      dspaces_check_put(ds_client, dspaces_put_req_list[i], 1);
    }
    dspaces_wait_time += MPI_Wtime() - clk;
    free(dspaces_put_req_list);
  }
  free(iosum);

	if(xcol_rank == 0 && yz_rank == 0) {
		io_file_log->tau->iter[io_file_log->tau->count] = Run.globiter;
		io_file_log->tau->api_time[io_file_log->tau->count] = file_time;
		io_file_log->tau->time[io_file_log->tau->count] = file_time;
		io_file_log->tau->count++ ;
		if(Run.use_dspaces_io) {
			io_dspaces_log->tau->iter[io_dspaces_log->tau->count] = Run.globiter;
      io_dspaces_log->tau->api_time[io_dspaces_log->tau->count] = dspaces_time;
      io_dspaces_log->tau->wait_time[io_dspaces_log->tau->count] = dspaces_wait_time;
      io_dspaces_log->tau->time[io_dspaces_log->tau->count] = dspaces_time+dspaces_wait_time;
			io_dspaces_log->tau->count++ ;
		}
		if(Run.verbose >0) {
			std::cout << "File Output (TAU_SLICE) in " << file_time << " seconds" << std::endl;
			if(Run.use_dspaces_io) {
    		std::cout << "DataSpaces API Call (TAU_SLICE) in " << dspaces_time
									<< " seconds" << std::endl;
    		std::cout << "DataSpaces Wait (TAU_SLICE) in " << dspaces_wait_time
									<< " seconds" << std::endl;
    		std::cout << "DataSpaces Output (TAU_SLICE) in " << dspaces_time+dspaces_wait_time
									<< " seconds" << std::endl;
			}
		}
	}
}

