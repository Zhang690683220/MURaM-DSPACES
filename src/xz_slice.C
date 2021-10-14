#include <mpi.h>
#include <math.h>
#include <stdlib.h>
#include <fstream>
#include <stdio.h>
#include "physics.H"
#include "grid.H"
#include "run.H"
#include "comm_split.H"
#include "io.H"

using namespace std;

extern est_total_slice_iters;
extern struct log *io_file_log, *io_dspaces_log;

extern void slice_write(const GridData&,const int,float*,int,int,const int,
			const int,FILE*);

extern dspaces_put_req_t* slice_write_dspaces(const GridData& Grid, const int iroot,
                                float* vloc, int nloc, int nvar,int n0,
                                int n1, char* filename, const int iter,
                                const int ndim);

int xz_dspaces_bufnum = 1;
float **xzslice_buf = NULL;
int xzslice_nslice;
int xzslice_nslvar;
dspaces_put_req_t*** xzslice_dspaces_put_req_list = NULL;

//======================================================================
void xz_slice(const RunData&  Run, const GridData& Grid, 
	      const PhysicsData& Physics) {

  static int ini_flag = 1;
	static int xzslice_ref_count = 0;

  register int i, j, k, node, ind, nsl, v;

  int ibeg = Grid.lbeg[0], iend = Grid.lend[0];
  int kbeg = Grid.lbeg[2], kend = Grid.lend[2];

  int localsize  = Grid.lsize[0]*Grid.lsize[2];
   
  float* iobuf;

  char filename[128];

  static int nslice; 
  static int* ixpos;
  static int nslvar;

  FILE* fhandle=NULL;
	int gsize[2];
	double clk, file_time, dspaces_time, dspaces_wait_time;
	file_time = 0.0;
	int bufind;
	// dspaces_put_req_t* dspaces_put_req_list;
	if(Run.use_dspaces_io) {
		dspaces_time = 0.0;
    dspaces_wait_time = 0.0;
    // dspaces_put_req_list = NULL;
		bufind = xzslice_ref_count % xz_dspaces_bufnum;
	}

  //MPI_File fhandle_mpi;
  //int offset;  

  if (ini_flag ==1){
    nslice=Physics.slice[i_sl_xz];
    ixpos = (int*) malloc(nslice*sizeof(int));
    for (i=0;i<nslice;i++){
      ixpos[i] = Physics.xz_lev[i];
    }

    if (Run.rank == 0) {
      if(Run.verbose > 1)cout << "xz_slice: " << nslice << endl;
      for (i=0;i<nslice;i++){
      if(Run.verbose > 1) cout << "xz_slice: " << ixpos[i]<< endl;
      }     
    }

    nslvar = 0;
    for (v=0;v<12;v++){
      if (Physics.xz_var[v] == 1) nslvar+=1;
    }

    gsize[0] = Grid.gsize[2];
    gsize[1] = Grid.gsize[0];
		io_file_log->xz = (struct log_entry*) malloc(sizeof(struct log_entry));
    log_entry_init(io_file_log->xz, "XZ", est_total_slice_iters, 2, gsize, nslice*nslvar);
		if(Run.use_dspaces_io) {
      io_dspaces_log->xz = (struct log_entry*) malloc(sizeof(struct log_entry));
      log_entry_init(io_dspaces_log->xz, "XZ", est_total_slice_iters, 2, gsize, nslice*nslvar);

			xzslice_nslice = nslice;
			xzslice_nslvar = nslvar;
			xzslice_dspaces_put_req_list = (dspaces_put_req_t***) malloc(xz_dspaces_bufnum *
																															sizeof(dspaces_put_req_t**));
			xzslice_buf = (float**) malloc(xz_dspaces_bufnum*sizeof(float*));
			for(int j=0; j<xz_dspaces_bufnum; j++) {
				xzslice_dspaces_put_req_list[j] = (dspaces_put_req_t**) malloc(nslice*sizeof(dspaces_put_req_t*));
			// prevent non-NULL pointer exists when the rank is not in the selected domain
      	for(int i=0; i<nslice; i++) {
        	xzslice_dspaces_put_req_list[j][i] = NULL;
      	}
    		xzslice_buf[j] = (float*) malloc(nslice*nslvar*localsize*sizeof(float));
			}
    }

    ini_flag = 0;
  }

  iobuf = (float*) malloc(nslvar*localsize*sizeof(float));

  for (nsl = 0; nsl<nslice; nsl++){

    if ( (Grid.beg[1] <= ixpos[nsl]+Grid.gbeg[1] ) and 
	 (Grid.end[1] >= ixpos[nsl]+Grid.gbeg[1] )){

		 // check dspaces_iput() except for the first iter
      // if(Run.use_dspaces_io && nsl > 0) {
      //   for(int i=0; i<nslvar; i++) {
      //     dspaces_check_put(ds_client, dspaces_put_req_list[i], 1);
      //   }
      //   dspaces_wait_time += MPI_Wtime() - clk;
      //   free(dspaces_put_req_list);
      // }


			if(Run.use_dspaces_io && xzslice_ref_count > xz_dspaces_bufnum-1) {
				// int reqind = (xzslice_ref_count-1) % dspaces_bufnum;
				clk = MPI_Wtime();
        for(int i=0; i<nslvar; i++) {
          dspaces_check_put(ds_client, xzslice_dspaces_put_req_list[bufind][nsl][i], 1);
        }
				double dspaces_check_time = MPI_Wtime() - clk;
				if(dspaces_check_time > nslvar*dspaces_check_overhead) {
					dspaces_wait_time += MPI_Wtime() - clk - nslvar*dspaces_check_overhead;
				}
        free(xzslice_dspaces_put_req_list[bufind][nsl]);
      }

			if(Run.use_dspaces_io) {

			// update iobuf values
      for (i=ibeg; i<=iend; i++){
	for (k=kbeg; k<=kend; k++){
	  ind  = k-kbeg + (i-ibeg)*Grid.lsize[2];
	  j    = Grid.lbeg[1]+ixpos[nsl]+Grid.gbeg[1]-Grid.beg[1];
	  node = Grid.node(i,j,k);
	  
	  if (Physics.xz_var[0] == 1){
	    iobuf[ind] = (float) Grid.U[node].d;
			xzslice_buf[bufind][ind+nsl*nslvar*localsize] = (float) Grid.U[node].d;
	    ind += localsize;
	  }
	  if (Physics.xz_var[1] == 1){
	    iobuf[ind] = (float) Grid.U[node].M.x;
			xzslice_buf[bufind][ind+nsl*nslvar*localsize] = (float) Grid.U[node].M.x;
	    ind += localsize;
	  }
	  if (Physics.xz_var[2] == 1){
	    iobuf[ind] = (float) Grid.U[node].M.y; 
			xzslice_buf[bufind][ind+nsl*nslvar*localsize] = (float) Grid.U[node].M.y; 
	    ind += localsize;
	  }
	  if (Physics.xz_var[3] == 1){
	    iobuf[ind] = (float) Grid.U[node].M.z;
			xzslice_buf[bufind][ind+nsl*nslvar*localsize] = (float) Grid.U[node].M.z;
	    ind += localsize;
	  }
	  if (Physics.xz_var[4] == 1){
	    iobuf[ind] = (float) (Grid.U[node].e/Grid.U[node].d);
			xzslice_buf[bufind][ind+nsl*nslvar*localsize] = (float) (Grid.U[node].e/Grid.U[node].d);
	    ind += localsize;
	  }
	  if (Physics.xz_var[5] == 1){
	    iobuf[ind] = (float) Grid.U[node].B.x;
			xzslice_buf[bufind][ind+nsl*nslvar*localsize] = (float) Grid.U[node].B.x;
	    ind += localsize;
	  }
	  if (Physics.xz_var[6] == 1){
	    iobuf[ind] = (float) Grid.U[node].B.y;  
			xzslice_buf[bufind][ind+nsl*nslvar*localsize] = (float) Grid.U[node].B.y;
	    ind += localsize;
	  }
	  if (Physics.xz_var[7] == 1){
	    iobuf[ind] = (float) Grid.U[node].B.z;
			xzslice_buf[bufind][ind+nsl*nslvar*localsize] = (float) Grid.U[node].B.z;
	    ind += localsize;
	  }
	  if (Physics.xz_var[8] == 1){
	    iobuf[ind] = (float) sqrt(Grid.U[node].M.sqr());
			xzslice_buf[bufind][ind+nsl*nslvar*localsize] = (float) sqrt(Grid.U[node].M.sqr());
	    ind += localsize;
	  }
	  if (Physics.xz_var[9] == 1){
	    iobuf[ind] = (float) sqrt(Grid.U[node].B.sqr());
			xzslice_buf[bufind][ind+nsl*nslvar*localsize] = (float) sqrt(Grid.U[node].B.sqr());
	    ind += localsize;
	  }
	  if (Physics.xz_var[10] == 1){
	    iobuf[ind] = (float) Grid.temp[node];
			xzslice_buf[bufind][ind+nsl*nslvar*localsize] = (float) Grid.temp[node];
	    ind += localsize;
	  }
	  if (Physics.xz_var[11] == 1){
	    iobuf[ind] = (float) Grid.pres[node];
			xzslice_buf[bufind][ind+nsl*nslvar*localsize] = (float) Grid.pres[node];
	  }
	}
      }
	} else {
		// update iobuf values
      for (i=ibeg; i<=iend; i++){
	for (k=kbeg; k<=kend; k++){
	  ind  = k-kbeg + (i-ibeg)*Grid.lsize[2];
	  j    = Grid.lbeg[1]+ixpos[nsl]+Grid.gbeg[1]-Grid.beg[1];
	  node = Grid.node(i,j,k);
	  
	  if (Physics.xz_var[0] == 1){
	    iobuf[ind] = (float) Grid.U[node].d;
	    ind += localsize;
	  }
	  if (Physics.xz_var[1] == 1){
	    iobuf[ind] = (float) Grid.U[node].M.x;
	    ind += localsize;
	  }
	  if (Physics.xz_var[2] == 1){
	    iobuf[ind] = (float) Grid.U[node].M.y; 
	    ind += localsize;
	  }
	  if (Physics.xz_var[3] == 1){
	    iobuf[ind] = (float) Grid.U[node].M.z;
	    ind += localsize;
	  }
	  if (Physics.xz_var[4] == 1){
	    iobuf[ind] = (float) (Grid.U[node].e/Grid.U[node].d);
	    ind += localsize;
	  }
	  if (Physics.xz_var[5] == 1){
	    iobuf[ind] = (float) Grid.U[node].B.x;
	    ind += localsize;
	  }
	  if (Physics.xz_var[6] == 1){
	    iobuf[ind] = (float) Grid.U[node].B.y;  
	    ind += localsize;
	  }
	  if (Physics.xz_var[7] == 1){
	    iobuf[ind] = (float) Grid.U[node].B.z;
	    ind += localsize;
	  }
	  if (Physics.xz_var[8] == 1){
	    iobuf[ind] = (float) sqrt(Grid.U[node].M.sqr());
	    ind += localsize;
	  }
	  if (Physics.xz_var[9] == 1){
	    iobuf[ind] = (float) sqrt(Grid.U[node].B.sqr());
	    ind += localsize;
	  }
	  if (Physics.xz_var[10] == 1){
	    iobuf[ind] = (float) Grid.temp[node];
	    ind += localsize;
	  }
	  if (Physics.xz_var[11] == 1){
	    iobuf[ind] = (float) Grid.pres[node];
	  }
	}
      }
	}

      if(Physics.slice[i_sl_collect] == 0) {
	if(xz_rank == 0) { 
	  sprintf(filename,"%s%s_%04d.%06d",Run.path_2D,"xz_slice",ixpos[nsl],
		  Run.globiter);  
	  fhandle=fopen(filename,"w");
	  
	  float header[4];            
	  header[0] = (float) nslvar;
	  header[1] = (float) Grid.gsize[2];
	  header[2] = (float) Grid.gsize[0];
	  header[3] = (float) Run.time;
	  fwrite(header,sizeof(float),4,fhandle);
	}
	
	clk = MPI_Wtime();
	slice_write(Grid,0,iobuf,localsize,nslvar,2,0,fhandle);
	file_time += MPI_Wtime() - clk;
	
	if(xz_rank == 0) 
	  fclose(fhandle);
	
      } else {
	if(xz_rank == 0){
	  sprintf(filename,"%s_%04d.dat","xz_slice",ixpos[nsl]);
	  fhandle=fopen(filename,"a");
	}
	
	clk = MPI_Wtime();
	slice_write(Grid,0,iobuf,localsize,nslvar,2,0,fhandle);
	file_time += MPI_Wtime() - clk;
	
	if(xz_rank == 0){ 
	  fclose(fhandle);
	  
	  sprintf(filename,"%s_%04d.log","xz_slice",ixpos[nsl]);
	  
	  fstream fptr;
	  int newfile = 0;
	  fptr.open(filename,ios::in);
	  if (!fptr) newfile = 1;
	  fptr.close();
	  
	  fptr.open(filename,ios::out|ios::app);
	  fptr.precision(10);
	  if (newfile) {      
	    fptr <<  nslvar << ' ' <<  Grid.gsize[0] << ' ' 
		 << Grid.gsize[2] << endl;
	    fptr << Physics.xz_var[0]  << ' ' 
		 << Physics.xz_var[1]  << ' ' 
		 << Physics.xz_var[2]  << ' ' 
		 << Physics.xz_var[3]  << ' ' 
		 << Physics.xz_var[4]  << ' ' 
		 << Physics.xz_var[5]  << ' ' 
		 << Physics.xz_var[6]  << ' ' 
		 << Physics.xz_var[7]  << ' ' 
		 << Physics.xz_var[8]  << ' ' 
		 << Physics.xz_var[9]  << ' ' 
		 << Physics.xz_var[10] << ' '
		 << Physics.xz_var[11] << endl;
	  }
	  fptr << Run.globiter << ' ' << Run.time << endl;
	  fptr.close();
	}
      }

	  if(Run.use_dspaces_io) {
        char ds_var_name[128];
        sprintf(ds_var_name, "%s%s_%04d", Run.path_2D,"xz_slice",ixpos[nsl]);
				if(xzslice_ref_count == 0) {
					uint64_t gdim[2];
          for(int i=0; i<2; i++) {
            gdim[i] = gsize[i];
          }
					char vname[128];
					for(int i=0; i<nslvar; i++) {
						sprintf(vname, "%s_%d", ds_var_name, i);
						dspaces_define_gdim(ds_client, vname, 2, gdim);
					}
				}
        clk = MPI_Wtime();
        xzslice_dspaces_put_req_list[bufind][nsl] = slice_write_dspaces(Grid, 0, 
																													&xzslice_buf[bufind][nsl*nslvar*localsize],
																													 			localsize, nslvar, 2, 0, ds_var_name,
																									 							Run.globiter, 2);
        dspaces_time += MPI_Wtime() - clk;
        char header_filename[128];
        if(xz_rank == 0) {
          sprintf(header_filename, "%s.header", ds_var_name);
          fstream fptr;
          int newfile = 0;
          fptr.open(filename,ios::in);
          if (!fptr) newfile = 1;
          fptr.close();
      
          fptr.open(filename,ios::out|ios::app);
          fptr.precision(10);
          if (newfile) {      
	    			fptr <<  nslvar << ' ' <<  Grid.gsize[0] << ' ' 
		 				<< Grid.gsize[2] << endl;
	    			fptr << Physics.xz_var[0]  << ' ' 
		 				<< Physics.xz_var[1]  << ' ' 
		 				<< Physics.xz_var[2]  << ' ' 
		 				<< Physics.xz_var[3]  << ' ' 
		 				<< Physics.xz_var[4]  << ' ' 
		 				<< Physics.xz_var[5]  << ' ' 
		 				<< Physics.xz_var[6]  << ' ' 
		 				<< Physics.xz_var[7]  << ' ' 
		 				<< Physics.xz_var[8]  << ' ' 
		 				<< Physics.xz_var[9]  << ' ' 
		 				<< Physics.xz_var[10] << ' '
		 				<< Physics.xz_var[11] << endl;
	  			}
          fptr << Run.globiter << ' ' << Run.time << endl;
          fptr.close(); 
        }
      }
    }	
  }

	// check if put finish for the last dspaces_iput() before iobuf free
	// if ( (Grid.beg[1] <= ixpos[nsl]+Grid.gbeg[1] ) and 
	//  			(Grid.end[1] >= ixpos[nsl]+Grid.gbeg[1] )){
	// 	if(Run.use_dspaces_io) {
  //     for(int i=0; i<nslvar; i++) {
  //       dspaces_check_put(ds_client, dspaces_put_req_list[i], 1);
  //     }
  //     dspaces_wait_time += MPI_Wtime() - clk;
  //     free(dspaces_put_req_list);
  //   }	
	// }
  free(iobuf);

	if(Run.rank == 0) {
		io_file_log->xz->iter[io_file_log->xz->count] = Run.globiter;
		io_file_log->xz->api_time[io_file_log->xz->count] = file_time;
		io_file_log->xz->time[io_file_log->xz->count] = file_time;
		io_file_log->xz->count++;
		if(Run.use_dspaces_io) {
			io_dspaces_log->xz->iter[io_dspaces_log->xz->count] = Run.globiter;
      io_dspaces_log->xz->api_time[io_dspaces_log->xz->count] = dspaces_time;
      if(io_dspaces_log->xz->count > xz_dspaces_bufnum-1) {
        io_dspaces_log->xz->wait_time[io_dspaces_log->xz->count-xz_dspaces_bufnum] = dspaces_wait_time;
        io_dspaces_log->xz->time[io_dspaces_log->xz->count-xz_dspaces_bufnum] = dspaces_wait_time
                                + io_dspaces_log->xz->api_time[io_dspaces_log->xz->count-xz_dspaces_bufnum];
      }
			io_dspaces_log->xz->count++;
		}
		if(Run.verbose > 0) {
			std::cout << "File Output (XZ_SLICE) in " << file_time << " seconds" << std::endl;
			if(Run.use_dspaces_io) {
				std::cout << "DataSpaces API Call (XZ_SLICE) in " << dspaces_time
									<< " seconds" <<  " Bin: " << bufind << std::endl;
    		// std::cout << "DataSpaces Wait (XZ_SLICE) in " << dspaces_wait_time
				// 					<< " seconds" << std::endl;
    		// std::cout << "DataSpaces Output (XZ_SLICE) in " << dspaces_time
				// 					<< " seconds" << std::endl;
			}
		}
	}
	xzslice_ref_count++;
}


