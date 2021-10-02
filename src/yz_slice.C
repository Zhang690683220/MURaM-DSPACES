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

extern est_total_slice_iters;
extern struct log *io_file_log, *io_dspaces_log;

extern void slice_write(const GridData&,const int,float*,int,int,const int,
            const int,FILE*);

extern dspaces_put_req_t* slice_write_dspaces(const GridData& Grid, const int iroot,
                                float* vloc, int nloc, int nvar,int n0,
                                int n1, char* filename, const int iter,
                                const int ndim);

//======================================================================
void yz_slice(const RunData&  Run, const GridData& Grid, 
          const PhysicsData& Physics,RTS *rts) {

  static int ini_flag = 1;

  register int i, j, k, node, ind, nsl, v;

  int jbeg = Grid.lbeg[1], jend = Grid.lend[1];
  int kbeg = Grid.lbeg[2], kend = Grid.lend[2];

  int localsize  = Grid.lsize[1]*Grid.lsize[2];
   
  float* iobuf;

  char filename[128];

  static int nslice; 
  static int* ixpos;
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
    nslice=Physics.slice[i_sl_yz];
    ixpos = (int*) malloc(nslice*sizeof(int));
    for (i=0;i<nslice;i++){
      ixpos[i] = Physics.yz_lev[i];
    }

    if (Run.rank == 0) {
      if(Run.verbose > 1)cout << "yz_slice: " << nslice << endl;
      for (i=0;i<nslice;i++){
      if(Run.verbose > 1) cout << "yz_slice: " << ixpos[i]<< endl;
      }     
    }

    nslvar = 0;
    for (v=0;v<13;v++)
      if (Physics.yz_var[v] == 1) nslvar+=1;

    io_file_log->yz = (struct log_entry*) malloc(sizeof(struct log_entry));
    log_entry_init(io_file_log->yz, "YZ", est_total_slice_iters);
		if(Run.use_dspaces_io) {
      io_dspaces_log->yz = (struct log_entry*) malloc(sizeof(struct log_entry));
      log_entry_init(io_dspaces_log->yz, "YZ", est_total_slice_iters);
    }

    ini_flag = 0;
  }

  iobuf = (float*) malloc(nslvar*localsize*sizeof(float));

  for (nsl = 0; nsl<nslice; nsl++){

    if ( (Grid.beg[0] <= ixpos[nsl]+Grid.gbeg[0] ) and 
         (Grid.end[0] >= ixpos[nsl]+Grid.gbeg[0] )){

      // check dspaces_iput() except for the first iter
      if(Run.use_dspaces_io && nsl > 0) {
        for(int i=0; i<nslvar; i++) {
          dspaces_check_put(ds_client, dspaces_put_req_list[i], 1);
        }
        dspaces_wait_time += MPI_Wtime() - clk;
        free(dspaces_put_req_list);
      }

      // update iobuf values
      for (j=jbeg; j<=jend; j++)
        for (k=kbeg; k<=kend; k++){
          ind  = j-jbeg + (k-kbeg)*Grid.lsize[1];
          i    = Grid.lbeg[0]+ixpos[nsl]+Grid.gbeg[0]-Grid.beg[0];
          node = Grid.node(i,j,k);  
       
      if (Physics.yz_var[0] == 1){
        iobuf[ind] = (float) Grid.U[node].d;
        ind += localsize;
      }
      if (Physics.yz_var[1] == 1){
        iobuf[ind] = (float) Grid.U[node].M.x;
        ind += localsize;
      }
      if (Physics.yz_var[2] == 1){
        iobuf[ind] = (float) Grid.U[node].M.y; 
        ind += localsize;
      }
      if (Physics.yz_var[3] == 1){
        iobuf[ind] = (float) Grid.U[node].M.z;
        ind += localsize;
      }
      if (Physics.yz_var[4] == 1){
        iobuf[ind] = (float) (Grid.U[node].e/Grid.U[node].d);
        ind += localsize;
      }
      if (Physics.yz_var[5] == 1){
        iobuf[ind] = (float) Grid.U[node].B.x;
        ind += localsize;
      }
      if (Physics.yz_var[6] == 1){
        iobuf[ind] = (float) Grid.U[node].B.y;  
        ind += localsize;
      }
      if (Physics.yz_var[7] == 1){
        iobuf[ind] = (float) Grid.U[node].B.z;
        ind += localsize;
      }
      if (Physics.yz_var[8] == 1){
        iobuf[ind] = (float) sqrt(Grid.U[node].M.sqr());
        ind += localsize;
      }
      if (Physics.yz_var[9] == 1){
        iobuf[ind] = (float) sqrt(Grid.U[node].B.sqr());
        ind += localsize;
      }
      if (Physics.yz_var[10] == 1){
        iobuf[ind] = (float) Grid.temp[node];
        ind += localsize;
      }
      if (Physics.yz_var[11] == 1){
        iobuf[ind] = (float) Grid.pres[node];
        ind += localsize;
      }
      if (Physics.yz_var[12] == 1){
        iobuf[ind] = (float) rts->Iout(j,k);
      }
      }
      

      if(Physics.slice[i_sl_collect] == 0) {
    if(yz_rank == 0) { 
      sprintf(filename,"%s%s_%04d.%06d",Run.path_2D,"yz_slice",ixpos[nsl],
          Run.globiter);
      fhandle=fopen(filename,"w");
      
      float header[4];            
      header[0] = (float) nslvar;
      header[1] = (float) Grid.gsize[1];
      header[2] = (float) Grid.gsize[2];
      header[3] = (float) Run.time;
      fwrite(header,sizeof(float),4,fhandle);
    }
    
    clk = MPI_Wtime();
    slice_write(Grid,0,iobuf,localsize,nslvar,1,2,fhandle);
    file_time += MPI_Wtime() - clk;

    if(yz_rank == 0)
      fclose(fhandle);
    
      } else {
    // collective means write all iterations into the same file
    if(yz_rank == 0) { 
      sprintf(filename,"%s_%04d.dat","yz_slice",ixpos[nsl]);
      fhandle=fopen(filename,"a");
    }
    
    clk = MPI_Wtime();
    slice_write(Grid,0,iobuf,localsize,nslvar,1,2,fhandle);
    file_time += MPI_Wtime() - clk;
    
    if(yz_rank == 0){
      fclose(fhandle);
      
      sprintf(filename,"%s_%04d.log","yz_slice",ixpos[nsl]);
      
      fstream fptr;
      int newfile = 0;
      fptr.open(filename,ios::in);
      if (!fptr) newfile = 1;
      fptr.close();
      
      fptr.open(filename,ios::out|ios::app);
      fptr.precision(10);
      if (newfile) {      
        fptr <<  nslvar << ' ' <<  Grid.gsize[2] << ' ' 
         << Grid.gsize[1] << endl;
        fptr << Physics.yz_var[0]  << ' ' 
         << Physics.yz_var[1]  << ' ' 
         << Physics.yz_var[2]  << ' ' 
         << Physics.yz_var[3]  << ' ' 
         << Physics.yz_var[4]  << ' ' 
         << Physics.yz_var[5]  << ' ' 
         << Physics.yz_var[6]  << ' ' 
         << Physics.yz_var[7]  << ' ' 
         << Physics.yz_var[8]  << ' ' 
         << Physics.yz_var[9]  << ' ' 
         << Physics.yz_var[10] << ' ' 
         << Physics.yz_var[11] << ' '
         << Physics.yz_var[12] << endl;
      }
      fptr << Run.globiter << ' ' << Run.time << endl;
      fptr.close(); 
    }
      }

      if(Run.use_dspaces_io) {
        char ds_var_name[128];
        sprintf(ds_var_name, "%s%s_%04d", Run.path_2D,"yz_slice",ixpos[nsl]);
        clk = MPI_Wtime();
        dspaces_put_req_list = slice_write_dspaces(Grid, 0, iobuf, localsize, nslvar, 1, 2, ds_var_name,
                                                   Run.globiter, 2);
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
            fptr <<  nslvar << ' ' <<  Grid.gsize[2] << ' ' 
            << Grid.gsize[1] << endl;
            fptr << Physics.yz_var[0]  << ' ' 
            << Physics.yz_var[1]  << ' ' 
            << Physics.yz_var[2]  << ' ' 
            << Physics.yz_var[3]  << ' ' 
            << Physics.yz_var[4]  << ' ' 
            << Physics.yz_var[5]  << ' ' 
            << Physics.yz_var[6]  << ' ' 
            << Physics.yz_var[7]  << ' ' 
            << Physics.yz_var[8]  << ' ' 
            << Physics.yz_var[9]  << ' ' 
            << Physics.yz_var[10] << ' ' 
            << Physics.yz_var[11] << ' '
            << Physics.yz_var[12] << endl;
          }
          fptr << Run.globiter << ' ' << Run.time << endl;
          fptr.close(); 
        }
        clk = MPI_Wtime();
      }
    }

         
  }

  // check if put finish for the last dspaces_iput() before iobuf free
  if ( (Grid.beg[0] <= ixpos[nsl]+Grid.gbeg[0] ) and 
         (Grid.end[0] >= ixpos[nsl]+Grid.gbeg[0] )){
    if(Run.use_dspaces_io) {
      for(int i=0; i<nslvar; i++) {
        dspaces_check_put(ds_client, dspaces_put_req_list[i], 1);
      }
      dspaces_wait_time += MPI_Wtime() - clk;
      free(dspaces_put_req_list);
    }
  }
  free(iobuf);

  if(Run.rank == 0 && Run.verbose >0) {
    io_file_log->yz->iter[io_file_log->yz->count] = Run.globiter;
		io_file_log->yz->api_time[io_file_log->yz->count] = file_time;
		io_file_log->yz->time[io_file_log->yz->count] = file_time;
    io_file_log->yz->count++;
		if(Run.use_dspaces_io) {
			io_dspaces_log->yz->iter[io_dspaces_log->yz->count] = Run.globiter;
      io_dspaces_log->yz->api_time[io_dspaces_log->yz->count] = dspaces_time;
      io_dspaces_log->yz->wait_time[io_dspaces_log->yz->count] = dspaces_wait_time;
      io_dspaces_log->yz->time[io_dspaces_log->yz->count] = dspaces_time+dspaces_wait_time;
      io_dspaces_log->yz->count++;
		}
    if(Run.verbose > 0) {
		  std::cout << "File Output (YZ_SLICE) in " << file_time << " seconds" << std::endl;
      if(Run.use_dspaces_io) {
        std::cout << "DataSpaces API Call (YZ_SLICE) in " << dspaces_time
                  << " seconds" << std::endl;
        std::cout << "DataSpaces Wait (YZ_SLICE) in " << dspaces_wait_time
                  << " seconds" << std::endl;
        std::cout << "DataSpaces Output (YZ_SLICE) in " << dspaces_time+dspaces_wait_time
                  << " seconds" << std::endl;
      }
    }
  }
}

