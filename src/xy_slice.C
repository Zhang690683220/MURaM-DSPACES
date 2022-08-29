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

extern int est_total_slice_iters;
extern struct log *io_file_log, *io_dspaces_log;

extern void slice_write(const GridData&,const int,float*,int,int,const int,
			const int,FILE*);

extern void slice_write_dspaces(const GridData& Grid, const int iroot,
                                float* vloc, int nloc, int nvar,int n0,
                                int n1, char* filename, const int iter,
                                const int ndim);

void xy_slice(const RunData&  Run, const GridData& Grid, 
	      const PhysicsData& Physics) {

	static int ini_flag = 1;

	register int i, j, k, node, ind, nsl, v;

	int ibeg = Grid.lbeg[0], iend = Grid.lend[0];
	int jbeg = Grid.lbeg[1], jend = Grid.lend[1];

	int localsize  = Grid.lsize[0]*Grid.lsize[1];
	
	float* iobuf;

	char filename[128];

	static int nslice; 
	static int* ixpos;
	static int nslvar;

	FILE* fhandle=NULL;
	int gsize[2];
	double clk, file_time, dspaces_time;
	file_time = 0.0;
	if(Run.use_dspaces_io) {
		dspaces_time = 0.0;
	}

	//MPI_File fhandle_mpi;
	//int offset;
  
	if(ini_flag ==1) {
		nslice=Physics.slice[i_sl_xy];
		ixpos = (int*) malloc(nslice*sizeof(int));
		for(i=0;i<nslice;i++) {
			ixpos[i] = Physics.xy_lev[i];
		}

		if(Run.rank == 0) {
			if(Run.verbose > 1) {
				cout << "xy_slice: " << nslice << endl;
			}
			for(i=0;i<nslice;i++) {
				if(Run.verbose > 1) {
					cout << "xy_slice: " << ixpos[i]<< endl;
				}
			}     
		}

		nslvar = 0;
		for(v=0;v<12;v++) {
			if(Physics.xy_var[v] == 1) {
				nslvar+=1;
			}
		}

		gsize[0] = Grid.gsize[1];
    	gsize[1] = Grid.gsize[0];

		io_file_log->xy = (struct log_entry*) malloc(sizeof(struct log_entry));
		log_entry_init(io_file_log->xy, (char*)"XY", est_total_slice_iters);
		if(Run.use_dspaces_io) {
			io_dspaces_log->xy = (struct log_entry*) malloc(sizeof(struct log_entry));
			log_entry_init(io_dspaces_log->xy, (char*)"XY", est_total_slice_iters);
		}
	}

  	iobuf = (float*) malloc(nslvar*localsize*sizeof(float));

  	for(nsl = 0; nsl<nslice; nsl++) {
    	if((Grid.beg[2] <= ixpos[nsl]+Grid.gbeg[2]) and 
         	(Grid.end[2] >= ixpos[nsl]+Grid.gbeg[2])) {
      		for(i=ibeg; i<=iend; i++) {
        		for(j=jbeg; j<=jend; j++) {
					ind  = j-jbeg + (i-ibeg)*Grid.lsize[1];
					k    = Grid.lbeg[2]+ixpos[nsl]+Grid.gbeg[2]-Grid.beg[2];
					node = Grid.node(i,j,k);

					if(Physics.xy_var[0] == 1) {
						iobuf[ind] = (float) Grid.U[node].d;
						ind += localsize;
					}
					if(Physics.xy_var[1] == 1) {
						iobuf[ind] = (float) Grid.U[node].M.x;
						ind += localsize;
					}
					if(Physics.xy_var[2] == 1) {
						iobuf[ind] = (float) Grid.U[node].M.y; 
						ind += localsize;
					}
					if(Physics.xy_var[3] == 1) {
						iobuf[ind] = (float) Grid.U[node].M.z;
						ind += localsize;
					}
					if(Physics.xy_var[4] == 1) {
						iobuf[ind] = (float) (Grid.U[node].e/Grid.U[node].d);
						ind += localsize;
					}
					if(Physics.xy_var[5] == 1) {
						iobuf[ind] = (float) Grid.U[node].B.x;
						ind += localsize;
					}
					if(Physics.xy_var[6] == 1) {
						iobuf[ind] = (float) Grid.U[node].B.y;  
						ind += localsize;
					}
					if(Physics.xy_var[7] == 1) {
						iobuf[ind] = (float) Grid.U[node].B.z;
						ind += localsize;
					}
					if(Physics.xy_var[8] == 1) {
						iobuf[ind] = (float) sqrt(Grid.U[node].M.sqr());
						ind += localsize;
					}
					if(Physics.xy_var[9] == 1) {
						iobuf[ind] = (float) sqrt(Grid.U[node].B.sqr());
						ind += localsize;
					}
					if(Physics.xy_var[10] == 1) {
						iobuf[ind] = (float) Grid.temp[node];
						ind += localsize;
					}
					if(Physics.xy_var[11] == 1) {
						iobuf[ind] = (float) Grid.pres[node];
					}
				}
      		}

      		if(Physics.slice[i_sl_collect] == 0) {
				if(xy_rank == 0) { 
	  				sprintf(filename,"%s%s_%04d.%06d",Run.path_2D,"xy_slice",ixpos[nsl],
		  			Run.globiter);
	    			fhandle=fopen(filename,"w");
	    
					float header[4];            
					header[0] = (float) nslvar;
					header[1] = (float) Grid.gsize[1];
					header[2] = (float) Grid.gsize[0];
					header[3] = (float) Run.time;
					fwrite(header,sizeof(float),4,fhandle);
				}
				clk = MPI_Wtime();
				slice_write(Grid,0,iobuf,localsize,nslvar,1,0,fhandle);
				file_time += MPI_Wtime() - clk;
				
				if(xy_rank == 0) 
				fclose(fhandle);
      		} else {
				if(xy_rank == 0) { 
					sprintf(filename,"%s_%04d.dat","xy_slice",ixpos[nsl]);
					fhandle=fopen(filename,"a");
				}
				
				clk = MPI_Wtime();
				slice_write(Grid,0,iobuf,localsize,nslvar,1,0,fhandle);
				file_time += MPI_Wtime() - clk;
				
				if(xy_rank == 0) { 
					fclose(fhandle);
					
					sprintf(filename,"%s_%04d.log","xy_slice",ixpos[nsl]);
					
					fstream fptr;
					int newfile = 0;
					fptr.open(filename,ios::in);
					if(!fptr) {
						newfile = 1;
					}
					fptr.close();
					
					fptr.open(filename,ios::out|ios::app);
					fptr.precision(10);
					if(newfile) {      
						fptr <<  nslvar << ' ' <<  Grid.gsize[0] << ' ' 
						<< Grid.gsize[1] << endl;
						fptr << Physics.xy_var[0]  << ' ' 
						<< Physics.xy_var[1]  << ' ' 
						<< Physics.xy_var[2]  << ' ' 
						<< Physics.xy_var[3]  << ' ' 
						<< Physics.xy_var[4]  << ' ' 
						<< Physics.xy_var[5]  << ' ' 
						<< Physics.xy_var[6]  << ' ' 
						<< Physics.xy_var[7]  << ' ' 
						<< Physics.xy_var[8]  << ' ' 
						<< Physics.xy_var[9]  << ' ' 
						<< Physics.xy_var[10] << ' ' 
						<< Physics.xy_var[11] << endl;
					}
					fptr << Run.globiter << ' ' << Run.time << endl;
					fptr.close();
				}
      		}

			if(Run.use_dspaces_io) {
				char ds_var_name[128];
				sprintf(ds_var_name, "%s%s_%04d", Run.path_2D,"xy_slice",ixpos[nsl]);
				if(ini_flag) {
					uint64_t gdim[2];
					for(int d=0; d<2; d++) {
						gdim[d] = gsize[d];
					}
					char vname[128];
					for(int idx=0; idx<nslvar; idx++) {
						sprintf(vname, "%s_%d", ds_var_name, idx);
						dspaces_define_gdim(ds_client, vname, 2, gdim);
					}
				}
				clk = MPI_Wtime();
				slice_write_dspaces(Grid, 0, iobuf, localsize, nslvar, 1, 0, ds_var_name,
									Run.globiter, 2);
				dspaces_time += MPI_Wtime() - clk;
				char header_filename[128];
				if(xy_rank == 0) {
					sprintf(header_filename, "%s.header", ds_var_name);
					fstream fptr;
					int newfile = 0;
					fptr.open(filename,ios::in);
					if(!fptr) {
						newfile = 1;
					}
					fptr.close();
				
					fptr.open(filename,ios::out|ios::app);
					fptr.precision(10);
					if(newfile) {      
						fptr <<  nslvar << ' ' <<  Grid.gsize[0] << ' ' 
							<< Grid.gsize[1] << endl;
						fptr << Physics.xy_var[0]  << ' ' 
							<< Physics.xy_var[1]  << ' ' 
							<< Physics.xy_var[2]  << ' ' 
							<< Physics.xy_var[3]  << ' ' 
							<< Physics.xy_var[4]  << ' ' 
							<< Physics.xy_var[5]  << ' ' 
							<< Physics.xy_var[6]  << ' ' 
							<< Physics.xy_var[7]  << ' ' 
							<< Physics.xy_var[8]  << ' ' 
							<< Physics.xy_var[9]  << ' ' 
							<< Physics.xy_var[10] << ' ' 
							<< Physics.xy_var[11] << endl;
					}
					fptr << Run.globiter << ' ' << Run.time << endl;
					fptr.close(); 
				}
				clk = MPI_Wtime();
			}
    	}   
  	}


  	free(iobuf);

	if(Run.rank == 0) {
		io_file_log->xy->iter[io_file_log->xy->count] = Run.globiter;
		io_file_log->xy->api_time[io_file_log->xy->count] = file_time;
		io_file_log->xy->time[io_file_log->xy->count] = file_time;
		io_file_log->xy->count++ ;
		if(Run.use_dspaces_io) {
			io_dspaces_log->xy->iter[io_dspaces_log->xy->count] = Run.globiter;
			io_dspaces_log->xy->api_time[io_dspaces_log->xy->count] = dspaces_time;
			io_dspaces_log->xy->wait_time[io_dspaces_log->xy->count] = 0.0;
			io_dspaces_log->xy->time[io_dspaces_log->xy->count] = dspaces_time;
			io_dspaces_log->xy->count++ ;
		}
		if(Run.verbose > 0) {
			std::cout << "File Output (XY_SLICE) in " << file_time << " seconds" << std::endl;
			if(Run.use_dspaces_io) {
    			std::cout << "DataSpaces Output (XY_SLICE) in " << dspaces_time
							<< " seconds" << std::endl;
			}
		}
	}

	if(ini_flag) {
		ini_flag = 0;
	}
}

