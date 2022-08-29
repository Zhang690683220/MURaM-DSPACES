#include <mpi.h>
#include <fstream>
#include <stdio.h>
#include <string.h>
#include "io.H"
#include "precision.h"
#include "grid.H"
#include "run.H"
#include <math.h>
#include "comm_split.H"
#include "rt/rt.h"
#include "eos.H"
#include "limit_va.H"
#include "dfparser.h"
// dataspaces header
#include "dspaces.h"

static int blocksize     = 8;

static int nblocks,blsz;
static MPI_Comm io_xy_comm,io_z_comm,io_comm;

static int ds_terminate = 0;
static dspaces_client_t ndcl = dspaces_CLIENT_NULL;
static uint64_t *lb, *ub, *op_lb, *op_ub;
static int io_rank; // for print msg for optimized IO

static double clk, ds_io_time;

int ds_IO_Init(const GridData& Grid, const RunData& Run) {
    int i;
    int gsz[3]; int lsz[3]; int str[3];

    MPI_Comm_dup(MPI_COMM_WORLD,&io_comm);
    MPI_Comm_dup(XY_COMM,&io_xy_comm);
    MPI_Comm_dup(ZCOL_COMM,&io_z_comm);

    for (i=0;i<2;i++){
        gsz[i]=Grid.gsize[i];
        lsz[i]=Grid.lsize[i];
        str[i]=Grid.beg[i]-Grid.gbeg[i];
    }
    gsz[2] = Grid.gsize[2];
    lsz[2] = Grid.gsize[2];
    str[2] = 0;

    if(Run.use_dspaces_io) {
        char listen_addr_str[128];
        ds_io_time = 0.0;
        if(Run.dspaces_terminate) {
            ds_terminate = 1;
        }
        lb = (uint64_t*) malloc(3*sizeof(uint64_t));
        ub = (uint64_t*) malloc(3*sizeof(uint64_t));
        for(int ii=0; ii<3; ii++) {
            lb[ii] = str[ii];
            ub[ii] = str[ii]+lsz[ii]-1;
        }

        int dspaces_rank = xy_rank;

        if(Run.dspaces_optimized) {
            MPI_Comm_rank(io_comm, &io_rank);
            dspaces_rank = io_rank;
            op_lb = (uint64_t*) malloc(3*sizeof(uint64_t));
            op_ub = (uint64_t*) malloc(3*sizeof(uint64_t));
            for(int ii=0; ii<2; ii++) {
                op_lb[ii] = str[ii];
                op_ub[ii] = str[ii]+lsz[ii]-1;
            }
            op_lb[2] = Grid.beg[2]-Grid.gbeg[2];
            op_ub[2] = op_lb[2] + Grid.lsize[2] - 1;
        }

        dspaces_init(dspaces_rank, &ndcl);

    } else {
        return -1;
    }

    if ( Grid.gsize[2]%blocksize == 0 ){
        nblocks = Grid.gsize[2]/blocksize;
        blsz    = blocksize;
    } else {
        nblocks = Grid.gsize[2];
        blsz    = 1;
    }

    if(xy_rank+zcol_rank == 0)
      std::cout << "xy_slice_io: " 
        << nblocks << ' ' << blsz 
        << std::endl;

    return 0;
}

void ds_IO_Finalize() {
    MPI_Comm_free(&io_comm);
    MPI_Comm_free(&io_xy_comm);
    MPI_Comm_free(&io_z_comm);
    if(ds_terminate) {
        dspaces_kill(ndcl);
    }
    dspaces_fini(ndcl);
}

//=======================================================================
static void z_scatter_io(const GridData& Grid, const int iroot,float* vglo,int nglo,
          float* vloc, int nloc){
  
  int np;
  const int nprocs = Grid.procs[2];
  
  static int* sendcounts;
  static int* offsets;
  static int sendbufsize;
  static int ini_flag = 1;

  if (ini_flag == 1){
    sendcounts  = (int*) malloc(nprocs*sizeof(int));
    offsets     = (int*) malloc(nprocs*sizeof(int)); 
 
    MPI_Allgather(&nloc,1,MPI_INT,sendcounts,1,MPI_INT,io_z_comm);
    
    offsets[0]=0;
    for (np=1;np<nprocs;np++){  
      offsets[np]=offsets[np-1]+sendcounts[np-1];
    }
    sendbufsize=0;
    for (np=0;np<nprocs;np++){
      sendbufsize+=sendcounts[np];
    }

    ini_flag = 0;
  }
    
  if(zcol_rank == iroot){
    if(nglo != sendbufsize){
      std::cout << "z_scatter_io: nglo != sendbufsize " << zcol_rank << ' ' 
       << nglo << ' ' << sendbufsize << std::endl;
      MPI_Abort(io_comm,1);
    }
  }
  
  if(nloc != sendcounts[zcol_rank]){
    std::cout << "z_scatter_io: nloc != sendcounts " << zcol_rank << ' ' 
         << nloc << ' ' << sendcounts[zcol_rank] << std::endl;
    MPI_Abort(io_comm,1);
  }
  
  MPI_Scatterv(vglo,sendcounts,offsets,MPI_FLOAT,vloc,nloc,MPI_FLOAT,
           iroot,io_z_comm);
}

void eos_ds_read(const RunData& Run, const GridData& Grid, const PhysicsData& Physics) {

    char ds_var_name[128];

    register int i,j,k,loc;

    float* iobuf_loc;
    float* iobuf_glo=NULL;

    int sizex=Grid.lend[0]-Grid.lbeg[0]+1;
    int sizey=Grid.lend[1]-Grid.lbeg[1]+1;
    int sizez=Grid.lend[2]-Grid.lbeg[2]+1;

    int lsize=sizex*sizey*sizez; 
    int gsize=sizex*sizey*Grid.gsize[2];

    int v1_max,v2_max,v1,v2, var;
    int max_vars = 14;
    char eos_names[max_vars][128];
    double* eos_vars[max_vars];

    //Only Read in variables that re allocated based on phyasics configuration
    int var_init[max_vars];
    for(var=0;var<max_vars;var++) {
        var_init[var] = 0;
    }
  
    var_init[0] = 1;
    var_init[1] = 1;
    var_init[2] = 1;
    var_init[3] = 1;
    var_init[4] = 1;
    var_init[5] = 1;
    var_init[6] = 1;
    var_init[7] = 1;
    var_init[8] = 1;

    if(Physics.rt_ext[i_ext_cor] >= 1) {
        var_init[9] = 1;
    }
    if(Physics.rt_ext[i_ext_cor] == 2) {
        var_init[10] = 1;
        var_init[11] = 1;
        var_init[12] = 1;
        var_init[13] = 1;
    }

    int var_index[max_vars];

    int tot_vars = 0;
    for(var=0;var<max_vars;var++) {
        if((Run.eos_output[var] == 1) && (var_init[var] == 1)) {
            var_index[tot_vars]=var;
            tot_vars +=1;
        }
    }
      
    if(Grid.procs[2]<tot_vars){
        v1_max = tot_vars;
        v2_max = 1;
    } else {
        v1_max = 1;
        v2_max = tot_vars;
    }

    // This is ugly, but for now it works so I will stick with it -- DP 
    sprintf(eos_names[0],"%s","eosT");
    sprintf(eos_names[1],"%s","eosP");
    sprintf(eos_names[2],"%s","eosne");
    sprintf(eos_names[3],"%s","eosrhoi");
    sprintf(eos_names[4],"%s","eosamb");
    sprintf(eos_names[5],"%s","Qtot");
    sprintf(eos_names[6],"%s","tau");
    sprintf(eos_names[7],"%s","Jtot");
    sprintf(eos_names[8],"%s","Stot");
    sprintf(eos_names[9],"%s","QxCor");
    sprintf(eos_names[10],"%s","QxH");
    sprintf(eos_names[11],"%s","QxMg");
    sprintf(eos_names[12],"%s","QxCa");
    sprintf(eos_names[13],"%s","QxChr");

    eos_vars[0] = Grid.temp;
    eos_vars[1] = Grid.pres;
    eos_vars[2] = Grid.ne;
    eos_vars[3] = Grid.rhoi;
    eos_vars[4] = Grid.amb;
    eos_vars[5] = Grid.Qtot;
    eos_vars[6] = Grid.Tau;
    eos_vars[7] = Grid.Jtot;
    eos_vars[8] = Grid.Stot;
    eos_vars[9] = Grid.Qthin;
    eos_vars[10] = Grid.QH;
    eos_vars[11] = Grid.QMg;
    eos_vars[12] = Grid.QCa;
    eos_vars[13] = Grid.QChr;

    iobuf_loc = (float*)malloc(lsize*sizeof(float));
    for(v2=0;v2<v2_max;v2++) {
        if (zcol_rank == v2) iobuf_glo = (float*)malloc(gsize*sizeof(float));
    }
    

    for(v1=0; v1<v1_max; v1++) {

        for(v2=0; v2<v2_max;v2++) {
            var = var_index[v1*v2_max + v2];
            if(zcol_rank == v2) {
                sprintf(ds_var_name, "%s%s", Run.path_3D,eos_names[var]);
                if(xy_rank == 0) {
                    std::cout << "Read " << ds_var_name
                            << " Version:" << Run.globiter <<std::endl;
                }
                clk = MPI_Wtime();
                int ret = dspaces_get(ndcl, ds_var_name, Run.globiter, 
                                        sizeof(float), 3, lb, ub, iobuf_glo, -1);
                if(ret != 0) {
                    std::cout << "Error Reading " << ds_var_name << "Version: " << Run.globiter
                            << "to DataSpaces Server. Aborting ... " << std::endl;
                MPI_Abort(MPI_COMM_WORLD,1);
                }
                ds_io_time += MPI_Wtime() - clk;

                /*
                if(xy_rank == 0) {
                    std::cout << ds_var_name << " Version:" << Run.globiter << "Data: " <<std::endl;
                    int format =0 ;
                    for(k=0; k<Grid.gsize[2]; k++) {
                        for(j=0; j<sizey; j++) {
                            for(i=0; i<sizex; i++) {
                                std::cout.width(10);
                                std::cout<< std::left << iobuf_glo[i+j*sizex+k*sizex*sizey];
                                format++;
                                if(format%8==0)
                                    std::cout << std::endl;
                            }
                        }
                    }
                }
                */
            }
        }

        for(v2=0; v2<v2_max; v2++) {
            var = var_index[v1*v2_max + v2];
            z_scatter_io(Grid,v2,iobuf_glo,gsize,iobuf_loc,lsize);
            for(k=0; k<sizez; k++) {
                for(j=0; j<sizey; j++) {
                    for(i=0; i<sizex; i++) {
                        loc=Grid.node(i+Grid.ghosts[0],j+Grid.ghosts[1],k+Grid.ghosts[2]);
                        eos_vars[var][loc] = (double) iobuf_loc[i+j*sizex+k*sizex*sizey];
                    }
                }
            }
        }
    }

    if(xy_rank == 0) {
        if(Run.verbose >0) {
            cout << "DataSpaces Output (EOS) in " << ds_io_time << " seconds" << endl;
        } 
    }

    free(iobuf_loc); 
    for(v2=0;v2<v2_max;v2++) {
        if (zcol_rank == v2) {
            free(iobuf_glo);
        }
    }
    
}

void eos_ds_read_optimized(const RunData& Run, const GridData& Grid, const PhysicsData& Physics) {

    char ds_var_name[128];

    register int i,j,k,loc;

    float* iobuf_loc;
    //float* iobuf_glo=NULL;

    int sizex=Grid.lend[0]-Grid.lbeg[0]+1;
    int sizey=Grid.lend[1]-Grid.lbeg[1]+1;
    int sizez=Grid.lend[2]-Grid.lbeg[2]+1;

    int lsize=sizex*sizey*sizez; 
    int gsize=sizex*sizey*Grid.gsize[2];

    int v_max, vi, var;
    int max_vars = 14;
    char eos_names[max_vars][128];
    double* eos_vars[max_vars];

    //Only Read in variables that re allocated based on phyasics configuration
    int var_init[max_vars];
    for(var=0;var<max_vars;var++) {
        var_init[var] = 0;
    }
  
    var_init[0] = 1;
    var_init[1] = 1;
    var_init[2] = 1;
    var_init[3] = 1;
    var_init[4] = 1;
    var_init[5] = 1;
    var_init[6] = 1;
    var_init[7] = 1;
    var_init[8] = 1;

    if(Physics.rt_ext[i_ext_cor] >= 1) {
        var_init[9] = 1;
    }
    if(Physics.rt_ext[i_ext_cor] == 2) {
        var_init[10] = 1;
        var_init[11] = 1;
        var_init[12] = 1;
        var_init[13] = 1;
    }

    int var_index[max_vars];

    int tot_vars = 0;
    for(var=0;var<max_vars;var++) {
        if((Run.eos_output[var] == 1) && (var_init[var] == 1)) {
            var_index[tot_vars]=var;
            tot_vars +=1;
        }
    }
      
    v_max = tot_vars;

    // This is ugly, but for now it works so I will stick with it -- DP 
    sprintf(eos_names[0],"%s","eosT");
    sprintf(eos_names[1],"%s","eosP");
    sprintf(eos_names[2],"%s","eosne");
    sprintf(eos_names[3],"%s","eosrhoi");
    sprintf(eos_names[4],"%s","eosamb");
    sprintf(eos_names[5],"%s","Qtot");
    sprintf(eos_names[6],"%s","tau");
    sprintf(eos_names[7],"%s","Jtot");
    sprintf(eos_names[8],"%s","Stot");
    sprintf(eos_names[9],"%s","QxCor");
    sprintf(eos_names[10],"%s","QxH");
    sprintf(eos_names[11],"%s","QxMg");
    sprintf(eos_names[12],"%s","QxCa");
    sprintf(eos_names[13],"%s","QxChr");

    eos_vars[0] = Grid.temp;
    eos_vars[1] = Grid.pres;
    eos_vars[2] = Grid.ne;
    eos_vars[3] = Grid.rhoi;
    eos_vars[4] = Grid.amb;
    eos_vars[5] = Grid.Qtot;
    eos_vars[6] = Grid.Tau;
    eos_vars[7] = Grid.Jtot;
    eos_vars[8] = Grid.Stot;
    eos_vars[9] = Grid.Qthin;
    eos_vars[10] = Grid.QH;
    eos_vars[11] = Grid.QMg;
    eos_vars[12] = Grid.QCa;
    eos_vars[13] = Grid.QChr;

    iobuf_loc = (float*)malloc(lsize*sizeof(float));
    
    for(vi=0; vi<v_max; vi++) {
        var = var_index[vi];
        sprintf(ds_var_name, "%s%s", Run.path_3D,eos_names[var]);
        if(io_rank == 0) {
            std::cout << "Read " << ds_var_name
                    << " Version:" << Run.globiter <<std::endl;
        }
        clk = MPI_Wtime();
        int ret = dspaces_get(ndcl, ds_var_name, Run.globiter, 
                                sizeof(float), Grid.NDIM, lb, ub, iobuf_loc, -1);
        if(ret != 0) {
            std::cout << "Error Reading" << ds_var_name << "Version: " << Run.globiter
                    << "to DataSpaces Server. Aborting ... " << std::endl;
            MPI_Abort(MPI_COMM_WORLD,1);
        }

        for(k=0; k<sizez; k++) {
            for(j=0; j<sizey; j++) {
                for(i=0; i<sizex; i++) {
                    loc=Grid.node(i+Grid.ghosts[0],j+Grid.ghosts[1],k+Grid.ghosts[2]);
                    eos_vars[var][loc] = (double) iobuf_loc[i+j*sizex+k*sizex*sizey];
                }
            }
        }

    }


    if(io_rank == 0) {
        if(Run.verbose >0) {
            cout << "DataSpaces Output (EOS) in " << ds_io_time << " seconds" << endl;
        } 
    }

    free(iobuf_loc); 
    
}

int eos_compare(const GridData& Grid1, const GridData& Grid2, const RunData& Run,
                const PhysicsData& Physics) {
    char ds_var_name[128];
    register int i,j,k,loc;

    int sizex=Grid1.lend[0]-Grid1.lbeg[0]+1;
    int sizey=Grid1.lend[1]-Grid1.lbeg[1]+1;
    int sizez=Grid1.lend[2]-Grid1.lbeg[2]+1;

    int v_max, vi, var;
    int max_vars = 14;
    char eos_names[max_vars][128];
    double* eos_vars1[max_vars], *eos_vars2[max_vars];

    //Only Read in variables that re allocated based on phyasics configuration
    int var_init[max_vars];
    for(var=0;var<max_vars;var++) {
        var_init[var] = 0;
    }
  
    var_init[0] = 1;
    var_init[1] = 1;
    var_init[2] = 1;
    var_init[3] = 1;
    var_init[4] = 1;
    var_init[5] = 1;
    var_init[6] = 1;
    var_init[7] = 1;
    var_init[8] = 1;

    if(Physics.rt_ext[i_ext_cor] >= 1) {
        var_init[9] = 1;
    }
    if(Physics.rt_ext[i_ext_cor] == 2) {
        var_init[10] = 1;
        var_init[11] = 1;
        var_init[12] = 1;
        var_init[13] = 1;
    }

    int var_index[max_vars];

    int tot_vars = 0;
    for(var=0;var<max_vars;var++) {
        if((Run.eos_output[var] == 1) && (var_init[var] == 1)) {
            var_index[tot_vars]=var;
            tot_vars +=1;
        }
    }
      
    v_max = tot_vars;

    // This is ugly, but for now it works so I will stick with it -- DP 
    sprintf(eos_names[0],"%s","eosT");
    sprintf(eos_names[1],"%s","eosP");
    sprintf(eos_names[2],"%s","eosne");
    sprintf(eos_names[3],"%s","eosrhoi");
    sprintf(eos_names[4],"%s","eosamb");
    sprintf(eos_names[5],"%s","Qtot");
    sprintf(eos_names[6],"%s","tau");
    sprintf(eos_names[7],"%s","Jtot");
    sprintf(eos_names[8],"%s","Stot");
    sprintf(eos_names[9],"%s","QxCor");
    sprintf(eos_names[10],"%s","QxH");
    sprintf(eos_names[11],"%s","QxMg");
    sprintf(eos_names[12],"%s","QxCa");
    sprintf(eos_names[13],"%s","QxChr");

    eos_vars1[0] = Grid1.temp;
    eos_vars1[1] = Grid1.pres;
    eos_vars1[2] = Grid1.ne;
    eos_vars1[3] = Grid1.rhoi;
    eos_vars1[4] = Grid1.amb;
    eos_vars1[5] = Grid1.Qtot;
    eos_vars1[6] = Grid1.Tau;
    eos_vars1[7] = Grid1.Jtot;
    eos_vars1[8] = Grid1.Stot;
    eos_vars1[9] = Grid1.Qthin;
    eos_vars1[10] = Grid1.QH;
    eos_vars1[11] = Grid1.QMg;
    eos_vars1[12] = Grid1.QCa;
    eos_vars1[13] = Grid1.QChr;

    eos_vars2[0] = Grid2.temp;
    eos_vars2[1] = Grid2.pres;
    eos_vars2[2] = Grid2.ne;
    eos_vars2[3] = Grid2.rhoi;
    eos_vars2[4] = Grid2.amb;
    eos_vars2[5] = Grid2.Qtot;
    eos_vars2[6] = Grid2.Tau;
    eos_vars2[7] = Grid2.Jtot;
    eos_vars2[8] = Grid2.Stot;
    eos_vars2[9] = Grid2.Qthin;
    eos_vars2[10] = Grid2.QH;
    eos_vars2[11] = Grid2.QMg;
    eos_vars2[12] = Grid2.QCa;
    eos_vars2[13] = Grid2.QChr;

    for(vi=0; vi<v_max; vi++) {
        var = var_index[vi];
        sprintf(ds_var_name, "%s%s", Run.path_3D,eos_names[var]);
        if(io_rank == 0) {
            std::cout << "Comparing " << ds_var_name
                    << " Version:" << Run.globiter <<std::endl;
        }

        for(k=0; k<sizez; k++) {
            for(j=0; j<sizey; j++) {
                for(i=0; i<sizex; i++) {
                    loc=Grid1.node(i+Grid1.ghosts[0],j+Grid1.ghosts[1],k+Grid1.ghosts[2]);
                    if(eos_vars1[var][loc]!=eos_vars2[var][loc]) {
                        std::cout << "Not same!!!! " << ds_var_name
                            << " Version:" << Run.globiter <<std::endl;
                        MPI_Abort(MPI_COMM_WORLD,1);
                    }   
                }
            }
        }
    }
    std::cout << "Comparison Succeed!!!! " <<std::endl;

    return 0;
}

int Initialize(RunData& Run,GridData& Grid,
	       PhysicsData& Physics) {
  char datafile[256] = "parameters.dat";
  char rtype[16] = "double";
  //int mode = newrun;
  int i,rank;
  if( sizeof(real)==sizeof(float) ) strcpy(rtype,"float");

  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  //ACCH::SetGPU(rank%ACCH::GetNumGPU());

  if( rank==0 ) {
    int flag = 0;
    char ratype[16];
    strcpy(ratype,rtype);

    getvar(Run.anlfile,"anlfile","char*",datafile);
    getvar(Run.resfile,"resfile","char*",datafile);
    getvar(Run.backfile,"backfile","char*",datafile);
    getvar(&Run.maxiter,"maxiter","int",datafile);
    getvar(&Run.anlfreq,"anlfreq","int",datafile);
    getvar(&Run.resfreq,"resfreq","int",datafile);
    getvar(&Run.backfreq,"backfreq","int",datafile);
    getvar(&Run.outcad,"outcad",rtype,datafile);
    getvar(&Run.slicefreq,"slicefreq","int",datafile);
    getvar(&Run.dt,"dt",rtype,datafile);
    getvar(&Run.Tmax,"Tmax",rtype,datafile);
    getvar(&Run.CFL,"CFL",rtype,datafile);
    getvar(&Run.CFL_tvd,"CFL_tvd",rtype,datafile);
    getvar(&Run.maxWtime,"maxWtime",rtype,datafile);
    getvar(&Run.comment,"comment","char*",datafile);
    getvar(&Run.verbose,"verbose","int",datafile);
    getvar(&Run.path_3D,"path_3D","char*",datafile);
    getvar(&Run.path_2D,"path_2D","char*",datafile);
    getvar(&Run.eos_name,"eos_name","char*",datafile);
    getvar(&Run.kap_name,"kap_name","char*",datafile);
    getvar(&Run.eos_output,"eos_output","int[14]",datafile);
    getvar(&Run.diagnostics,"diagnostics","int",datafile);
    getvar(&Run.HAVG,"HAVG","int",datafile);
    getvar(&Run.RT_HAVG,"RT_HAVG","int",datafile);
    getvar(&Run.DEM,"DEM","int",datafile);
    getvar(&Run.diag_output,"diag_output","int[11]",datafile);

    getvar(&Run.use_dspaces_io, "use_dspaces_io", "int", datafile);
    getvar(&Run.dspaces_terminate, "dspaces_terminate", "int", datafile);

    getvar_s(&Grid.NDIM,"NDIM","int",datafile);

    if( Grid.NDIM==3 ) {
      strcat(ratype,"[3]");
      getvar(Grid.gxmin,"gxmin",ratype,datafile);
      getvar(Grid.gxmax,"gxmax",ratype,datafile);
      getvar(Grid.gsize,"gsize","int[3]",datafile);
      getvar(Grid.pardim,"pardim","int[3]",datafile);
      getvar(Grid.periods,"periods","int[3]",datafile);
      getvar(Grid.procs,"procs","int[3]",datafile);
      flag = getvar(Grid.ghosts,"ghosts","int[3]",datafile);
    } else if (Grid.NDIM==2) {
      strcat(ratype,"[2]");
      getvar(Grid.gxmin,"gxmin",ratype,datafile);
      getvar(Grid.gxmax,"gxmax",ratype,datafile);
      getvar(Grid.gsize,"gsize","int[2]",datafile);
      getvar(Grid.pardim,"pardim","int[2]",datafile);
      getvar(Grid.periods,"periods","int[2]",datafile);
      getvar(Grid.procs,"procs","int[2]",datafile);
      flag = getvar(Grid.ghosts,"ghosts","int[2]",datafile);
    } else if (Grid.NDIM==1) {
      strcat(ratype,"[1]");
      getvar(Grid.gxmin,"gxmin",ratype,datafile);
      getvar(Grid.gxmax,"gxmax",ratype,datafile);
      getvar(Grid.gsize,"gsize","int[1]",datafile);
      getvar(Grid.pardim,"pardim","int[1]",datafile);
      getvar(Grid.periods,"periods","int[1]",datafile);
      getvar(Grid.procs,"procs","int[1]",datafile);
      flag = getvar(Grid.ghosts,"ghosts","int[1]",datafile);
    }

    if( flag )
      for(i=0;i<Grid.NDIM;i++) Grid.ghosts[i] = 2;
    
    getvar_s(&Physics.params[i_param_grav], "param_gravity",rtype,datafile);
    getvar_s(&Physics.params[i_param_va_max],"param_va_max",rtype,datafile);
    getvar_s(&Physics.params[i_param_va_adjust],"param_va_adjust",rtype,datafile);
    getvar_s(&Physics.params[i_param_spitzer],"param_spitzer",rtype,datafile);
    getvar(&Physics.params[i_param_eta],"param_eta",rtype,datafile);
    getvar(&Physics.params[i_param_max_fill],"param_max_fill",rtype,datafile);
    getvar(&Physics.params[i_param_ambipolar],"param_ambipolar",rtype,datafile);
    getvar(&Physics.params[i_param_ambfac_max],"param_ambfac_max",rtype,datafile);
    getvar(&Physics.params[i_param_ambvel_max],"param_ambvel_max",rtype,datafile);
    
    getvar_s(&Physics.bnd[i_bnd_top],   "bnd_top",rtype,datafile);
    getvar_s(&Physics.bnd[i_bnd_pot],   "bnd_pot",rtype,datafile);
    getvar_s(&Physics.bnd[i_bnd_bcrit], "bnd_bcrit",rtype,datafile);
    getvar(&Physics.bnd[i_bnd_eps_top], "bnd_eps_top",rtype,datafile);
    getvar(&Physics.bnd[i_bnd_fem],     "bnd_fem",rtype,datafile);

    getvar_s(&Physics.tvd_h,           "tvd_h","double[4]",datafile);
    getvar(&Physics.tvd_cs,            "tvd_cs","double[4]",datafile);
    getvar(&Physics.tvd[i_tvd_rholev], "tvd_rholev",rtype,datafile);
    getvar(&Physics.tvd[i_tvd_rholog], "tvd_rholog",rtype,datafile);
    getvar(&Physics.tvd[i_tvd_qrho],   "tvd_qrho",rtype,datafile);
    getvar(&Physics.tvd[i_tvd_Bpar],   "tvd_Bpar",rtype,datafile);
    getvar(&Physics.tvd[i_tvd_vhyp],   "tvd_vhyp",rtype,datafile);
    getvar(&Physics.tvd[i_tvd_Qdiff_bnd], "tvd_Qdiff_bnd",rtype,datafile);
    getvar(&Physics.tvd[i_tvd_pm_v],      "tvd_pm_v",rtype,datafile);
    getvar(&Physics.tvd[i_tvd_pm_B],      "tvd_pm_B",rtype,datafile);
    getvar(&Physics.tvd[i_tvd_vmax_lim],  "tvd_vmax_lim",rtype,datafile);
    getvar(&Physics.tvd[i_tvd_CME_thresh],"tvd_CME_thresh",rtype,datafile);

    getvar_s(&Physics.tvd_h_bnd,   "tvd_h_bnd",   "double[2]",datafile);
    getvar_s(&Physics.tvd_visc_bnd,"tvd_visc_bnd","double[2]",datafile);
    getvar_s(&Physics.tvd_eta_bnd, "tvd_eta_bnd", "double[2]",datafile);
    getvar(&Physics.tvd_coeff,   "tvd_coeff",   "double[4]",datafile);

    getvar_s(&Physics.divB[i_divB_switch],"divB_switch",rtype,datafile);
    getvar(&Physics.divB[i_divB_itmax], "divB_itmax",rtype,datafile);
    getvar(&Physics.divB[i_divB_err],   "divB_err",rtype,datafile);

    getvar_s(&Physics.tchk[i_tchk_eps_min], "tchk_eps_min",rtype,datafile);
    getvar_s(&Physics.tchk[i_tchk_rho_min], "tchk_rho_min",rtype,datafile);
    getvar_s(&Physics.tchk[i_tchk_eps_max], "tchk_eps_max",rtype,datafile);
    getvar_s(&Physics.tchk[i_tchk_vmax],    "tchk_vmax",rtype,datafile);

    getvar_s(&Physics.dmp[i_dmp_switch], "dmp_switch",rtype,datafile);
    getvar_s(&Physics.dmp[i_dmp_tau_ref],"dmp_tau_ref",rtype,datafile);
    getvar_s(&Physics.dmp[i_dmp_vel_ref],"dmp_vel_ref",rtype,datafile);
    getvar_s(&Physics.dmp[i_dmp_tau_min],"dmp_tau_min",rtype,datafile);

    getvar_s(&Physics.rt[i_rt_update],"rt_update",rtype,datafile);
    getvar(&Physics.rt[i_rt_tau_min],"rt_tau_min",rtype,datafile);
    getvar(&Physics.rt[i_rt_tr_tem],"rt_tr_tem",rtype,datafile);
    getvar(&Physics.rt[i_rt_tr_pre],"rt_tr_pre",rtype,datafile);
    getvar(&Physics.rt[i_rt_pre_cut],"rt_pre_cut",rtype,datafile);
    getvar_s(&Physics.rt[i_rt_tstep], "rt_tstep",rtype,datafile);
    getvar(&Physics.rt[i_rt_cfl],"rt_cfl",rtype,datafile);
    getvar_s(&Physics.rt[i_rt_type],"rt_type",rtype,datafile);
    getvar(&Physics.rt[i_rt_epsilon],"rt_epsilon",rtype,datafile);
    getvar(&Physics.rt[i_rt_iout],"rt_iout",rtype,datafile);
    
    getvar_s(&Physics.rt_ext[i_ext_chr],"ext_chr","int",datafile);
    getvar_s(&Physics.rt_ext[i_ext_cor],"ext_cor","int",datafile);
    getvar(&Physics.rt_ext[i_ext_hlines],"ext_hlines","int",datafile);
    getvar(&Physics.rt_ext[i_ext_mglines],"ext_mglines","int",datafile);
    getvar(&Physics.rt_ext[i_ext_calines],"ext_calines","int",datafile);

    // If set to off (0) or Chianti (1) then turn line losses off

    if (Physics.rt_ext[i_ext_cor]<2){
      Physics.rt_ext[i_ext_hlines]=0;
      Physics.rt_ext[i_ext_mglines]=0;
      Physics.rt_ext[i_ext_calines]=0;
    } else if (Physics.rt_ext[i_ext_cor]==2) {
      Physics.rt_ext[i_ext_hlines]=min(Physics.rt_ext[i_ext_hlines],1);
      Physics.rt_ext[i_ext_mglines]=min(Physics.rt_ext[i_ext_mglines],1);
      Physics.rt_ext[i_ext_calines]=min(Physics.rt_ext[i_ext_calines],1);
    }

    getvar(&Physics.slice[i_sl_collect],"sl_collect","int",datafile);
    getvar(&Physics.slice[i_sl_ic],     "sl_I_out","int",datafile); 

    int nsl;
    getvar(&nsl,"sl_tau","int",datafile); 
    if (nsl > 20) nsl = 20;
    Physics.slice[i_sl_tau] = nsl;
    getvar(&Physics.tau_lev,"tau_lev","double[20]",datafile);

    getvar(&nsl,"sl_xz","int",datafile); 
    if (nsl > 20) nsl = 20;

    Physics.slice[i_sl_xz] = nsl;    
    getvar(&Physics.xz_lev,"xz_lev","int[20]",datafile);

    getvar(&nsl,"sl_xy","int",datafile); 
    if (nsl > 20) nsl = 20;
    Physics.slice[i_sl_xy] = nsl;    
    getvar(&Physics.xy_lev,"xy_lev","int[20]",datafile);

    getvar(&nsl,"sl_yz","int",datafile); 
    if (nsl > 20) nsl = 20;
    Physics.slice[i_sl_yz] = nsl;    
    getvar(&Physics.yz_lev,"yz_lev","int[20]",datafile);

    getvar(&Physics.tau_var,"tau_var","int[14]",datafile); 
    getvar(&Physics.xz_var,"xz_var","int[13]",datafile); 
    getvar(&Physics.xy_var,"xy_var","int[12]",datafile); 
    getvar(&Physics.yz_var,"yz_var","int[12]",datafile);    
  }

  MPI_Bcast(&Run,sizeof(Run),MPI_BYTE,0,MPI_COMM_WORLD);

  MPI_Bcast(&Grid.NDIM,  1,MPI_INT, 0,MPI_COMM_WORLD);
  MPI_Bcast(Grid.gxmin,  3,REALTYPE,0,MPI_COMM_WORLD);
  MPI_Bcast(Grid.gxmax,  3,REALTYPE,0,MPI_COMM_WORLD);
  MPI_Bcast(Grid.gsize,  3,MPI_INT, 0,MPI_COMM_WORLD);
  MPI_Bcast(Grid.ghosts, 3,MPI_INT, 0,MPI_COMM_WORLD);
  MPI_Bcast(Grid.pardim, 3,MPI_INT, 0,MPI_COMM_WORLD);
  MPI_Bcast(Grid.periods,3,MPI_INT, 0,MPI_COMM_WORLD);
  MPI_Bcast(Grid.procs,3,MPI_INT, 0,MPI_COMM_WORLD);

  MPI_Bcast(&Physics,sizeof(Physics),MPI_BYTE,0,MPI_COMM_WORLD);
  
  int tot_procs=1;
  for (int dim = 0; dim<3;dim++)
    tot_procs*=Grid.procs[dim];

  int nprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  if (tot_procs != nprocs) {

    cout << "!!! ERROR procs set in parameters.dat != number of processors!!" << endl;
    MPI_Abort(MPI_COMM_WORLD,1);

  }

  Run.Init(rank);
  Physics.Init();
  Grid.Init(Run,Physics);
  comm_split_init(Run,Grid);

  if( rank==0 ) {
    /*
    if( HasBackupFile(Run.backfile) )
      mode = restart;
    else
      mode = newrun;
    */
    Run.Show();
    Grid.Show();
    Physics.Show();
  }
  

  //MPI_Bcast(&mode,1,MPI_INT,0,MPI_COMM_WORLD);

  int ret = ds_IO_Init(Grid, Run);

  //return mode;
  return ret;
}

int main(int argc, char* argv[]) {

    RunData Run;
    GridData Grid, GridOp;
    PhysicsData Physics;

    MPI_Init(&argc,&argv);

    if(Initialize(Run, Grid, Physics) != 0) {
        std::cout<<"Init failed. Aborting ... "<<std::endl;
        MPI_Abort(MPI_COMM_WORLD,1);
    }

    if(Run.dspaces_optimized) {
        if(Initialize(Run, GridOp, Physics) != 0) {
            std::cout<<"GridOp Init failed. Aborting ... "<<std::endl;
            MPI_Abort(MPI_COMM_WORLD,1);
        }
    }

    while(Run.IsValid()) {
        if (Run.rank ==0)
            cout << "***** Start new iteration ***********************" << endl;

        int needread = Run.NeedsOutput();

        if(needread) {
            eos_ds_read(Run,Grid,Physics);
            if(Run.dspaces_optimized) {
                eos_ds_read_optimized(Run,GridOp,Physics);
                int res = eos_compare(Grid, GridOp, Run, Physics);
            }
            /*
            if(Run.diagnostics) {
                diag_output();
            }
            */
        }
    }

    comm_split_finalize();
    ds_IO_Finalize();
    int ierr = MPI_Finalize();

    return 0;
}


