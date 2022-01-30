#include <dspaces.h>
#include <margo.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>

#include "dfparser.h"
#include "physics.H"
#include "grid.H"
#include "run.H"
#include "ds_common.H"
#include "ds_file_writer.H"

DSGridData::DSGridData() {
    ndim = 3;
    for(int i=0; i<3; i++) {
        gsize[i] = 1;
        periods[i] = 0;
        procs[i] = 1;
    }
}

void DSGridData::Init(MPI_Comm comm) {
    int reorder = 0;
    int remain[3];
    MPI_Comm_rank(gcomm, &grank);
    MPI_Cart_create(gcomm, ndim, procs, periods, reorder, &cart_comm);
    MPI_Cart_coords(cart_comm, grank, ndim, coord);

    for(int d=0; d<ndim; d++) {
        // Calculate size of local grid nx/ncores_x etc
        lsize[d] = (int) gsize[d]/procs[d];

        // My Starting point in the grid
        start[d] = lsize[d] * coord[d];

        // any remaining points
        remain[d] = gsize[d] % procs[d];

        // distribute remaining points
        if(coord[d] < remain[d]) {
            start[d] += coord[d];
            lsize[d] += 1;
        } else {
            start[d] += remain[d];
        }

        // My end point in the grid
        end[d] = start[d] + lsize[d] - 1;
    }

}

void DSGridData::Show() const {
    std::cout << " ------------DataSpaces_File_Writer Grid Parameter Settings -------------" << std::endl;
    std::cout << "Decomposition = " << procs[0] << 'x' << procs[1] << 'x' << procs[2] << std::endl
              << "gsize         = " << gsize[0] << ' ' << gsize[1] << ' ' << gsize[2] << std::endl
              << "lsize         = " << lsize[0] << ' ' << lsize[1] << ' ' << lsize[2] << std::endl
              << "start         = " << start[0] << ' ' << start[1] << ' ' << start[2] << std::endl
              << "end           = " << end[0] << ' ' << end[1] << ' ' << end[2] << std::endl;
}

void write_eos(dspaces_client_t client, const RunData& Run, const GridData& Grid,const PhysicsData& Physics,
                const DSGridData & DSGrid, int globiter, MPI_Comm comm)
{
    double clk, time_get = 0, time_mpi_file = 0;
    char ds_var_name[128];
    char filename[128];
    uint64_t lb[3], ub[3];
    int var;
    MPI_File mfh;
    int max_vars = 14;
    char eos_names[max_vars][128];

    //Only write out variables that re allocated based on physics configuration
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
    for(int v=0; v<max_vars; v++){
        if((Run.eos_output[v] == 1) && (var_init[v] == 1)) {
            var_index[tot_vars] = v;
            tot_vars += 1;
        }
    }

    // This is ugly, but for now it works so I will stick with it -- DP 
    sprintf(eos_names[0], "%s", "eosT");
    sprintf(eos_names[1], "%s", "eosP");
    sprintf(eos_names[2], "%s", "eosne");
    sprintf(eos_names[3], "%s", "eosrhoi");
    sprintf(eos_names[4], "%s", "eosamb");
    sprintf(eos_names[5], "%s", "Qtot");
    sprintf(eos_names[6], "%s", "tau");
    sprintf(eos_names[7], "%s", "Jtot");
    sprintf(eos_names[8], "%s", "Stot");
    sprintf(eos_names[9], "%s", "QxCor");
    sprintf(eos_names[10], "%s", "QxH");
    sprintf(eos_names[11], "%s", "QxMg");
    sprintf(eos_names[12], "%s", "QxCa");
    sprintf(eos_names[13], "%s", "QxChr");

    uint64_t vol = 1;
    for(int d=0; d<3; d++) {
        lb[d] = DSGrid.start[d];
        ub[d] = DSGrid.end[d];
        vol *= DSGrid.lsize[d];
    }

    void* buffer = (void*) malloc(vol*sizeof(float));

    for(int v=0; v<tot_vars; v++) {
        var = var_index[v];
        sprintf(ds_var_name, "%s%s", Run.path_3D, eos_names[var]);
        sprintf(filename,"%s%s.%06d",Run.path_3D, eos_names[var], globiter);

        clk = MPI_Wtime();
        dspaces_get(client, ds_var_name, globiter, sizeof(float), 3, lb, ub, buffer, -1);
        time_get += MPI_Wtime() - clk;

        MPI_Datatype io_subarray;
        MPI_Type_create_subarray(3, DSGrid.gsize, DSGrid.lsize, DSGrid.start, MPI_ORDER_FORTRAN, MPI_FLOAT, &io_subarray);
        MPI_Type_commit(&io_subarray);

        clk = MPI_Wtime();
        MPI_File_open(comm, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &mfh);
        MPI_File_set_view(mfh, 0, MPI_FLOAT, io_subarray, (char *) "native", MPI_INFO_NULL);
        MPI_File_write_all(mfh, buffer, vol, MPI_FLOAT, MPI_STATUS_IGNORE);
        MPI_File_close(&mfh);
        MPI_File_close(&mfh);
    }
    free(buffer);
    if(DSGrid.grank == 0) {
        fprintf(stdout, "Write EOS: GlobalIter = %d, Time of dspaces_get() = %lf, Time of MPI_File_write() = %lf.\n",
                globiter, time_get, time_mpi_file);
    }
}

void Initialize(RunData& Run,GridData& Grid, PhysicsData& Physics, DSGridData ds_Grid, MPI_Comm gcomm) {
    char datafile[256] = "parameters.dat";
    char rtype[16] = "double";
    int i, rank;
    if( sizeof(real)==sizeof(float) ) strcpy(rtype,"float");

    MPI_Comm_rank(gcomm, &rank);


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
        getvar(&Run.dspaces_put_local, "dspaces_put_local", "int", datafile);
        getvar(&Run.dspaces_gpu, "dspaces_gpu", "int", datafile);
        getvar(&Run.dspaces_manual_listen_addr, "dspaces_manual_listen_addr", "int", datafile);
        if(Run.dspaces_manual_listen_addr) {
            getvar(&Run.dspaces_client_listen_addr, "dspaces_client_listen_addr", "char*", datafile);
        }
        getvar(&Run.io_log_path, "io_log_path", "char*", datafile);
        getvar_s(&ds_Grid.ndim,"NDIM","int",datafile);
        getvar_s(&Grid.NDIM,"NDIM","int",datafile);

        if( Grid.NDIM==3 ) {
            strcat(ratype,"[3]");
            getvar(Grid.gxmin,"gxmin",ratype,datafile);
            getvar(Grid.gxmax,"gxmax",ratype,datafile);
            getvar(Grid.gsize,"gsize","int[3]",datafile);
            getvar(Grid.pardim,"pardim","int[3]",datafile);
            getvar(Grid.periods,"periods","int[3]",datafile);
            getvar(Grid.procs,"procs","int[3]",datafile);

            // Only works for 3D now
            getvar(ds_Grid.gsize,"gsize","int[3]",datafile);
            getvar(ds_Grid.procs, "dspaces_file_writer_procs", "int[3]", datafile);

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

    MPI_Bcast(ds_Grid.procs,3,MPI_INT, 0,MPI_COMM_WORLD);

    MPI_Bcast(&Physics,sizeof(Physics),MPI_BYTE,0,MPI_COMM_WORLD);

    ds_Grid.Init(gcomm);
    //Run.Init(rank);
    //Physics.Init();
    //Grid.Init(Run,Physics);
    //   comm_split_init(Run,Grid);

    if(rank == 0) {
        // Run.Show();
        // Grid.Show();
        // Physics.Show();
        ds_Grid.Show();
    }
  
}

int main(int argc, char** argv) {
    dspaces_client_t client;
    int rank, color;
    int ret;

    RunData Run;
    GridData Grid;
    PhysicsData Physics;
    DSGridData DSGrid;

    struct ds_meta *mdata;
    unsigned int mdatalen;
    int mver, mverlast;

    double clk;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm gcomm = MPI_COMM_WORLD;

    color = 1;
    MPI_Comm_split(MPI_COMM_WORLD, color, rank, &gcomm);

    ret = dspaces_init(rank, &client);
    if(ret != 0) {
        return ret;
    }

    fprintf(stdout, "DEBUG0\n");

    Initialize(Run, Grid, Physics, DSGrid, gcomm);

    // make sure MPI IO errors leed to program termination
    MPI_File_set_errhandler(MPI_FILE_NULL,MPI_ERRORS_ARE_FATAL); 

    if(rank != 0) {
        mdata = (struct ds_meta *)malloc(sizeof(*mdata));
    }
    mverlast = -1;
    do {

        if(rank == 0) {
            clk = MPI_Wtime();
            dspaces_get_meta(client, "muram_meta", META_MODE_NEXT, mverlast, &mver, (void **)&mdata, &mdatalen);
            fprintf(stdout, "Time of dspaces_get_meta = %lf, Meta Version = %d.\n", MPI_Wtime()-clk, mver);
            if(mdatalen != sizeof(*mdata)) {
                fprintf(stderr, "ERROR: corrupt metadata of size %d\n", mdatalen);
                free(mdata);
                mdata = (struct ds_meta *)malloc(sizeof(*mdata));
                mdata->globiter = STEP_ERR; 
            } else if(mdata->globiter == STEP_DONE) {
                fprintf(stdout, "MURaM is complete. Shutting down...\n");
            }
        }

        MPI_Bcast(mdata, sizeof(*mdata), MPI_BYTE, 0, gcomm);

        if(mdata->globiter == STEP_DONE || mdata->globiter == STEP_ERR) {
            break;
        }

        switch (mdata->var)
        {
        // 3D vars
        case EOS:
            if(rank == 0) {
                fprintf(stdout, "Rank: %d: Write EOS: GlobalIter = %d ...\n", rank, mdata->globiter);
            }
            write_eos(client, Run, Grid, Physics, DSGrid,mdata->globiter, gcomm);
            if(rank == 0) {
                fprintf(stdout, "Rank: %d: Write EOS Done...\n", rank);
            }
            break;
        case DIAG:
            // use nslvar as DIAG_flag only for DIAG
            // if(rank == 0) {
            //     fprintf(stdout, "Rank: %d: Write DIAG: GlobalIter = %d ...\n", rank, mdata->globiter);
            // }
            // write_diag(s, Run, Grid, Physics, mdata->globiter, gcomm, mdata->nslvar);
            // if(rank == 0) {
            //     fprintf(stdout, "Rank: %d: Write DIAG Done...\n", rank);
            // }
            break;
        case SOLUTION:
            // write_solution(s, Run, Grid, Physics, mdata->globiter, gcomm, metadata->nslvar);
            break;
        // 2D vars
        case TAU_SLICE:
            // write_tau_slice(s, Run, Grid, mdata->globiter, gcomm);
            break;
        case YZ_SLICE:
            // write_yz_slice(s, Run, Grid, mdata->globiter, gcomm);
            break;
        case XY_SLICE:
            // write_xy_slice(s, Run, Grid, mdata->globiter, gcomm);
            break;
        case XZ_SLICE:
            // write_xz_slice(s, Run, Grid, mdata->globiter, gcomm);
            break;
        case CORONA:
            // write_corona(s, Run, Grid, mdata->globiter, gcomm);
            break;
        // 1D vars
        case ANALYZE_VP:
            // write_analyze_vp(s, Grid, Run, mdata->globiter, gcomm);
            break;
        default:
            break;
        }

        mverlast = mver;

    } while (1);

    if(Run.dspaces_terminate && rank == 0) {
        dspaces_kill(client);
    }

    MPI_Barrier(gcomm);

    dspaces_fini(client);
    
    if(rank == 0) {
        fprintf(stdout, "File Writer is all done!\n");
    }

    MPI_Finalize();
    return 0;
}