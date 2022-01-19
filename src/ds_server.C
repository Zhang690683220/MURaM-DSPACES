/*
 * Copyright (c) 2020, Rutgers Discovery Informatics Institute, Rutgers
 * University
 *
 * See COPYRIGHT in top-level directory.
 */

#include <dspaces-server.h>
#include <dspaces.h>
#include <margo.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "comm_split.H"
#include "dfparser.h"
#include "physics.H"
#include "grid.H"
#include "run.H"
#include "ds_common.H"

void write_analyze_vp(dspaces_provider_t s, GridData Grid, RunData Run, int globiter, MPI_Comm gcomm)
{
}

void Initialize(RunData& Run,GridData& Grid, PhysicsData& Physics, MPI_Comm gcomm) {
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

  //Run.Init(rank);
  //Physics.Init();
  //Grid.Init(Run,Physics);
  comm_split_init(Run,Grid);

  if(rank == 0) {
    Run.Show();
    Grid.Show();
    Physics.Show();
  }
  
}

void write_eos(dspaces_provider_t server, const RunData& Run, const GridData& Grid,const PhysicsData& Physics, int globiter, MPI_Comm comm)
{
    char ds_var_name[128];
    char filename[128];
    int gsz[3], lsz[3], str[3];
    int var;
    MPI_File mfh;
    int max_vars = 14;
  
    char eos_names[max_vars][128];
    double* eos_vars[max_vars];

    //Only write out variables that re allocated based on phyasics configuration
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
    for(var=0;var<max_vars;var++){
        if((Run.eos_output[var] == 1) && (var_init[var] == 1)) {
            var_index[tot_vars]=var;
            tot_vars +=1;
        }
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

    for(int d=0; d<3; d++) {
        gsz[d] = Gird.gsize[d];
    }

    for(int v=0; v<tot_vars; v++) {
        var = var_index[v];
        sprintf(ds_var_name, "%s%s", Run.path_3D, eos_names[var]);
        sprintf(filename,"%s%s.%06d",Run.path_3D, eos_names[var], globiter);

        struct dspaces_data_obj **objs;
        int obj_num = dspaces_server_find_objs(server, ds_var_name, globiter, objs);

        MPI_File_open(comm, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &mfh);

        MPI_Datatype *io_subarray = (MPI_Datatype*) malloc(obj_num*sizeof(*io_subarray));


        for(int i=0; i<obj_num; i++) {
            uint64_t vol = 1;
            int elem_size = (*objs+i)->size;
            for(int d=0; d<(*objs+i)->ndim; d++) {
                lsz[d] = (*objs+i)->ub[d] - (*objs+i)->lb[d] + 1;
                str[d] = (*objs+i)->lb[d];
                vol = vol * lsz[d];
            }
            void* buffer = (void*) malloc(elem_size*vol);
            dspaces_server_get_objdata(server, (*objs+i), buffer);

            MPI_Type_create_subarray(3, gsz, lsz, str, MPI_ORDER_FORTRAN, MPI_FLOAT, &io_subarray[i]);
            MPI_Type_commit(&io_subarray[i]);

	        MPI_File_set_view(mfh, 0, MPI_FLOAT, io_subarray[i], NULL, MPI_INFO_NULL);
	        MPI_File_write_all(mfh, buffer, vol, MPI_FLOAT, MPI_STATUS_IGNORE);
	        
            MPI_Type_free(&io_subarray[i]);
            free(buffer);
        }

        free(io_subarray);
        MPI_File_close(&mfh);
    }
}

int main(int argc, char **argv)
{
    char *listen_addr_str;
    const char *conf_file;
    int rank, color;
    int mver, mverlast;
    dspaces_client_t dsp;
    dspaces_provider_t s; 
    struct ds_meta *mdata;
    unsigned int mdatalen;
    RunData Run;
    GridData Grid;
    PhysicsData Physics;
    int ret;

    if(argc < 2 || argc > 3) {
        fprintf(stderr, "Usage: %s <listen-address> [<conffile>]\n", argv[0]);
        return -1;
    }

    listen_addr_str = argv[1];
    if(argc == 3) {
        conf_file = argv[2];
    } else {
        conf_file = "dataspaces.conf";
    }

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm gcomm = MPI_COMM_WORLD;

    color = 1;
    MPI_Comm_split(MPI_COMM_WORLD, color, rank, &gcomm);

    ret = dspaces_server_init(listen_addr_str, gcomm, conf_file, &s);
    if(ret != 0)
        return ret;

    ret = dspaces_init(rank, &dsp);
    if(ret != 0)
        return ret;

    if(rank != 0) {
        mdata = (struct ds_meta *)malloc(sizeof(*mdata));
    }


    Initialize(Run, Grid, Physics, gcomm);

    // make sure MPI IO errors leed to program termination
    MPI_File_set_errhandler(MPI_FILE_NULL,MPI_ERRORS_ARE_FATAL); 

    mverlast = -1;
    do {
        
        if(rank == 0) {
            dspaces_get_meta(dsp, "muram_meta", META_MODE_NEXT, mverlast, &mver, (void **)&mdata, &mdatalen);
            if(mdatalen != sizeof(*mdata)) {
                fprintf(stderr, "ERROR: corrupt metadata of size %d\n", mdatalen);
                free(mdata);
                mdata = (struct ds_meta *)malloc(sizeof(*mdata));
                mdata->globiter = STEP_ERR; 
            } else if(mdata->globiter == STEP_DONE) {
                fprintf(stdout, "MURaM is complete. Shutting down...\n");
            }
        }

        MPI_Bcast(mdata, sizeof(mdata), MPI_BYTE, 0, gcomm);

        if(mdata->globiter == STEP_DONE || mdata->globiter == STEP_ERR) {
            break;
        }

        if(mdata->f_anl) {
            // Nothing writes at anlfreq yet.
        }
        if(mdata->f_res) {
            // TODO: need to wait the entire EOS domain to finish then start
            write_eos(Run, Grid, Physics, mdata->globiter, gcomm);
            write_diag(Run, Grid, Physics, mdata->globiter, gcomm);
        }
        if(mdata->f_slice) {
            if(Grid.NDIM > 1) {
               if(Physics.slice[i_sl_tau] > 0) {
                   write_tau_slice(Run, Grid, mdata->globiter, gcomm);
               }
               if(Physics.slice[i_sl_yz] > 0) {
                   write_yz_slice(Run, Grid, mdata->globiter, gcomm);
               }
            }
            if(Grid.NDIM == 3) {
                if(Physics.slice[i_sl_xy] > 0) {
                    write_yz_slice(Run, Grid, mdata->globiter, gcomm);
                }
                if(Physics.slice[i_sl_xz] > 0) {
                    write_xz_slice(Run, Grid, mdata->globiter, gcomm);
                }
            }
            if(Run.DEM){
                write_corona(Run, Grid, mdata->globiter, gcomm);
            }
            if(mdata->globiter > 0 && Run.HAVG) {
                write_analyze_vp(s, Grid, Run, mdata->globiter, gcomm);
            }
        }
        if(mdata->f_back) {
            write_solution(Run, Grid, mdata->globiter, gcomm);
        }
    } while(1);

    MPI_Barrier(gcomm);
    if(rank == 0) {
        dspaces_kill(dsp);
    }

    dspaces_fini(dsp);

    // make margo wait for finalize
    dspaces_server_fini(s);

    if(rank == 0) {
        fprintf(stderr, "Server is all done!\n");
    }

    MPI_Finalize();
    return 0;
}
