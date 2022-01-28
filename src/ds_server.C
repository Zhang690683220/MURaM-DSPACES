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
#include <string>
#include <iostream>

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
//   comm_split_init(Run,Grid);

  if(rank == 0) {
    // Run.Show();
    // Grid.Show();
    // Physics.Show();
  }
  
}

void write_eos(dspaces_client_t client, const RunData& Run, const GridData& Grid,const PhysicsData& Physics, int globiter, MPI_Comm comm)
{
    double clk, time_find_objs = 0, time_get_objs = 0, time_mpi_file = 0, time_get = 0;
    char ds_var_name[128];
    char filename[128];
    int gsz[3], lsz[3], str[3];
    int lb[3], ub[3];
    int var;
    int io_rank, nprocs;
    MPI_File mfh;
    int max_vars = 14;

    MPI_Comm_size(comm, &nprocs);
    MPI_Comm_rank(comm, &io_rank);
  
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
            var_index[tot_vars]=v;
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
        gsz[d] = Grid.gsize[d];
    }

    lb[0] = (gsz[0] / nprocs)*io_rank;
    lb[1] = 0;
    lb[2] = 0;

    if(io_rank < (gsz[0] % nprocs)) {
        lb[0] += io_rank;
        ub[0] = lb[0] + (gsz[0] / nprocs);
    } else {
        lb[0] += gsz[0] % nprocs;
        ub[0] = lb[0] + (gsz[0] / nprocs) -1;
    }
    ub[1] = lb[1] + lsz[1];
    ub[2] = lb[2] + lsz[2];



    uint64_t vol = 1;
    for(int d=0; d<3; d++) {
        lsz = ub[d] - lb[d] + 1;
        vol = vol * lsz[d];
    }

    void* buffer = (void*) malloc(vol*sizeof(float));

    for(int v=0; v<tot_vars; v++) {
        var = var_index[v];
        sprintf(ds_var_name, "%s%s", Run.path_3D, eos_names[var]);
        sprintf(filename,"%s%s.%06d",Run.path_3D, eos_names[var], globiter);

        // struct dspaces_data_obj *objs;
        // clk = MPI_Wtime();
        // int obj_num = dspaces_server_find_objs(server, ds_var_name, globiter, &objs);
        // time_find_objs += MPI_Wtime() - clk;
        // fprintf(stdout, "Rank %d: dspaces_server_find_objs() Find %d objs.\n", io_rank, obj_num);

        clk = MPI_Wtime();
        dspaces_get(client, ds_var_name, globiter, sizeof(float), 3, lb, ub, buffer, -1);
        time_get += MPI_Wtime() - clk;

        int obj_num_max;
        MPI_Allreduce(&obj_num, &obj_num_max, 1, MPI_INT, MPI_MAX, comm);

        MPI_File_open(comm, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &mfh);

        // MPI_Datatype *io_subarray = (MPI_Datatype*) malloc(obj_num_max*sizeof(*io_subarray));

        MPI_Datatype io_subarray;

        clk = MPI_Wtime();
        MPI_Type_create_subarray(3, gsz, lsz, lb, MPI_ORDER_FORTRAN, MPI_FLOAT, &io_subarray);
        MPI_Type_commit(&io_subarray);
        MPI_File_set_view(mfh, 0, MPI_FLOAT, io_subarray, (char *) "native", MPI_INFO_NULL);
        MPI_File_write_all(mfh, NULL, 0, MPI_FLOAT, MPI_STATUS_IGNORE);
        time_mpi_file += MPI_Wtime() - clk;

        // for(int i=0; i<obj_num_max; i++) {
        //     if(i < obj_num) {
        //         std::string lszmsg, strmsg;
        //         lszmsg = "lsz = {";
        //         strmsg = "str = {";
        //         uint64_t vol = 1;
        //         int elem_size = objs[i].size;
        //         for(int d=0; d<objs[i].ndim; d++) {
        //             lsz[d] = objs[i].ub[d] - objs[i].lb[d] + 1;
        //             str[d] = objs[i].lb[d];
        //             vol = vol * lsz[d];
        //             lszmsg += " " + lsz[d];
        //             strmsg += " " + str[d];
        //         }

        //         lszmsg += "}";
        //         strmsg += "}";

        //         std::cout<< "Rank" << io_rank << ": EOS Write: Obj index = " << i << ", " << lszmsg << ", " << strmsg << std::endl; 

        //         void* buffer = (void*) malloc(elem_size*vol);

        //         clk = MPI_Wtime();
        //         int ret=dspaces_server_get_objdata(server, &objs[i], buffer);
        //         if(ret !=0) {
        //             fprintf(stdout, "Rank %d: EOS Write: dspaces_server_get_objdata() falied!"
        //                             "Total Objs = %d, Retriving Obj index = %d.\n", io_rank, obj_num, i);
        //         }
        //         time_get_objs += MPI_Wtime() - clk;

        //         clk = MPI_Wtime();
        //         MPI_Type_create_subarray(3, gsz, lsz, str, MPI_ORDER_FORTRAN, MPI_FLOAT, &io_subarray[i]);
        //         MPI_Type_commit(&io_subarray[i]);

	    //         MPI_File_set_view(mfh, 0, MPI_FLOAT, io_subarray[i], (char *) "native", MPI_INFO_NULL);

        //         // if(i < obj_num_min) {
        //             MPI_File_write_all(mfh, buffer, vol, MPI_FLOAT, MPI_STATUS_IGNORE);
        //         // } else {
	    //         //    MPI_File_write(mfh, buffer, vol, MPI_FLOAT, MPI_STATUS_IGNORE);
        //         // }
        //         MPI_Type_free(&io_subarray[i]);
        //         time_mpi_file += MPI_Wtime() - clk;
        //         free(buffer);
        //         // fprintf(stdout, "Rank %d: Total objs = %d, Write objs %d.\n", io_rank, obj_num,i);
        //     } else {
        //         // ! MPI_File_write_all() writes nothing here, Just to avoid collective stuck.
        //         // ! Fake code, doesn't do anything useful
        //         for(int d=0; d<3; d++) {
        //             lsz[d] = 1;
        //             str[d] = 0;
        //         }
        //         MPI_Type_create_subarray(3, gsz, lsz, str, MPI_ORDER_FORTRAN, MPI_FLOAT, &io_subarray[i]);
        //         MPI_Type_commit(&io_subarray[i]);
        //         MPI_File_set_view(mfh, 0, MPI_FLOAT, io_subarray[i], (char *) "native", MPI_INFO_NULL);

        //         MPI_File_write_all(mfh, NULL, 0, MPI_FLOAT, MPI_STATUS_IGNORE);

        //         MPI_Type_free(&io_subarray[i]);
        //     }
        // }

        // int bigbuf_size = 0;

        // for(int i=0; i<obj_num; i++) {
        //     uint64_t vol = 1;
        //     int elem_size = objs[i].size;
        //     for(int d=0; d<objs[i].ndim; d++) {
        //         lsz[d] = objs[i].ub[d] - objs[i].lb[d] + 1;
        //         str[d] = objs[i].lb[d];
        //         vol = vol * lsz[d];
        //     }
        //     bigbuf_size += elem_size*vol;
        // }

        // void *buffer = (void*) malloc(bigbuf_size);

        // int offset=0;
        // for(int i=0; i<obj_num; i++) {
        //     uint64_t vol = 1;
        //     int elem_size = objs[i].size;
        //     for(int d=0; d<objs[i].ndim; d++) {
        //         lsz[d] = objs[i].ub[d] - objs[i].lb[d] + 1;
        //         str[d] = objs[i].lb[d];
        //         vol = vol * lsz[d];
        //     }

        //     dspaces_server_get_objdata(server, &objs[i], buffer[offset]);

        //     offset = elem_size*vol;
        // }
        
        // int *procoff = (int*) malloc(nprocs*sizeof(int));
        // MPI_Allgather()

        // clk = MPI_Wtime();
        // MPI_Type_create_subarray(1, gsz, lsz, str, MPI_ORDER_FORTRAN, MPI_FLOAT, &io_subarray[i]);
        // MPI_Type_commit(&io_subarray[i]);

        // MPI_File_set_view(mfh, 0, MPI_FLOAT, io_subarray[i], (char *) "native", MPI_INFO_NULL);

        // // if(i < obj_num_min) {
        //     MPI_File_write_all(mfh, buffer, vol, MPI_FLOAT, MPI_STATUS_IGNORE);
        // // } else {
        // //    MPI_File_write(mfh, buffer, vol, MPI_FLOAT, MPI_STATUS_IGNORE);
        // // }
        // MPI_Type_free(&io_subarray[i]);
        // time_mpi_file += MPI_Wtime() - clk;
        // free(io_subarray);
        MPI_Barrier(comm);
        MPI_File_close(&mfh);

    }
    
    if(io_rank == 0) {
        // fprintf(stdout, "Time of find_objs() = %lf, Time of get_objs() = %lf, Time of MPI_File_write() = %lf.\n",
        //         time_find_objs, time_get_objs, time_mpi_file);
        fprintf(stdout, "Time of dspaces_get() = %lf, Time of MPI_File_write() = %lf.\n",
                time_get, time_mpi_file);
    }
}

// Physics is redundant in arguments
void write_diag(dspaces_provider_t server, const RunData& Run, const GridData& Grid,const PhysicsData& Physics, int globiter, MPI_Comm comm, int var_flag)
{
    double clk, time_find_objs = 0, time_get_objs = 0, time_mpi_file = 0;
    char ds_var_name[128];
    char filename[128];
    int gsz[3], lsz[3], str[3];
    int var;
    int io_rank;
    MPI_File mfh;
    int max_vars = 11;

    MPI_Comm_rank(comm, &io_rank);
  
    int var_index[max_vars];
    char diag_names[max_vars][128];

    int tot_vars = 0;
    // first 10 vars are determined  in param.dat
    for(int v=0; v<max_vars-1; v++) {
        if(Run.diag_output[v] == 1) {
	        var_index[tot_vars]=v;
	        tot_vars +=1;
        }
    }

    // get var_flag from meta sent by MURaM
    if(var_flag) {
        var_index[tot_vars] = max_vars-1;
        tot_vars +=1;
    }

    sprintf(diag_names[0],"%s","tvar1");
    sprintf(diag_names[1],"%s","tvar2");
    sprintf(diag_names[2],"%s","tvar3");
    sprintf(diag_names[3],"%s","tvar4");
    sprintf(diag_names[4],"%s","tvar5");
    sprintf(diag_names[5],"%s","tvar6");
    sprintf(diag_names[6],"%s","tvar7");
    sprintf(diag_names[7],"%s","tvar8");
    sprintf(diag_names[8],"%s","Qres");
    sprintf(diag_names[9],"%s","Qvis");
    sprintf(diag_names[10],"%s","Qamb");

    for(int d=0; d<3; d++) {
        gsz[d] = Grid.gsize[d];
    }

    for(int v=0; v<tot_vars; v++) {
        var = var_index[v];
        sprintf(ds_var_name, "%s%s", Run.path_3D, diag_names[var]);
        sprintf(filename,"%s%s.%06d",Run.path_3D, diag_names[var], globiter);

        struct dspaces_data_obj *objs;
        clk = MPI_Wtime();
        int obj_num = dspaces_server_find_objs(server, ds_var_name, globiter, &objs);
        time_find_objs += MPI_Wtime() - clk;

        int obj_num_max;
        MPI_Allreduce(&obj_num, &obj_num_max, 1, MPI_INT, MPI_MAX, comm);

        MPI_File_open(comm, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &mfh);

        MPI_Datatype *io_subarray = (MPI_Datatype*) malloc(obj_num_max*sizeof(*io_subarray));


        for(int i=0; i<obj_num_max; i++) {
            if(i < obj_num) {
                std::string lszmsg, strmsg;
                lszmsg = "lsz = {";
                strmsg = "str = {";
                uint64_t vol = 1;
                int elem_size = objs[i].size;
                for(int d=0; d<objs[i].ndim; d++) {
                    lsz[d] = objs[i].ub[d] - objs[i].lb[d] + 1;
                    str[d] = objs[i].lb[d];
                    vol = vol * lsz[d];
                    lszmsg += " " + lsz[d];
                    strmsg += " " + str[d];
                }

                lszmsg += "}";
                strmsg += "}";

                std::cout<< "Rank" << io_rank << ": EOS Write: Obj index = " << i << ", " << lszmsg << ", " << strmsg << std::endl; 
                void* buffer = (void*) malloc(elem_size*vol);

                clk = MPI_Wtime();
                int ret=dspaces_server_get_objdata(server, &objs[i], buffer);
                if(ret !=0) {
                    fprintf(stdout, "Rank %d: DIAG Write: dspaces_server_get_objdata() falied!"
                                    "Total Objs = %d, Retriving Obj index = %d.\n", io_rank, obj_num, i);
                }
                time_get_objs += MPI_Wtime() - clk;

                clk = MPI_Wtime();
                MPI_Type_create_subarray(3, gsz, lsz, str, MPI_ORDER_FORTRAN, MPI_FLOAT, &io_subarray[i]);
                MPI_Type_commit(&io_subarray[i]);

                MPI_File_set_view(mfh, 0, MPI_FLOAT, io_subarray[i], (char *) "native", MPI_INFO_NULL);
                MPI_File_write_all(mfh, buffer, vol, MPI_FLOAT, MPI_STATUS_IGNORE);
                
                MPI_Type_free(&io_subarray[i]);
                time_mpi_file += MPI_Wtime() - clk;
                free(buffer);
            } else {
                // ! MPI_File_write_all() writes nothing here, Just to avoid collective stuck.
                // ! Fake code, doesn't do anything useful
                for(int d=0; d<3; d++) {
                    lsz[d] = 1;
                    str[d] = 0;
                }
                MPI_Type_create_subarray(3, gsz, lsz, str, MPI_ORDER_FORTRAN, MPI_FLOAT, &io_subarray[i]);
                MPI_Type_commit(&io_subarray[i]);
                MPI_File_set_view(mfh, 0, MPI_FLOAT, io_subarray[i], (char *) "native", MPI_INFO_NULL);

                MPI_File_write_all(mfh, NULL, 0, MPI_FLOAT, MPI_STATUS_IGNORE);

                MPI_Type_free(&io_subarray[i]);
            }
        }

        free(io_subarray);
        MPI_File_close(&mfh);
    }

    
    if(io_rank == 0) {
        fprintf(stdout, "Time of find_objs() = %lf, Time of get_objs() = %lf, Time of MPI_File_write() = %lf.\n",
                time_find_objs, time_get_objs, time_mpi_file);
    }

}

void write_tau_slice(dspaces_provider_t server, const RunData& Run, const GridData& Grid,const PhysicsData& Physics, int globiter, MPI_Comm comm)
{
    // char ds_var_name[128];
    // char filename[128];
    // int gsz[2], lsz[2], str[2];
    // MPI_File mfh;

    // int nslice;
    // double* tau_lev;
    // int nslvar = 0;

    // nslice = Physics.slice[i_sl_tau];
    // tau_lev = (double*) malloc(nslice*sizeof(double));
    // for(int i=0; i<nslice; i++) {
    //     tau_lev[i] = Physics.tau_lev[i];
    // }

    // for(int v=0; v<14; v++) {
    //     if(Physics.tau_var[v] == 1) {
    //         nslvar+=1;
    //     }
    // }

    // gsz[0] = Grid.gsize[1];
    // gsz[1] = Grid.gsize[2];

    // for(int nsl=0; nsl<nslice; nsl++) {
    //     for()
    //     if(tau_lev[nsl] >= 1e-3) {
    //         sprintf(ds_var_name, "%s%s_%.3f_%d", Run.path_2D, "tau_slice", tau_lev[nsl]);
    //     } else {
	// 	    sprintf(ds_var_name, "%s%s_%.6f_%d", Run.path_2D, "tau_slice", tau_lev[nsl]);
    //     }
        
    //     // flag for whether collecting slices in different iters into 1 file
    //     if(Physics.slice[i_sl_collect] == 0) {
    //         if(tau_lev[nsl] >= 1e-3) {
    //             sprintf(filename, "%s%s_%.3f.%06d", Run.path_2D, "tau_slice", tau_lev[nsl], globiter);
    //         } else {
	// 	        sprintf(filename, "%s%s_%.6f.%06d", Run.path_2D, "tau_slice", tau_lev[nsl], globiter);
    //         }
    //         MPI_File_open(comm, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &mfh);
    //     } else {
    //         sprintf(filename, "%s_%.3f.dat", "tau_slice", tau_lev[nsl]);
    //         MPI_File_open(comm, filename, MPI_MODE_CREATE | MPI_MODE_APPEND, MPI_INFO_NULL, &mfh);
    //     }

    //     struct dspaces_data_obj *objs;
    //     int obj_num = dspaces_server_find_objs(server, ds_var_name, globiter, &objs);

    //     MPI_Datatype *io_subarray = (MPI_Datatype*) malloc(obj_num*sizeof(*io_subarray));

    //     for(int i=0; i<obj_num; i++) {
    //         uint64_t vol = 1;
    //         int elem_size = objs[i].size;
    //         for(int d=0; d<objs[i].ndim; d++) {
    //             lsz[d] = objs[i].ub[d] - objs[i].lb[d] + 1;
    //             str[d] = objs[i].lb[d];
    //             vol = vol * lsz[d];
    //         }
    //         void* buffer = (void*) malloc(elem_size*vol);
    //         dspaces_server_get_objdata(server, &objs[i], buffer);

    //         MPI_Type_create_subarray(2, gsz, lsz, str, MPI_ORDER_FORTRAN, MPI_FLOAT, &io_subarray[i]);
    //         MPI_Type_commit(&io_subarray[i]);

    //         MPI_File_set_view(mfh, 0, MPI_FLOAT, io_subarray[i], NULL, MPI_INFO_NULL);
	//         MPI_File_write(mfh, buffer, vol, MPI_FLOAT, MPI_STATUS_IGNORE);
	        
    //         MPI_Type_free(&io_subarray[i]);
    //         free(buffer);
    //     }

    //     free(io_subarray);
    //     MPI_File_close(&mfh);
    // }

}

void write_yz_slice(dspaces_provider_t server, const RunData& Run, const GridData& Grid,const PhysicsData& Physics, int globiter, MPI_Comm comm)
{
    // char ds_var_name[128];
    // char filename[128];
    // int gsz[2], lsz[2], str[2];
    // MPI_File mfh;

    // int nslice;
    // int *ixpos;
    // int nslvar = 0;

    // nslice = Physics.slice[i_sl_yz];
    // ixpos = (int*) malloc(nslice*sizeof(int));
    // for(int i=0; i<nslice; i++) {
    //     ixpos[i] = Physics.yz_lev[i];
    // }

    // for(int v=0; v<13; v++) {
    //     if(Physics.yz_var[v] == 1) {
    //         nslvar+=1;
    //     }
    // }

    // for(int nsl=0; nsl<nslice; nsl++) {
    //     sprintf(ds_var_name, "%s%s_%04d", Run.path_2D,"yz_slice",ixpos[nsl]);
    //     // flag for whether collecting slices in different iters into 1 file
    //     if(Physics.slice[i_sl_collect] == 0) {
    //         sprintf(filename, "%s%s_%04d.%06d", Run.path_2D, "yz_slice", ixpos[nsl], globiter);
    //         MPI_File_open(comm, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &mfh);
    //     } else {
    //         sprintf(filename,"%s_%04d.dat","yz_slice",ixpos[nsl]);
    //         MPI_File_open(comm, filename, MPI_MODE_CREATE | MPI_MODE_APPEND, MPI_INFO_NULL, &mfh);
    //     }

    //     struct dspaces_data_obj *objs;
    //     int obj_num = dspaces_server_find_objs(server, ds_var_name, globiter, &objs);

    //     MPI_Datatype *io_subarray = (MPI_Datatype*) malloc(obj_num*sizeof(*io_subarray));

    //     for(int i=0; i<obj_num; i++) {
    //         uint64_t vol = 1;
    //         int elem_size = objs[i].size;
    //         for(int d=0; d<objs[i].ndim; d++) {
    //             lsz[d] = objs[i].ub[d] - objs[i].lb[d] + 1;
    //             str[d] = objs[i].lb[d];
    //             vol = vol * lsz[d];
    //         }
    //         void* buffer = (void*) malloc(elem_size*vol);
    //         dspaces_server_get_objdata(server, &objs[i], buffer);

    //         MPI_Type_create_subarray(2, gsz, lsz, str, MPI_ORDER_FORTRAN, MPI_FLOAT, &io_subarray[i]);
    //         MPI_Type_commit(&io_subarray[i]);

    //         MPI_File_set_view(mfh, 0, MPI_FLOAT, io_subarray[i], NULL, MPI_INFO_NULL);
	//         MPI_File_write(mfh, buffer, vol, MPI_FLOAT, MPI_STATUS_IGNORE);
	        
    //         MPI_Type_free(&io_subarray[i]);
    //         free(buffer);

    //     }

    //     free(io_subarray);
    //     MPI_File_close(&mfh);
    // }

}

void write_xy_slice(dspaces_provider_t server, const RunData& Run, const GridData& Grid,const PhysicsData& Physics, int globiter, MPI_Comm comm)
{
    // char ds_var_name[128];
    // char filename[128];
    // int gsz[2], lsz[2], str[2];
    // MPI_File mfh;

    // int nslice;
    // int *ixpos;
    // int nslvar = 0;

    // nslice = Physics.slice[i_sl_xy];
    // ixpos = (int*) malloc(nslice*sizeof(int));
    // for(int i=0; i<nslice; i++) {
    //     ixpos[i] = Physics.xy_lev[i];
    // }

    // for(int v=0; v<12; v++) {
    //     if(Physics.xy_var[v] == 1) {
    //         nslvar+=1;
    //     }
    // }

    // for(int nsl=0; nsl<nslice; nsl++) {
    //     sprintf(ds_var_name, "%s%s_%04d", Run.path_2D,"xy_slice",ixpos[nsl]);
    //     // flag for whether collecting slices in different iters into 1 file
    //     if(Physics.slice[i_sl_collect] == 0) {
    //         sprintf(filename, "%s%s_%04d.%06d", Run.path_2D, "xy_slice",ixpos[nsl], globiter);
    //         MPI_File_open(comm, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &mfh);
    //     } else {
    //         sprintf(filename,"%s_%04d.dat","xy_slice",ixpos[nsl]);
    //         MPI_File_open(comm, filename, MPI_MODE_CREATE | MPI_MODE_APPEND, MPI_INFO_NULL, &mfh);
    //     }

    //     struct dspaces_data_obj *objs;
    //     int obj_num = dspaces_server_find_objs(server, ds_var_name, globiter, &objs);

    //     MPI_Datatype *io_subarray = (MPI_Datatype*) malloc(obj_num*sizeof(*io_subarray));

    //     for(int i=0; i<obj_num; i++) {
    //         uint64_t vol = 1;
    //         int elem_size = objs[i].size;
    //         for(int d=0; d<objs[i].ndim; d++) {
    //             lsz[d] = objs[i].ub[d] - objs[i].lb[d] + 1;
    //             str[d] = objs[i].lb[d];
    //             vol = vol * lsz[d];
    //         }
    //         void* buffer = (void*) malloc(elem_size*vol);
    //         dspaces_server_get_objdata(server, &objs[i], buffer);

    //         MPI_Type_create_subarray(2, gsz, lsz, str, MPI_ORDER_FORTRAN, MPI_FLOAT, &io_subarray[i]);
    //         MPI_Type_commit(&io_subarray[i]);

    //         MPI_File_set_view(mfh, 0, MPI_FLOAT, io_subarray[i], NULL, MPI_INFO_NULL);
	//         MPI_File_write(mfh, buffer, vol, MPI_FLOAT, MPI_STATUS_IGNORE);
	        
    //         MPI_Type_free(&io_subarray[i]);
    //         free(buffer);

    //     }

    //     free(io_subarray);
    //     MPI_File_close(&mfh);
    // }

}

void write_xz_slice(dspaces_provider_t server, const RunData& Run, const GridData& Grid,const PhysicsData& Physics, int globiter, MPI_Comm comm)
{
    // char ds_var_name[128];
    // char filename[128];
    // int gsz[2], lsz[2], str[2];
    // MPI_File mfh;

    // int nslice;
    // int *ixpos;
    // int nslvar = 0;

    // nslice = Physics.slice[i_sl_xz];
    // ixpos = (int*) malloc(nslice*sizeof(int));
    // for(int i=0; i<nslice; i++) {
    //     ixpos[i] = Physics.xz_lev[i];
    // }

    // for(int v=0; v<12; v++) {
    //     if(Physics.xz_var[v] == 1) {
    //         nslvar+=1;
    //     }
    // }

    // for(int nsl=0; nsl<nslice; nsl++) {
    //     // TODO: add nslvar to ds_var_name
    //     sprintf(ds_var_name, "%s%s_%04d", Run.path_2D,"xz_slice",ixpos[nsl]);
    //     // flag for whether collecting slices in different iters into 1 file
    //     if(Physics.slice[i_sl_collect] == 0) {
    //         sprintf(filename, "%s%s_%04d.%06d", Run.path_2D, "xz_slice", ixpos[nsl], globiter);
    //         MPI_File_open(comm, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &mfh);
    //     } else {
    //         sprintf(filename,"%s_%04d.dat","xz_slice",ixpos[nsl]);
    //         MPI_File_open(comm, filename, MPI_MODE_CREATE | MPI_MODE_APPEND, MPI_INFO_NULL, &mfh);
    //     }

    //     struct dspaces_data_obj *objs;
    //     int obj_num = dspaces_server_find_objs(server, ds_var_name, globiter, &objs);

    //     MPI_Datatype *io_subarray = (MPI_Datatype*) malloc(obj_num*sizeof(*io_subarray));

    //     for(int i=0; i<obj_num; i++) {
    //         uint64_t vol = 1;
    //         int elem_size = objs[i].size;
    //         for(int d=0; d<objs[i].ndim; d++) {
    //             lsz[d] = objs[i].ub[d] - objs[i].lb[d] + 1;
    //             str[d] = objs[i].lb[d];
    //             vol = vol * lsz[d];
    //         }
    //         void* buffer = (void*) malloc(elem_size*vol);
    //         dspaces_server_get_objdata(server, &objs[i], buffer);

    //         MPI_Type_create_subarray(2, gsz, lsz, str, MPI_ORDER_FORTRAN, MPI_FLOAT, &io_subarray[i]);
    //         MPI_Type_commit(&io_subarray[i]);

    //         MPI_File_set_view(mfh, 0, MPI_FLOAT, io_subarray[i], NULL, MPI_INFO_NULL);
	//         MPI_File_write(mfh, buffer, vol, MPI_FLOAT, MPI_STATUS_IGNORE);
	        
    //         MPI_Type_free(&io_subarray[i]);
    //         free(buffer);

    //     }

    //     free(io_subarray);
    //     MPI_File_close(&mfh);
    // }
}

void write_corona(dspaces_provider_t server, const RunData& Run, const GridData& Grid, int globiter, MPI_Comm comm)
{
    // char ds_var_name[128];
    // char filename[128];
    // const int nout = 4;
    // // YZ Slice

    // for(int v=0; v<nout; v++) {
    //     if(v == 0) {
	//         sprintf(ds_var_name, "%s%s", Run.path_2D, "corona_emission_adj_fil_x");
    //         sprintf(filename, "%s%s.%06d",Run.path_2D, "corona_emission_adj_fil_x", globiter);
    //     }
	//     if(v == 1) {
	//         sprintf(ds_var_name, "%s%s", Run.path_2D, "corona_emission_adj_dem_x");
    //         sprintf(filename, "%s%s.%06d", Run.path_2D, "corona_emission_adj_dem_x", globiter);
    //     }
	//     if(v == 2) {
	//         sprintf(ds_var_name, "%s%s",Run.path_2D, "corona_emission_adj_vlos_x");
    //         sprintf(filename, "%s%s.%06d", Run.path_2D, "corona_emission_adj_vlos_x", globiter);
    //     }
	//     if(v == 3) {
	//       sprintf(ds_var_name, "%s%s", Run.path_2D, "corona_emission_adj_vrms_x");
    //       sprintf(filename, "%s%s.%06d", Run.path_2D, "corona_emission_adj_vrms_x", globiter);
    //     }
    //     MPI_File_open(comm, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &mfh);

    //     struct dspaces_data_obj *objs;
    //     int obj_num = dspaces_server_find_objs(server, ds_var_name, globiter, &objs);

    //     MPI_Datatype *io_subarray = (MPI_Datatype*) malloc(obj_num*sizeof(*io_subarray));

    //     for(int i=0; i<obj_num; i++) {
    //         uint64_t vol = 1;
    //         int elem_size = objs[i].size;
    //         for(int d=0; d<objs[i].ndim; d++) {
    //             lsz[d] = objs[i].ub[d] - objs[i].lb[d] + 1;
    //             str[d] = objs[i].lb[d];
    //             vol = vol * lsz[d];
    //         }
    //         void* buffer = (void*) malloc(elem_size*vol);

    //         dspaces_server_get_objdata(server, &objs[i], buffer);

    //         // TODO: aggregate put should be reorganized

    //         MPI_Type_create_subarray(2, gsz, lsz, str, MPI_ORDER_FORTRAN, MPI_FLOAT, &io_subarray[i]);
    //         MPI_Type_commit(&io_subarray[i]);

    //         MPI_File_set_view(mfh, 0, MPI_FLOAT, io_subarray[i], NULL, MPI_INFO_NULL);
	//         MPI_File_write(mfh, buffer, vol, MPI_FLOAT, MPI_STATUS_IGNORE);
	        
    //         MPI_Type_free(&io_subarray[i]);
    //         free(buffer);
    //     }

    //     free(io_subarray);
    //     MPI_File_close(&mfh);
    // }

}

void write_analyze_vp(dspaces_provider_t server, const RunData& Run, const GridData& Grid, int globiter, MPI_Comm comm)
{
    // char ds_var_name[128];
    // char filename[128];
    // static int nvar = 50;

    // for(int v=0; v<nvar; v++) {
    //     sprintf(ds_var_name, "%s%s_%d", Run.path_2D, "hmean1D", v);
    //     sprintf(filename, "%s%s.%06d", Run.path_2D, "hmean1D", globiter);

    //     MPI_File_open(comm, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &mfh);

    //     // TODO: need a header

    //     struct dspaces_data_obj *objs;
    //     int obj_num = dspaces_server_find_objs(server, ds_var_name, globiter, &objs);

    //     MPI_Datatype *io_subarray = (MPI_Datatype*) malloc(obj_num*sizeof(*io_subarray));

    //     for(int i=0; i<obj_num; i++) {
    //         uint64_t vol = 1;
    //         int elem_size = objs[i].size;
    //         for(int d=0; d<objs[i].ndim; d++) {
    //             lsz[d] = objs[i].ub[d] - objs[i].lb[d] + 1;
    //             str[d] = objs[i].lb[d];
    //             vol = vol * lsz[d];
    //         }
    //         void* buffer = (void*) malloc(elem_size*vol);

    //         dspaces_server_get_objdata(server, &objs[i], buffer);
    //         MPI_Type_create_subarray(1, gsz, lsz, str, MPI_ORDER_FORTRAN, MPI_FLOAT, &io_subarray[i]);
    //         MPI_Type_commit(&io_subarray[i]);

    //         MPI_File_set_view(mfh, 4*sizeof(float), MPI_FLOAT, io_subarray[i], NULL, MPI_INFO_NULL);
    //         MPI_File_write(mfh, buffer, vol, MPI_FLOAT, MPI_STATUS_IGNORE);

    //         MPI_Type_free(&io_subarray[i]);
    //         free(buffer);
    //     }
    //     free(io_subarray);
    //     MPI_File_close(&mfh);
    // }
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
    double clk;

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
            clk = MPI_Wtime();
            dspaces_get_meta(dsp, "muram_meta", META_MODE_NEXT, mverlast, &mver, (void **)&mdata, &mdatalen);
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
            write_eos(s, Run, Grid, Physics, mdata->globiter, gcomm);
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

        // if(mdata->f_anl) {
        //     // Nothing writes at anlfreq yet.
        // }
        // if(mdata->f_res) {
        //     // TODO: need to wait the entire EOS domain to finish then start
        //     write_eos(s, Run, Grid, Physics, mdata->globiter, gcomm);
        //     // MURaM configures the proper diag var nums after 1st iter, so here we need to get metadata instead of reading param.dat
        //     // TODO: get tot_vars and var_index from get_meta()
        //     write_diag(s, Run, Grid, Physics, mdata->globiter, gcomm, var_flag);
        // }
        // if(mdata->f_slice) {
        //     if(Grid.NDIM > 1) {
        //        if(Physics.slice[i_sl_tau] > 0) {
        //            write_tau_slice(s, Run, Grid, mdata->globiter, gcomm);
        //        }
        //        if(Physics.slice[i_sl_yz] > 0) {
        //            write_yz_slice(s, Run, Grid, mdata->globiter, gcomm);
        //        }
        //     }
        //     if(Grid.NDIM == 3) {
        //         if(Physics.slice[i_sl_xy] > 0) {
        //             write_xy_slice(s, Run, Grid, mdata->globiter, gcomm);
        //         }
        //         if(Physics.slice[i_sl_xz] > 0) {
        //             write_xz_slice(s, Run, Grid, mdata->globiter, gcomm);
        //         }
        //     }
        //     if(Run.DEM){
        //         write_corona(s, Run, Grid, mdata->globiter, gcomm);
        //     }
        //     if(mdata->globiter > 0 && Run.HAVG) {
        //         write_analyze_vp(s, Grid, Run, mdata->globiter, gcomm);
        //     }
        // }
        // if(mdata->f_back) {
        //     write_solution(Run, Grid, mdata->globiter, gcomm);
        // }
    } while(1);

    MPI_Barrier(gcomm);

    dspaces_fini(dsp);

    // make margo wait for finalize
    dspaces_server_fini(s);

    if(rank == 0) {
        fprintf(stderr, "Server is all done!\n");
    }

    MPI_Finalize();
    return 0;
}
