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
#include <netcdf.h>
#include <netcdf_par.h>

DSGridData::DSGridData() {
    ndim = 3;
    for(int i=0; i<3; i++) {
        gsize[i] = 1;
        periods[i] = 0;
        procs[i] = 1;
    }

    for(int i=0; i<2; i++) {
        procs2d[i] = 1;
        xygsize[i] = 1;
        xzgsize[i] = 1;
        yzgsize[i] = 1;
    }
}

void DSGridData::Init(MPI_Comm comm) {
    gcomm = comm;
    int reorder = 0;
    int remain[3];
    MPI_Comm_rank(gcomm, &grank);

    /* Create 3D domain decomposition */
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

    /* Create 2D domain decomposition */
    MPI_Cart_create(gcomm, 2, procs2d, periods, reorder, &comm2d);
    MPI_Cart_coords(comm2d, grank, 2, coord2d);

    xygsize[0] = gsize[0];
    xygsize[1] = gsize[1];

    xzgsize[0] = gsize[0];
    xzgsize[1] = gsize[2];

    yzgsize[0] = gsize[1];
    yzgsize[1] = gsize[2];

    int xyremain[2], xzremain[2], yzremain[2];

    for(int d=0; d<2; d++) {
        /* XY Plane */
        xylsize[d] = (int) xygsize[d]/procs2d[d];
        xystart[d] = xylsize[d] * coord2d[d];
        xyremain[d] = xygsize[d] % procs2d[d];

        if(coord2d[d] < xyremain[d]) {
            xystart[d] += coord2d[d];
            xylsize[d] += 1;
        } else {
            xystart[d] += xyremain[d];
        }

        xyend[d] = xystart[d] + xylsize[d] - 1;

        /* XZ Plane */
        xzlsize[d] = (int) xzgsize[d]/procs2d[d];
        xzstart[d] = xzlsize[d] * coord2d[d];
        xzremain[d] = xzgsize[d] % procs2d[d];

        if(coord2d[d] < xzremain[d]) {
            xzstart[d] += coord2d[d];
            xzlsize[d] += 1;
        } else {
            xzstart[d] += xzremain[d];
        }

        xzend[d] = xzstart[d] + xzlsize[d] - 1;

        /* YZ Plane */
        yzlsize[d] = (int) yzgsize[d]/procs2d[d];
        yzstart[d] = yzlsize[d] * coord2d[d];
        yzremain[d] = yzgsize[d] % procs2d[d];

        if(coord2d[d] < yzremain[d]) {
            yzstart[d] += coord2d[d];
            yzlsize[d] += 1;
        } else {
            yzstart[d] += yzremain[d];
        }

        yzend[d] = yzstart[d] + yzlsize[d] - 1;
    }

}

void DSGridData::Show() const {
    std::cout << " ------------DataSpaces_File_Writer Grid Parameter Settings -------------" << std::endl;
    std::cout << " ----------------------------------3D------------------------------------" << std::endl;
    std::cout << "Decomposition = " << procs[0] << 'x' << procs[1] << 'x' << procs[2] << std::endl
              << "gsize         = " << gsize[0] << ' ' << gsize[1] << ' ' << gsize[2] << std::endl
              << "lsize         = " << lsize[0] << ' ' << lsize[1] << ' ' << lsize[2] << std::endl
              << "start         = " << start[0] << ' ' << start[1] << ' ' << start[2] << std::endl
              << "end           = " << end[0] << ' ' << end[1] << ' ' << end[2] << std::endl;
    std::cout << " ----------------------------------2D------------------------------------" << std::endl;
    std::cout << "Decomposition = " << procs2d[0] << 'x' << procs2d[1] << std::endl;
    std::cout << " -------------------------------XY Plane---------------------------------" << std::endl;
    std::cout << "gsize         = " << xygsize[0] << ' ' << xygsize[1] << std::endl
              << "lsize         = " << xylsize[0] << ' ' << xylsize[1] << std::endl
              << "start         = " << xystart[0] << ' ' << xystart[1] << std::endl
              << "end           = " << xyend[0] << ' ' << xyend[1] << std::endl;
    std::cout << " -------------------------------XZ Plane---------------------------------" << std::endl;
    std::cout << "gsize         = " << xzgsize[0] << ' ' << xzgsize[1] << std::endl
              << "lsize         = " << xzlsize[0] << ' ' << xzlsize[1] << std::endl
              << "start         = " << xzstart[0] << ' ' << xzstart[1] << std::endl
              << "end           = " << xzend[0] << ' ' << xzend[1] << std::endl;
    std::cout << " -------------------------------YZ Plane---------------------------------" << std::endl;
    std::cout << "gsize         = " << yzgsize[0] << ' ' << yzgsize[1] << std::endl
              << "lsize         = " << yzlsize[0] << ' ' << yzlsize[1] << std::endl
              << "start         = " << yzstart[0] << ' ' << yzstart[1] << std::endl
              << "end           = " << yzend[0] << ' ' << yzend[1] << std::endl;
    std::cout << " ------------------------------------------------------------------------" << std::endl;
}

void write_eos(dspaces_client_t client, const RunData& Run, const PhysicsData& Physics,
                const DSGridData & DSGrid, int globiter)
{
    static int f_gdim=1;

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

    MPI_Datatype io_subarray;
    MPI_Type_create_subarray(3, DSGrid.gsize, DSGrid.lsize, DSGrid.start, MPI_ORDER_FORTRAN, MPI_FLOAT, &io_subarray);
    MPI_Type_commit(&io_subarray);

    void* buffer = (void*) malloc(vol*sizeof(float));

    for(int v=0; v<tot_vars; v++) {
        var = var_index[v];
        sprintf(ds_var_name, "%s%s", Run.path_3D, eos_names[var]);
        sprintf(filename,"%s%s.%06d",Run.path_3D, eos_names[var], globiter);

        if(f_gdim) {
            uint64_t gdim[3];
            for(int d=0; d<3; d++) {
                gdim[d] = DSGrid.gsize[d];
            }
            dspaces_define_gdim(client, ds_var_name, 3, gdim);
        }

        clk = MPI_Wtime();
        dspaces_get(client, ds_var_name, globiter, sizeof(float), 3, lb, ub, buffer, -1);
        time_get += MPI_Wtime() - clk;

        clk = MPI_Wtime();
        MPI_File_open(DSGrid.gcomm, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &mfh);
        MPI_File_set_view(mfh, 0, MPI_FLOAT, io_subarray, (char *) "native", MPI_INFO_NULL);
        MPI_File_write_all(mfh, buffer, vol, MPI_FLOAT, MPI_STATUS_IGNORE);
        MPI_File_close(&mfh);
        time_mpi_file += MPI_Wtime() - clk;
    }
    f_gdim = 0;
    free(buffer);
    if(DSGrid.grank == 0) {
        fprintf(stdout, "Write EOS: GlobalIter = %d, Time of dspaces_get() = %lf, Time of MPI_File_write() = %lf.\n",
                globiter, time_get, time_mpi_file);
    }
}

void nc_write_eos(dspaces_client_t client, const RunData& Run, const PhysicsData& Physics,
                    const DSGridData & DSGrid, int globiter)
{
    static int f_gdim = 1;

    double clk, time_get = 0, time_nc_file = 0;
    char ds_var_name[128];
    char filename[128];
    uint64_t lb[3], ub[3];
    int var;
    int max_vars = 14;
    char eos_names[max_vars][128];
    char nc_varname[128];
    int nc_fid;
    int nc_ret;
    int nc_dimid[3];
    char nc_dimname[3][128];
    size_t nc_str[3], nc_lsize[3];

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

    int* nc_varid = (int*) malloc(tot_vars*sizeof(int));

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
        nc_str[d] = DSGrid.start[d];
        nc_lsize[d] = DSGrid.lsize[d];
        vol *= DSGrid.lsize[d];
    }


    /* create parallel NetCDF file */
    sprintf(filename, "%snc_eos.%06d", Run.path_3D, globiter);
    nc_ret = nc_create_par(filename, NC_CLOBBER | NC_NETCDF4, DSGrid.gcomm, MPI_INFO_NULL, &nc_fid);
    if(nc_ret != NC_NOERR) {
        fprintf(stderr, "ERROR: Rank %i: %s, line %i (%s): nc_create_par() failed ! Error code: %s\n",
                DSGrid.grank, __FILE__, __LINE__, __func__, nc_strerror(nc_ret));
    }

    /* Define the dimensions. */
    sprintf(nc_dimname[0], "x");
    sprintf(nc_dimname[1], "y");
    sprintf(nc_dimname[2], "z");
    for(int d=0; d<3; d++) {
        nc_ret = nc_def_dim(nc_fid, nc_dimname[d], DSGrid.gsize[d], &nc_dimid[d]);
        if(nc_ret != NC_NOERR) {
            fprintf(stderr, "ERROR: Rank %i: %s, line %i (%s): nc_def_dim() failed ! Error code: %s\n",
                    DSGrid.grank, __FILE__, __LINE__, __func__, nc_strerror(nc_ret));
        }
    }

    float* buffer = (float*) malloc(vol*sizeof(float));

    for(int v=0; v<tot_vars; v++) {
        var = var_index[v];
        sprintf(ds_var_name, "%s%s", Run.path_3D, eos_names[var]);

        if(f_gdim) {
            uint64_t gdim[3];
            for(int d=0; d<3; d++) {
                gdim[d] = DSGrid.gsize[d];
            }
            dspaces_define_gdim(client, ds_var_name, 3, gdim);
        }

        /* Define the netCDF variables. The dimids array is used to pass
            the dimids of the dimensions of the variables.*/
        nc_ret = nc_def_var(nc_fid, eos_names[var], NC_FLOAT, 3, nc_dimid, &nc_varid[v]);
        if(nc_ret != NC_NOERR) {
            fprintf(stderr, "ERROR: Rank %i: %s, line %i (%s): nc_def_var() failed ! Error code: %s\n",
                    DSGrid.grank, __FILE__, __LINE__, __func__, nc_strerror(nc_ret));
        }

        /* Define units attributes for vars. */
        // char nc_var_unit[128];
        // nc_ret = nc_put_att_text(nc_fid, nc_varid[v], "units", strlen(nc_var_unit), nc_var_unit)

        /* End define mode. */
        // TODO: maybe cannot be used here
        // nc_ret = nc_enddef(nc_fid);

        clk = MPI_Wtime();
        dspaces_get(client, ds_var_name, globiter, sizeof(float), 3, lb, ub, (void*) buffer, -1);
        time_get += MPI_Wtime() - clk;

        /* Write Data */
        clk = MPI_Wtime();
        nc_ret = nc_put_vara_float(nc_fid, nc_varid[v], nc_str, nc_lsize, buffer);   
        time_nc_file += MPI_Wtime() - clk;
        if(nc_ret != NC_NOERR) {
            fprintf(stderr, "ERROR: Rank %i: %s, line %i (%s): nc_put_vara_float() failed ! Error code: %s\n",
                    DSGrid.grank, __FILE__, __LINE__, __func__, nc_strerror(nc_ret));
        }
    }
    f_gdim = 0;
    free(buffer);
    free(nc_varid);
    nc_ret = nc_close(nc_fid);
    if(nc_ret != NC_NOERR) {
            fprintf(stderr, "ERROR: Rank %i: %s, line %i (%s): nc_close() failed ! Error code: %s\n",
                    DSGrid.grank, __FILE__, __LINE__, __func__, nc_strerror(nc_ret));
    }

    if(DSGrid.grank == 0) {
        fprintf(stdout, "Write EOS: GlobalIter = %d, Time of dspaces_get() = %lf, Time of NetCDF_File_write() = %lf.\n",
                globiter, time_get, time_nc_file);
    }
}

void write_yz_slice(dspaces_client_t client, const RunData& Run, const PhysicsData& Physics,
                    const DSGridData & DSGrid, int globiter)
{
    static int f_gdim = 1;

    double clk, time_get = 0, time_mpi_file = 0;
    char ds_var_name[128];
    char filename[128];
    uint64_t lb[2], ub[2];

    int var;
    int max_vars = 13;
    int var_index[max_vars];

    int nslice = Physics.slice[i_sl_yz];
    int *ixpos = (int*) malloc(nslice*sizeof(int));
    int nslvar = 0;

    MPI_File mfh;

    for(int i=0; i<nslice; i++) {
        ixpos[i] = Physics.yz_lev[i];
    }

    for (int v=0; v<max_vars; v++) {
        if(Physics.yz_var[v] == 1) {
            var_index[nslvar] = v;
            nslvar += 1;
        }
    }

    uint64_t vol = 1;
    for(int d=0; d<2; d++) {
        lb[d] = DSGrid.yzstart[d];
        ub[d] = DSGrid.yzend[d];

        vol *= DSGrid.yzlsize[d];
    }

    MPI_Datatype io_subarray;
    MPI_Type_create_subarray(2, DSGrid.yzgsize, DSGrid.yzlsize, DSGrid.yzstart, MPI_ORDER_FORTRAN, MPI_FLOAT, &io_subarray);
    MPI_Type_commit(&io_subarray);

    float* buffer = (float*) malloc(vol*sizeof(float));

    for(int nsl=0; nsl<nslice; nsl++) {
        if(Physics.slice[i_sl_collect] == 0) {
            sprintf(filename, "%s%s_%04d.%06d", Run.path_2D, "yz_slice", ixpos[nsl], globiter);
            MPI_File_open(DSGrid.gcomm, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &mfh);

            for(int v=0; v<nslvar; v++) {
                var = var_index[v];
                sprintf(ds_var_name, "%s%s_%04d_%d", Run.path_2D, "yz_slice", ixpos[nsl], v);

                if(f_gdim) {
                    uint64_t gdim[2];
                    for(int d=0; d<2; d++) {
                        gdim[d] = DSGrid.yzgsize[d];
                    }
                    dspaces_define_gdim(client, ds_var_name, 2, gdim);
                }

                clk = MPI_Wtime();
                dspaces_get(client, ds_var_name, globiter, sizeof(float), 2, lb, ub, (void*) buffer, -1);
                time_get += MPI_Wtime() - clk;

                clk = MPI_Wtime();
                MPI_File_set_view(mfh, 0, MPI_FLOAT, io_subarray, (char *) "native", MPI_INFO_NULL);
                MPI_File_write_all(mfh, buffer, vol, MPI_FLOAT, MPI_STATUS_IGNORE);
                time_mpi_file += MPI_Wtime() - clk;
            }
        } else {
            // collective means write all iterations into the same file
            sprintf(filename,"%s_%04d.dat","yz_slice",ixpos[nsl]);
            MPI_File_open(DSGrid.gcomm, filename, MPI_MODE_APPEND | MPI_MODE_WRONLY, MPI_INFO_NULL, &mfh);

            for(int v=0; v<nslvar; v++) {
                var = var_index[v];
                sprintf(ds_var_name, "%s%s_%04d_%d", Run.path_2D, "yz_slice", ixpos[nsl], v);

                if(f_gdim) {
                    uint64_t gdim[2];
                    for(int d=0; d<2; d++) {
                        gdim[d] = DSGrid.yzgsize[d];
                    }
                    dspaces_define_gdim(client, ds_var_name, 2, gdim);
                }

                clk = MPI_Wtime();
                dspaces_get(client, ds_var_name, globiter, sizeof(float), 2, lb, ub, (void*) buffer, -1);
                time_get += MPI_Wtime() - clk;

                clk = MPI_Wtime();
                MPI_File_set_view(mfh, 0, MPI_FLOAT, io_subarray, (char *) "native", MPI_INFO_NULL);
                MPI_File_write_all(mfh, buffer, vol, MPI_FLOAT, MPI_STATUS_IGNORE);
                time_mpi_file += MPI_Wtime() - clk;

            }
        }
        MPI_File_close(&mfh);
    }
    f_gdim = 0;
    free(buffer);
    free(ixpos);

    if(DSGrid.grank == 0) {
        fprintf(stdout, "Write YZ_Slice: GlobalIter = %d, Time of dspaces_get() = %lf, Time of MPI_File_write() = %lf.\n",
                globiter, time_get, time_mpi_file);
    }
}

void nc_write_yz_slice(dspaces_client_t client, const RunData& Run, const PhysicsData& Physics,
                        const DSGridData & DSGrid, int globiter)
{
    static int f_gdim = 1;

    double clk, time_get = 0, time_nc_file = 0;
    char ds_var_name[128];
    char filename[128];
    uint64_t lb[2], ub[2];

    int var;
    int max_vars = 13;
    char yz_slice_names[max_vars][128];
    int var_index[max_vars];

    int nc_fid;
    int nc_ret;
    char nc_varname[128];
    int nc_dimid[2];
    char nc_dimname[2][128];
    size_t nc_str[2], nc_lsize[2];

    int nslice = Physics.slice[i_sl_yz];
    int *ixpos = (int*) malloc(nslice*sizeof(int));
    int nslvar = 0;

    for(int i=0; i<nslice; i++) {
        ixpos[i] = Physics.yz_lev[i];
    }

    for (int v=0; v<max_vars; v++) {
        sprintf(yz_slice_names[v], "var%d", v);
        if(Physics.yz_var[v] == 1) {
            var_index[nslvar] = v;
            nslvar += 1;
        }
    }

    int* nc_varid = (int*) malloc(nslvar*sizeof(int));

    uint64_t vol = 1;
    for(int d=0; d<2; d++) {
        lb[d] = DSGrid.yzstart[d];
        ub[d] = DSGrid.yzend[d];
        nc_str[d] = DSGrid.yzstart[d];
        nc_lsize[d] = DSGrid.yzlsize[d];
        vol *= DSGrid.yzlsize[d];
    }

    float* buffer = (float*) malloc(vol*sizeof(float));

    /* Define the dimensions. */
    sprintf(nc_dimname[0], "y");
    sprintf(nc_dimname[1], "z");


    for(int nsl=0; nsl<nslice; nsl++) {

        if(Physics.slice[i_sl_collect] == 0) {
            /* create parallel NetCDF file */
            sprintf(filename, "%snc_yz_slice_%04d.%06d", Run.path_2D, ixpos[nsl], globiter);
            nc_ret = nc_create_par(filename, NC_CLOBBER | NC_NETCDF4, DSGrid.gcomm, MPI_INFO_NULL, &nc_fid);
            if(nc_ret != NC_NOERR) {
                fprintf(stderr, "ERROR: Rank %i: %s, line %i (%s): nc_create_par() failed ! Error code: %s\n",
                        DSGrid.grank, __FILE__, __LINE__, __func__, nc_strerror(nc_ret));
            }

            /* Define the dimensions. */
            for(int d=0; d<2; d++) {
                nc_ret = nc_def_dim(nc_fid, nc_dimname[d], DSGrid.yzgsize[d], &nc_dimid[d]);
                if(nc_ret != NC_NOERR) {
                    fprintf(stderr, "ERROR: Rank %i: %s, line %i (%s): nc_def_dim() failed ! Error code: %s\n",
                            DSGrid.grank, __FILE__, __LINE__, __func__, nc_strerror(nc_ret));
                }
            }

            for(int v=0; v<nslvar; v++) {
                var = var_index[v];
                sprintf(ds_var_name, "%s%s_%04d_%d", Run.path_2D, "yz_slice", ixpos[nsl], v);

                if(f_gdim) {
                    uint64_t gdim[2];
                    for(int d=0; d<2; d++) {
                        gdim[d] = DSGrid.yzgsize[d];
                    }
                    dspaces_define_gdim(client, ds_var_name, 2, gdim);
                }

                /* Define the netCDF variables. The dimids array is used to pass
                    the dimids of the dimensions of the variables.*/
                nc_ret = nc_def_var(nc_fid, yz_slice_names[var], NC_FLOAT, 2, nc_dimid, &nc_varid[v]);
                if(nc_ret != NC_NOERR) {
                    fprintf(stderr, "ERROR: Rank %i: %s, line %i (%s): nc_def_var() failed ! Error code: %s\n",
                            DSGrid.grank, __FILE__, __LINE__, __func__, nc_strerror(nc_ret));
                }

                clk = MPI_Wtime();
                dspaces_get(client, ds_var_name, globiter, sizeof(float), 2, lb, ub, (void*) buffer, -1);
                time_get += MPI_Wtime() - clk;

                /* Write Data */
                clk = MPI_Wtime();
                nc_ret = nc_put_vara_float(nc_fid, nc_varid[v], nc_str, nc_lsize, buffer);   
                time_nc_file += MPI_Wtime() - clk;
                if(nc_ret != NC_NOERR) {
                    fprintf(stderr, "ERROR: Rank %i: %s, line %i (%s): nc_put_vara_float() failed ! Error code: %s\n",
                            DSGrid.grank, __FILE__, __LINE__, __func__, nc_strerror(nc_ret));
                }
            }
        
        } else {
            // collective means write all iterations into the same file
            char nc_var_name[128];
            /* create parallel NetCDF file */
            sprintf(filename, "yz_slice_%04d.dat", ixpos[nsl]);
            nc_ret = nc_create_par(filename, NC_NOCLOBBER | NC_NETCDF4, DSGrid.gcomm, MPI_INFO_NULL, &nc_fid);
            if(nc_ret == NC_EEXIST) {
                // File Exist. Open it
                nc_ret = nc_open_par(filename, NC_WRITE, DSGrid.gcomm, MPI_INFO_NULL, &nc_fid);
                if(nc_ret != NC_NOERR) {
                    fprintf(stderr, "ERROR: Rank %i: %s, line %i (%s): nc_open_par() failed ! Error code: %s\n",
                            DSGrid.grank, __FILE__, __LINE__, __func__, nc_strerror(nc_ret));
                }
            } else if(nc_ret != NC_NOERR) {
                fprintf(stderr, "ERROR: Rank %i: %s, line %i (%s): nc_create_par() failed ! Error code: %s\n",
                        DSGrid.grank, __FILE__, __LINE__, __func__, nc_strerror(nc_ret));
            }

            /* Define the dimensions. */
            for(int d=0; d<2; d++) {
                nc_ret = nc_def_dim(nc_fid, nc_dimname[d], DSGrid.yzgsize[d], &nc_dimid[d]);
                if(nc_ret != NC_NOERR) {
                    fprintf(stderr, "ERROR: Rank %i: %s, line %i (%s): nc_def_dim() failed ! Error code: %s\n",
                            DSGrid.grank, __FILE__, __LINE__, __func__, nc_strerror(nc_ret));
                }
            }

            for(int v=0; v<nslvar; v++) {
                var = var_index[v];
                sprintf(ds_var_name, "%s%s_%04d_%d", Run.path_2D, "yz_slice", ixpos[nsl], v);

                /* Define the netCDF variables. The dimids array is used to pass
                    the dimids of the dimensions of the variables.*/
                sprintf(nc_var_name, "%s.%06d", yz_slice_names[var], globiter);
                nc_ret = nc_def_var(nc_fid, nc_var_name, NC_FLOAT, 2, nc_dimid, &nc_varid[v]);
                if(nc_ret != NC_NOERR) {
                    fprintf(stderr, "ERROR: Rank %i: %s, line %i (%s): nc_def_var() failed ! Error code: %s\n",
                            DSGrid.grank, __FILE__, __LINE__, __func__, nc_strerror(nc_ret));
                }

                clk = MPI_Wtime();
                dspaces_get(client, ds_var_name, globiter, sizeof(float), 2, lb, ub, (void*) buffer, -1);
                time_get += MPI_Wtime() - clk;

                /* Write Data */
                clk = MPI_Wtime();
                nc_ret = nc_put_vara_float(nc_fid, nc_varid[v], nc_str, nc_lsize, buffer);   
                time_nc_file += MPI_Wtime() - clk;
                if(nc_ret != NC_NOERR) {
                    fprintf(stderr, "ERROR: Rank %i: %s, line %i (%s): nc_put_vara_float() failed ! Error code: %s\n",
                            DSGrid.grank, __FILE__, __LINE__, __func__, nc_strerror(nc_ret));
                }
            }
            
        }

        nc_ret = nc_close(nc_fid);
        if(nc_ret != NC_NOERR) {
            fprintf(stderr, "ERROR: Rank %i: %s, line %i (%s): nc_close() failed ! Error code: %s\n",
                    DSGrid.grank, __FILE__, __LINE__, __func__, nc_strerror(nc_ret));
        }
    }
    f_gdim = 0;
    free(buffer);
    free(nc_varid);
    free(ixpos);

    if(DSGrid.grank == 0) {
        fprintf(stdout, "Write YZ_Slice: GlobalIter = %d, Time of dspaces_get() = %lf, Time of NetCDF_File_write() = %lf.\n",
                globiter, time_get, time_nc_file);
    }
}

void Initialize(RunData& Run,GridData& Grid, PhysicsData& Physics, DSGridData& ds_Grid, MPI_Comm gcomm) {
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
            getvar(ds_Grid.procs, "dspaces_file_writer_3Dprocs", "int[3]", datafile);

            // 2D only read procs now, use the same decomposition for xy, xz, yz plane
            getvar(ds_Grid.procs2d, "dspaces_file_writer_2Dprocs", "int[2]", datafile);

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
    MPI_Bcast(ds_Grid.gsize,3,MPI_INT, 0,MPI_COMM_WORLD);
    MPI_Bcast(ds_Grid.procs2d,2,MPI_INT, 0,MPI_COMM_WORLD);

    MPI_Bcast(&Physics,sizeof(Physics),MPI_BYTE,0,MPI_COMM_WORLD);

    int nprocs, tot_procs = 1;
    MPI_Comm_size(gcomm, &nprocs);
    for(int d=0; d<ds_Grid.ndim; d++) {
        tot_procs *= ds_Grid.procs[d];
    }
    if (tot_procs != nprocs) {
        fprintf(stdout, "ERROR !!! : dspaces_file_writer_3Dprocs set in parameters.dat != number of processors!!\n");
        MPI_Abort(MPI_COMM_WORLD,1);
    }
    tot_procs = 1;
    for(int d=0; d<2; d++) {
        tot_procs *= ds_Grid.procs2d[d];
    }
    if (tot_procs != nprocs) {
        fprintf(stdout, "ERROR !!! : dspaces_file_writer_2Dprocs set in parameters.dat != number of processors!!\n");
        MPI_Abort(MPI_COMM_WORLD,1);
    }

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
            write_eos(client, Run, Physics, DSGrid, mdata->globiter);
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
            if(rank == 0) {
                fprintf(stdout, "Rank: %d: Write YZ_Slice: GlobalIter = %d ...\n", rank, mdata->globiter);
            }
            write_yz_slice(client, Run, Physics, DSGrid, mdata->globiter);
            if(rank == 0) {
                fprintf(stdout, "Rank: %d: Write YZ_Slice Done...\n", rank);
            }
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