#ifndef __DS_FILE_WRITER_H__
#define __DS_FILE_WRITER_H__

#include "mpi.h"

class DSGridData{
    
public:
    MPI_Comm gcomm;
    MPI_Comm cart_comm;
    MPI_Comm comm2d;
    int grank;
    int ndim;
    int periods[3];

    /* 3D vars */
    int gsize[3];
    int lsize[3];
    int procs[3];
    int coord[3];
    int start[3];
    int end[3];

    /* 2D XY vars */
    int procs2d[2];
    int coord2d[2];

    int xygsize[2];
    int xylsize[2];
    int xyprocs[2];
    int xycoord[2];
    int xystart[2];
    int xyend[2];

    /* 2D XZ vars */
    int xzgsize[2];
    int xzlsize[2];
    int xzprocs[2];
    int xzcoord[2];
    int xzstart[2];
    int xzend[2];

    /* 2D YZ vars */
    int yzgsize[2];
    int yzlsize[2];
    int yzprocs[2];
    int yzcoord[2];
    int yzstart[2];
    int yzend[2];


    DSGridData();
    ~DSGridData() = default;
    void Init(MPI_Comm comm);
    void Show() const;

};

#endif
