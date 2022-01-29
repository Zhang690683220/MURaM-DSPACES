#ifndef __DS_FILE_WRITER_H__
#define __DS_FILE_WRITER_H__

#include "mpi.h"

class DSGridData{
    MPI_Comm gcomm;
    MPI_Comm cart_comm;
    int grank;
    int ndim;
    int gsize[3];
    int lsize[3];
    int periods[3];
    int procs[3];
    int coord[3];
    int start[3];
    int end[3];

    DSGridData();
    ~DSGridData() = default;
    void Init();
    void Show() const;

};

#endif