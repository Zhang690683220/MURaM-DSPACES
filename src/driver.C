/* 
 *   Copyright (C) 1998-2000 University of Chicago. 
 *   See COPYRIGHT notice in top-level directory.
 */

#include <mpi.h>
#include "run.H"
#include "grid.H"
#include <stdlib.h>
#include <stdio.h>
#include "physics.H"
#include "solver.H"
#include "init.H"
#include "icbc.H"
#include "io.H"
#include "comm_split.H"
#include "rt/rt.h"

extern void ELTE_decon();

double start_time,rst_time;

int main(int argc, char** argv) {

  RunData Run;
  GridData Grid;
  PhysicsData Physics;
  
  RTS * rts = 0;
  
  double clock;

  MPI_Init(&argc,&argv);

  start_time=MPI_Wtime();
 
  if( Initialize(Run,Grid,Physics,rts)==restart ) {
    clock=MPI_Wtime();
    RestoreSolution(Run,Grid,Physics); 
    rst_time = MPI_Wtime()-clock;
    if (Run.rank == 0){
      cout << "Restored in " << rst_time << " seconds" << endl;
    }
  } else {
    if (Run.rank == 0) cout << Run.backfile << " not found. Aborting ... " << endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  ComputeSolution(Run,Grid,Physics,rts);

  comm_split_finalize();
  IO_Finalize();
  int ierr = MPI_Finalize();

  if(rts) delete rts;
  if (Physics.rt_ext[i_ext_cor]==2)
    ELTE_decon();
  if (Run.rank==0) fprintf(stderr,"deleted rts\n");
       
  return 0;
}
