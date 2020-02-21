#include <math.h>
#include <string.h>
#include <fstream>
#include <cmath>
#include <errno.h>
#include "mem.h"
#include "exchange.H"
#include "rt.h"
#include "rt_scatter.h"
#include "comm_split.H"

RTS *rt_new(GridData &Grid,RunData &Run,PhysicsData &Physics)
{
  int rttype;
  if(Run.rank==0){ // check solver type only
    rttype = Physics.rt[i_rt_type];
  }

  MPI_Bcast(&rttype,1,MPI_INT,0,MPI_COMM_WORLD);

  // Switch for chromospheric extension, three options here.
  // reads in parameters file
  // RT_DEFAULT - 0 - LTE rt, no scattering, optically thin and NLTE line losses can be incorporated.
  // SCATTER - 1 - rt with multigroup scattering as in Skartlein (2000) and Hayek (2010)
  // if nothing else return error

  switch(rttype){
      case(RT_DEFAULT): return new RTS(Grid,Run,Physics);
      case(RT_SCATTER): return new RTS_SCATTER(Grid,Run,Physics);
      default: return 0;
  }
}


