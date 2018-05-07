#ifndef __RT_SCATTER_H__   // __RT_SCATTER_H__
#define __RT_SCATTER_H__

#include <mpi.h>
#include "grid.H"
#include "rt.h"
#include "run.H"
#include "physics.H"

class RTS_SCATTER:public RTS{
    protected:
        double ****S;
        double ** J_havg;
        double ***lambda_star,***lambda_old;
        double ** logtau_5000;
        double *** lgTau;
        int *** Tau_ind;
        int * ini_flag;
        virtual double wrapper(int,GridData&,RunData&,const PhysicsData&);
        void driver(double DZ, double DX, double DY, int band); 
        void get_lambdastar();
 
public:
  RTS_SCATTER(GridData & Grid,RunData & Run,PhysicsData & Physics);
  double Stot(int,int,int); 
  virtual ~RTS_SCATTER(void);
};

#endif                // __RT_SCATTER_H__
