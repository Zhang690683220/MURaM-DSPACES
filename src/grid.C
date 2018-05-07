#include <mpi.h>
#include <iostream>
#include "grid.H"
#include "run.H"
#include "comm_split.H"
using namespace std;

GridData::GridData() {
  NDIM = 3;

  v_nvar = 8;
  bufsize = 1;
  gnnodes = 1;
  cellvol = 1.;
  volume = 1.;
  vsize = 1;

  for(int i=0;i<3;i++) {
    gxmin[i] = 0.;
    gxmax[i] = 1.;
    gsize[i] = 1;
    ghosts[i] = 0;
    periods[i] = 0;
    pardim[i] = 0;
    procs[i] = 1;
  }
  
  U   = NULL;
  U0  = NULL;
  Res = NULL;
  
  //eos
  pres = NULL;
  temp = NULL;

  // RT
  Jtot = NULL;
  Stot = NULL;
  Qtot = NULL;
  Tau = NULL;

  // Thin Losses
  
  Qthin = NULL;
  QH = NULL;
  QMg = NULL;
  QCa = NULL;
  QCor = NULL;
  QChr = NULL;

  // Resistive and Viscous Heating
  Qres = NULL;
  Qvis = NULL;

  // Heat Conduction
  sflx0 = NULL;
  sflx = NULL;
  Rflx = NULL;
  BgradT = NULL;

  // DivB Cleaning
  divB = NULL;
  phi = NULL;

  // Current
  curlB = NULL;
  
  // ambipolar diffusion
  v_amb  = NULL;
  v0_amb = NULL;
  R_amb  = NULL;
  curlBxB = NULL;
  Qamb   = NULL;
  
  // Electron Number
  ne = NULL;

  // Collisional Frequencies
  amb = NULL;
  rhoi = NULL;

  // Temporary arrays for stuff
  tvar1 = NULL;
  tvar2 = NULL;
  tvar3 = NULL;
  tvar4 = NULL;
  tvar5 = NULL;
  tvar6 = NULL;
  tvar7 = NULL;
  tvar8 = NULL;
  
}

GridData::~GridData() {
    delete[] U;
    delete[] U0;
    delete[] Res;
    
    delete[] pres;
    delete[] temp;

    delete[] ne;
    delete[] amb;
    delete[] rhoi;

    delete[] Qtot;
    delete[] Stot;
    delete[] Jtot;
    delete[] Tau;

    delete[] Qthin;
    delete[] QH;
    delete[] QMg;
    delete[] QCa;
    delete[] QCor;
    delete[] QChr;

    delete[] Qres;
    delete[] Qvis;

    delete[] divB;
    delete[] phi;

    delete[] sflx0;
    delete[] sflx;
    delete[] Rflx;
    delete[] BgradT;

    delete[] curlB;

    delete[] v_amb;
    delete[] v0_amb;
    delete[] R_amb;
    delete[] curlBxB;
    delete[] Qamb;

    delete[] tvar1;
    delete[] tvar2;
    delete[] tvar3;
    delete[] tvar4;
    delete[] tvar5;
    delete[] tvar6;
    delete[] tvar7;
    delete[] tvar8;
}

void GridData::Init(const RunData &Run,const PhysicsData &Physics) {
  MPI_Datatype MPI_STATE;

  MPI_Type_contiguous(sizeof(cState),MPI_BYTE,&MPI_STATE);
  MPI_Type_commit(&MPI_STATE);
  
  int reorder = 0,ndim=3;
  MPI_Cart_create(MPI_COMM_WORLD,3,procs,periods,reorder,&cart_comm);
  MPI_Cart_coords(cart_comm,Run.rank,ndim,lrank);

  // who am I?
  rank = Run.rank;
 
  // find my neighbors
  int ierr[3];
  for (int dim=0 ; dim<3 ; dim++)
    ierr[dim]=MPI_Cart_shift(cart_comm,dim,1,&leftr[dim],&rightr[dim]);

  // Calculate size of local grid nx/ncores_x etc
  lsize[0] = (int) gsize[0]/procs[0];
  lsize[1] = (int) gsize[1]/procs[1];
  lsize[2] = (int) gsize[2]/procs[2];

  // My beginning point in the grid
  beg[0] = ghosts[0]+lrank[0]*lsize[0];
  beg[1] = ghosts[1]+lrank[1]*lsize[1];
  beg[2] = ghosts[2]+lrank[2]*lsize[2];

  // any remaining points
  int remx = gsize[0]-lsize[0]*procs[0];
  int remy = gsize[1]-lsize[1]*procs[1];
  int remz = gsize[2]-lsize[2]*procs[2];

  // distribute remaining points
  for (int ii=0;ii<min(remx,lrank[0]);ii++)
    beg[0]+=1;
  if (lrank[0]<remx)
    lsize[0]+=1;

  for (int ii=0;ii<min(remy,lrank[1]);ii++)
    beg[1]+=1;
  if (lrank[1]<remy)
    lsize[1]+=1;

  for (int ii=0;ii<min(remz,lrank[2]);ii++)
    beg[2]+=1;
  if (lrank[2]<remz)
    lsize[2]+=1;

  // Strides
  stride[0] = 1;
  stride[1] = (lsize[0]+2*ghosts[0]);
  stride[2] = stride[1]*(lsize[1]+2*ghosts[1]);
  
  // Beginning of local grid
  lbeg[0] = ghosts[0];
  lbeg[1] = ghosts[1];
  lbeg[2] = ghosts[2];

  // Beginning of full grid
  gbeg[0] = ghosts[0];
  gbeg[1] = ghosts[1];
  gbeg[2] = ghosts[2];

  // Am I the lower boundary
  if (lrank[0]==0)
    is_gbeg[0] = 1;
  else
    is_gbeg[0] = 0;

  if (lrank[1]==0)
    is_gbeg[1] = 1;
  else
    is_gbeg[1] = 0;

  if (lrank[2]==0)
    is_gbeg[2] = 1;
  else
    is_gbeg[2] = 0;
  
  // My end point in the grid
  end[0] = beg[0] + lsize[0] - 1;
  end[1] = beg[1] + lsize[1] - 1;
  end[2] = beg[2] + lsize[2] - 1;

  // My local end point
  lend[0] = lsize[0] + ghosts[0]-1;
  lend[1] = lsize[1] + ghosts[1]-1;
  lend[2] = lsize[2] + ghosts[2]-1;

  // Full grid end point
  gend[0] = gsize[0]+ghosts[0]-1;
  gend[1] = gsize[1]+ghosts[1]-1;
  gend[2] = gsize[2]+ghosts[2]-1;

  // Am I the upper boundary
  if (lrank[0]==procs[0]-1)
    is_gend[0] = 1;
  else
    is_gend[0] = 0;

  if (lrank[1]==procs[1]-1)
    is_gend[1] = 1;
  else
    is_gend[1] = 0;

  if (lrank[2]==procs[2]-1)
    is_gend[2] = 1;
  else
    is_gend[2] = 0;

  for(int d=0;d<3;d++) {
    int vs   = lsize[d]+2*ghosts[d];
    bufsize *= vs;
    gnnodes *= gsize[d];
    volume *= gxmax[d]-gxmin[d];
    vsize = vsize > vs ? vsize : vs;
  }

  U   = new cState[bufsize];
  U0  = new cState[bufsize];
  Res = new cState[bufsize];

  pres = new double[bufsize]();
  temp = new double[bufsize]();

  divB = new double[bufsize]();
  phi = new double[bufsize]();

  ne = new double[bufsize]();
  rhoi = new double[bufsize]();
  amb = new double[bufsize]();

  Qtot = new double[bufsize]();
  Jtot = new double[bufsize]();
  Stot = new double[bufsize]();
  Tau  = new double[bufsize]();

  if(Physics.rt_ext[i_ext_cor]>=1)
    Qthin = new double[bufsize]();
  if(Physics.rt_ext[i_ext_cor]>=2){
    QH   = new double[bufsize]();
    QMg  = new double[bufsize]();
    QCa  = new double[bufsize]();
    QCor = new double[bufsize]();
    QChr = new double[bufsize]();
  }
    
  if(Physics.params[i_param_spitzer] > 0.0){
    sflx0 = new double[bufsize]();
    sflx = new double[bufsize]();
    Rflx = new double[bufsize]();
    BgradT = new double[bufsize]();
  }

  if(Physics.params[i_param_eta] > 0.0){
    curlB = new Vector[bufsize]();
  }

  if(Physics.params[i_param_ambipolar] > 0.0){
    v_amb  = new Vector[bufsize]();
    v0_amb = new Vector[bufsize]();
    R_amb  = new Vector[bufsize]();
    curlBxB = new Vector[bufsize]();
  }

  if(Run.diagnostics){
    tvar1 = new double[bufsize]();
    tvar2 = new double[bufsize]();
    tvar3 = new double[bufsize]();
    tvar4 = new double[bufsize]();
    tvar5 = new double[bufsize]();
    tvar6 = new double[bufsize]();
    tvar7 = new double[bufsize]();
    tvar8 = new double[bufsize]();
    
    Qres = new double[bufsize]();
    Qvis = new double[bufsize]();
    
    if(Physics.params[i_param_ambipolar] > 0.0)
      Qamb    = new double[bufsize]();
  }
    
  for(int d=0;d<3;d++) {
    dx[d] = (gxmax[d]-gxmin[d])/gsize[d];
    lxmin[d] = gxmin[d]+(  beg[d]-gbeg[d])*dx[d];
    lxmax[d] = gxmin[d]+(1+end[d]-gbeg[d])*dx[d];
    cellvol *= dx[d];
  }
}

void GridData::Show() const {
  cout << " ------------ Grid Parameter Settings -------------" << endl;
  cout << "Decomposition = " << procs[0] << 'x' << procs[1] << 'x' << procs[2] << endl
       << "gsize         = " << gsize[0] << ' ' << gsize[1] << ' ' << gsize[2] << endl
       << "lsize         = " << lsize[0] << ' ' << lsize[1] << ' ' << lsize[2] << endl
       << "ghosts        = " << ghosts[0] << ' ' << ghosts[1] << ' ' << ghosts[2] << endl
       << "beg           = " << beg[0] << ' ' << beg[1] << ' ' << beg[2] << endl
       << "end           = " << end[0] << ' ' << end[1] << ' ' << end[2] << endl
       << "lbeg          = " << lbeg[0] << ' ' << lbeg[1] << ' ' << lbeg[2] << endl
       << "lend          = " << lend[0] << ' ' << lend[1] << ' ' << lend[2] << endl
       << "gbeg          = " << gbeg[0] << ' ' << gbeg[1] << ' ' << gbeg[2] << endl
       << "gend          = " << gend[0] << ' ' << gend[1] << ' ' << gend[2] << endl
       << "stride        = " << stride[0] << ' ' << stride[1] << ' ' << stride[2] << endl
       << "gxmin         = " << gxmin[0] << ' ' << gxmin[1] << ' ' << gxmin[2] << endl
       << "gxmax         = " << gxmax[0] << ' ' << gxmax[1] << ' ' << gxmax[2] << endl
       << "dx            = " << dx[0] << ' ' << dx[1] << ' ' << dx[2] << endl
       << "lxmin         = " << lxmin[0] << ' ' << lxmin[1] << ' ' << lxmin[2] << endl
       << "lxmax         = " << lxmax[0] << ' ' << lxmax[1] << ' ' << lxmax[2] << endl;
}
