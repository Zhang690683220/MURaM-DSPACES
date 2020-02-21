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
#include "ACCH.h"

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


RTS::~RTS(void)
{
  for(int j = 0;j < cart_sizes[2]; j++) delete [] comm_col[j];
  delete [] comm_col;

///*
  ACCH::Free(I_o, ny*nx*sizeof(double)); // RT Grid
  ACCH::Free(tr_switch, nx*ny*nz*sizeof(int)); // switch for the lowest point at which we have reached the transition region
  ACCH::Free(lgTe, nx*ny*nz*sizeof(double)); // temperature
  ACCH::Free(lgPe, nx*ny*nz*sizeof(double)); // temperature
  ACCH::Free(T_ind, nx*ny*nz*sizeof(int)); // temperature
  ACCH::Free(P_ind, nx*ny*nz*sizeof(int)); // temperature
  ACCH::Free(rho, nx*ny*nz*sizeof(double));  // RT grid

  ACCH::Free(Tau, nx*ny*nz*sizeof(double));  // RT grid
  ACCH::Free3D<double>(Qt, ny-yo, nx-xo, nz-zo); // MHD grid

  ACCH::Free(Jt, nx*ny*nz*sizeof(double));
  ACCH::Free(St, nx*ny*nz*sizeof(double)); // RT grid

// frequency dependent quantities
  ACCH::Free(B, nx*ny*nz*sizeof(double));
  ACCH::Free(kap, nx*ny*nz*sizeof(double));
  
  if (rttype ==1){
    ACCH::Free(sig, nx*ny*nz*sizeof(double));
    ACCH::Free(abn, nx*ny*nz*sizeof(double));
  }
  
  ACCH::Free(J_band, nx*ny*nz*sizeof(double));
  ACCH::Free(I_n, nx*ny*nz*sizeof(double));

  ACCH::Free(Fx, nx*ny*nz*sizeof(double));
  ACCH::Free(Fy, nx*ny*nz*sizeof(double));
  ACCH::Free(Fz, nx*ny*nz*sizeof(double));

  ACCH::Free2D<double>(coeff, nx*ny*nz, 2);

  ACCH::Free3D<double>(Col_out, Nbands, col_nz, col_nvar);
   
  del_i5dim(numits,0,Nbands-1,FWD,BWD,RIGHT,LEFT,UP,DOWN,0,NMU-1);
  
  if (NDIM>1){
    ACCH::Free6D<real>(x_sbuf, Nbands, 2, 2, 2, NMU, nz*ny);
    ACCH::Free6D<real>(x_rbuf, Nbands, 2, 2, 2, NMU, nz*ny);
    ACCH::Free6D<real>(x_oldbuf, Nbands, 2, 2, 2, NMU, nz*ny);
  }

  if (NDIM==3){
    ACCH::Free6D<real>(y_sbuf, Nbands, 2, 2, 2, NMU, nz*nx);
    ACCH::Free6D<real>(y_rbuf, Nbands, 2, 2, 2, NMU, nz*nx);
    ACCH::Free6D<real>(y_oldbuf, Nbands, 2, 2, 2, NMU, nz*nx);
  }

  ACCH::Free6D<real>(z_sbuf, Nbands, 2, 2, 2, NMU, ny*nx);
  ACCH::Free6D<real>(z_rbuf, Nbands, 2, 2, 2, NMU, ny*nx);
  ACCH::Free6D<real>(z_oldbuf, Nbands, 2, 2, 2, NMU, ny*nx);

  ACCH::Free(tab_T, NT*sizeof(double));
  ACCH::Free(tab_p, Np*sizeof(double));
  ACCH::Free(invT_tab, NT*sizeof(double));
  ACCH::Free(invP_tab, Np*sizeof(double));
  
  if(fullodf) {
    ACCH::Free(nu_tab, Nlam*sizeof(float));
  }

  if (N5000){
    ACCH::Free(B_5000_tab, NT*sizeof(float));
    ACCH::Free2D<float>(kap_5000_tab, NT, Np);
  }

  ACCH::Free3D<float>(kap_tab, Nbands, NT, Np);
  ACCH::Free2D<float>(B_tab, Nbands, NT);

  ACCH::Free(I_band, nx*ny);

  ACCH::Delete(this, sizeof(RTS));

}
 
double RTS::tau(int z,int x,int y){
  double Tau_local =
    Tau[y*nx*nz+x*nz+z] +
    Tau[(y-yo)*nx*nz+x*nz+z] +
    Tau[y*nx*nz+(x-xo)*nz+z] +
    Tau[y*nx*nz+x*nz+(z-zo)] +
    Tau[(y-yo)*nx*nz+(x-xo)*nz+z] +
    Tau[(y-yo)*nx*nz+x*nz+(z-zo)] +
    Tau[y*nx*nz+(x-xo)*nz+(z-zo)] +
    Tau[(y-yo)*nx*nz+(x-xo)*nz+(z-zo)];

  Tau_local *= 0.125;

  return Tau_local;
}

double RTS::Qtot(int z,int x,int y)
{  
  return Qt[y-yl-yo][x-xl-xo][z-zl-zo];
}


double RTS::Jtot(int z,int x,int y){
  double J_local = Jt[y*nx*nz+x*nz+z] +
                   Jt[(y-yo)*nx*nz+x*nz+z] +
                   Jt[y*nx*nz+(x-xo)*nz+z] +
                   Jt[y*nx*nz+x*nz+(z-zo)] +
                   Jt[(y-yo)*nx*nz+(x-xo)*nz+z] +
                   Jt[(y-yo)*nx*nz+x*nz+(z-zo)] +
                   Jt[y*nx*nz+(x-xo)*nz+(z-zo)] +
                   Jt[(y-yo)*nx*nz+(x-xo)*nz+(z-zo)];
  J_local *= 0.125;

  return J_local;
}

double RTS::Stot(int z,int x,int y){
  double S_local = St[y*nx*nz+x*nz+z] +
                   St[(y-yo)*nx*nz+x*nz+z] +
                   St[y*nx*nz+(x-xo)*nz+z] +
                   St[y*nx*nz+x*nz+(z-zo)] +
                   St[(y-yo)*nx*nz+(x-xo)*nz+z] +
                   St[(y-yo)*nx*nz+x*nz+(z-zo)] +
                   St[y*nx*nz+(x-xo)*nz+(z-zo)] +
                   St[(y-yo)*nx*nz+(x-xo)*nz+(z-zo)];
  S_local *= 0.125;

  return S_local;
}


double RTS::Iout(int x,int y)
{
  double io = I_o[(y-yl)*nx+(x-xl)] + I_o[(y-yo-yl)*nx+(x-xl)] + I_o[(y-yl)*nx+(x-xo-xl)] + I_o[(y-yo-yl)*nx+(x-xo-xl)];
  io *= 0.25;
  return io;
}

RTS::RTS(GridData&Grid,RunData &Run,PhysicsData &Physics){
 
  std::cout << CACHE_SIZE << std::endl;

  NDIM = Grid.NDIM;
  rttype = Physics.rt[i_rt_type];
  eps_const = Physics.rt[i_rt_epsilon];

  verbose=Run.verbose;

  call_count=0;

  xl=Grid.lbeg[1];
  xh=Grid.lend[1];
  yl=Grid.lbeg[2];
  yh=Grid.lend[2];
  zl=Grid.lbeg[0];
  zh=Grid.lend[0];
//
  xl-=(xo=(Grid.lbeg[1]-Grid.lend[1])?1:0);
  yl-=(yo=(Grid.lbeg[2]-Grid.lend[2])?1:0);
  zl-=(zo=1);

//
  nx=xh-xl+1;
  ny=yh-yl+1;
  nz=zh-zl+1;
  myrank=Run.rank;

  for(int v=0;v<3;v++)
    next[v] = Grid.stride[v];

  int ndim=3,cart_periods[3];  

  MPI_Cart_get(cart_comm,ndim,cart_sizes,cart_periods,lrank);
  MPI_Cart_coords(cart_comm, Run.rank, ndim, lrank);
 
  for(int nleft,nright,nd=0;nd<ndim;nd++){
    MPI_Cart_shift(cart_comm,nd,1,&nleft,&nright);
    leftr[nd] = nleft;
    rightr[nd] = nright;
    
    isgbeg[nd]=Grid.is_gbeg[nd];
    isgend[nd]=Grid.is_gend[nd];
  }

  int ***colranks=i3dim(0,cart_sizes[2]-1,0,cart_sizes[1]-1,0,cart_sizes[0]-1);
  for(int j=0; j<cart_sizes[2];j++)
    for(int k=0; k<cart_sizes[1];k++)
      for(int i=0; i<cart_sizes[0];i++){
        int coords[3]={i,k,j};
        MPI_Cart_rank(cart_comm,coords,&(colranks[j][k][i]));
      }

  char carlson_fp[16];

  if (NMU == 3)
    strcpy(carlson_fp,"./carlson3.dat");
  else if (NMU == 6)
    strcpy(carlson_fp,"./carlson6.dat");
  else if (NMU == 10)
    strcpy(carlson_fp,"./carlson10.dat");
  else {
    fprintf(stdout," NMU %i not supported, Carlson 3,6,10 currently included.\n ",NMU);
    MPI_Abort(MPI_COMM_WORLD,1);
  }
  ifstream fptr(carlson_fp,ios::in);

  if(fptr){
    fptr.precision(16);
    for (int i=0; i < NMU; i++)
      fptr >> xmu[0][i] >> xmu[1][i] >> xmu[2][i] >> wmu[i];
    fptr.close();
  } else {
    fprintf(stdout,"carlson%i.dat not found. Aborting \n",NMU);
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  // load kappa bins

  ACCH::Copyin(this, sizeof(RTS));

  load_bins(Run.kap_name);

// output intensity
  I_o = (double*) ACCH::Malloc(nx*ny*sizeof(double));        // RT grid

  Fr_mean = (double*) ACCH::Malloc(Nbands*sizeof(double));
  gFr_mean = (double*) ACCH::Malloc(Nbands*sizeof(double));
  memset(Fr_mean,0,Nbands*sizeof(double));
  memset(gFr_mean,0,Nbands*sizeof(double));

  lgTe = (double*) ACCH::Malloc(nx*ny*nz*sizeof(double)); // temperature
  lgPe = (double*) ACCH::Malloc(nx*ny*nz*sizeof(double)); // pressure
  T_ind = (int*) ACCH::Malloc(nx*ny*nz*sizeof(int)); // pressure
  P_ind = (int*) ACCH::Malloc(nx*ny*nz*sizeof(int)); // pressure
  rho = (double*) ACCH::Malloc(nx*ny*nz*sizeof(double)); // RT grid
  tr_switch = (int*) ACCH::Malloc(nx*ny*nz*sizeof(int));

  Tau = (double*) ACCH::Malloc(nx*ny*nz*sizeof(double));  // RT grid

  Qt = ACCH::Malloc3D<double>(ny-yo, nx-xo, nz-zo); // MHD grid
  Jt = (double*) ACCH::Malloc(nx*ny*nz*sizeof(double)); // MHD grid
  St = (double*) ACCH::Malloc(nx*ny*nz*sizeof(double)); // MHD grid

// frequency dependent quantities
  B = (double*) ACCH::Malloc(nx*ny*nz*sizeof(double));
  kap = (double*) ACCH::Malloc(nx*ny*nz*sizeof(double));
  I_n = (double*) ACCH::Malloc(nx*ny*nz*sizeof(double));

  if (rttype==0){
    sig=kap;
    abn=kap;
  } else {
    sig = (double*) ACCH::Malloc(nx*ny*nz*sizeof(double));
    abn = (double*) ACCH::Malloc(nx*ny*nz*sizeof(double));
  }

  J_band = (double*) ACCH::Malloc(nx*ny*nz*sizeof(double));

  Fx = (double*) ACCH::Malloc(nx*ny*nz*sizeof(double));
  Fy = (double*) ACCH::Malloc(nx*ny*nz*sizeof(double));
  Fz = (double*) ACCH::Malloc(nx*ny*nz*sizeof(double));

  coeff = ACCH::Malloc2D<double>(nx*ny*nz, 2);
 
  /* Column outputs
   * 0: J_col - Angle averaged intensity
   * 1: S_col - Total Source function (=B in LTE)
   * 2: kap_col - total opacity
   * 3: abs_col - absorption opacity (= kap in LTE)
   * 4: sig_col - scattering opacity (= sig in LTE)
   * 5: B_col - Planck function (=S in LTE)
   * 6: tau_col - band dependant tau
   * 7: Qj_col - Q from radiative energy imbalance
   * 8: Qf_col - Q from flux divergence
   */

  col_offz=Grid.beg[0]-Grid.gbeg[0];
  col_nz = Grid.lsize[0];
  col_nzt = Grid.gsize[0];
  col_nvar = 9;

  Col_out = ACCH::Malloc3D<double>(Nbands, col_nz, col_nvar);

  memset(Col_out[0][0],0,col_nz*Nbands*col_nvar*sizeof(double));

  numits = i5dim(0,Nbands-1,FWD,BWD,RIGHT,LEFT,UP,DOWN,0,NMU-1);

  if (NDIM==3){
    y_sbuf = ACCH::Malloc6D<real>(Nbands, 2, 2, 2, NMU, nx*nz);
    y_rbuf = ACCH::Malloc6D<real>(Nbands, 2, 2, 2, NMU, nx*nz);
    y_oldbuf = ACCH::Malloc6D<real>(Nbands, 2, 2, 2, NMU, nx*nz);
    //y_sbuf=r6dim(0,Nbands-1,FWD,BWD,RIGHT,LEFT,UP,DOWN,0,NMU-1,0,nz*nx-1);
    //y_rbuf=r6dim(0,Nbands-1,FWD,BWD,RIGHT,LEFT,UP,DOWN,0,NMU-1,0,nz*nx-1);
    //y_oldbuf=r6dim(0,Nbands-1,FWD,BWD,RIGHT,LEFT,UP,DOWN,0,NMU-1,0,nz*nx-1);
  }

  if (NDIM>1){
    x_sbuf = ACCH::Malloc6D<real>(Nbands, 2, 2, 2, NMU, ny*nz);
    x_rbuf = ACCH::Malloc6D<real>(Nbands, 2, 2, 2, NMU, ny*nz);
    x_oldbuf = ACCH::Malloc6D<real>(Nbands, 2, 2, 2, NMU, ny*nz);
    //x_sbuf=r6dim(0,Nbands-1,FWD,BWD,RIGHT,LEFT,UP,DOWN,0,NMU-1,0,nz*ny-1);
    //x_rbuf=r6dim(0,Nbands-1,FWD,BWD,RIGHT,LEFT,UP,DOWN,0,NMU-1,0,nz*ny-1);
    //x_oldbuf=r6dim(0,Nbands-1,FWD,BWD,RIGHT,LEFT,UP,DOWN,0,NMU-1,0,nz*ny-1);
  }
  z_sbuf = ACCH::Malloc6D<real>(Nbands, 2, 2, 2, NMU, ny*nx);
  z_rbuf = ACCH::Malloc6D<real>(Nbands, 2, 2, 2, NMU, ny*nx);
  z_oldbuf = ACCH::Malloc6D<real>(Nbands, 2, 2, 2, NMU, ny*nx);
  //z_sbuf=r6dim(0,Nbands-1,FWD,BWD,RIGHT,LEFT,UP,DOWN,0,NMU-1,0,nx*ny-1);
  //z_rbuf=r6dim(0,Nbands-1,FWD,BWD,RIGHT,LEFT,UP,DOWN,0,NMU-1,0,nx*ny-1);
  //z_oldbuf=r6dim(0,Nbands-1,FWD,BWD,RIGHT,LEFT,UP,DOWN,0,NMU-1,0,nx*ny-1);
  //
 
  double DX=Grid.dx[1];
  double DY=Grid.dx[2];
  double DZ=Grid.dx[0];
  int perm[3][3]={{0,1,2},{1,2,0},{2,0,1}};
// Winkelabhaengige Koeff.
  for(int l=0;l<NMU;++l){ //  Koeff. fuer bilineare Grundflaecheninterp.
    double aM[]={DX/xmu[0][l],DY/xmu[1][l],DZ/xmu[2][l]};
    
    ibase[l]=0;
    ds_upw[l]=DX/xmu[0][l];

    if(((aM[2]<=aM[0])||(NDIM==1))&&((aM[2]<=aM[1])||(NDIM<3))){
      ibase[l]=2;
      ds_upw[l]=DZ/xmu[2][l];  
    }else{
      if((aM[1]<=aM[0])&&(NDIM==3)){
        ibase[l]=1;
        ds_upw[l]=DY/xmu[1][l];
      }
    }
    
    for(int i=0;i<=2;++i){
      a_00[i][l]=(1.0-(aM[perm[i][0]]/aM[perm[i][1]]))*(1.0-(aM[perm[i][0]]/aM[perm[i][2]]));
      a_01[i][l]=(aM[perm[i][0]]/aM[perm[i][2]])*(1.0-(aM[perm[i][0]]/aM[perm[i][1]]));
      a_10[i][l]=(aM[perm[i][0]]/aM[perm[i][1]])*(1.0-(aM[perm[i][0]]/aM[perm[i][2]]));
      a_11[i][l]=(aM[perm[i][0]]/aM[perm[i][1]])*(aM[perm[i][0]]/aM[perm[i][2]]);
    }
  } 

  if (NDIM>1){
    memset(x_sbuf[0][UP][RIGHT][FWD][0],0,Nbands*2*2*2*NMU*nz*ny*sizeof(real));
    memset(x_rbuf[0][UP][RIGHT][FWD][0],0,Nbands*2*2*2*NMU*nz*ny*sizeof(real));
    memset(x_oldbuf[0][UP][RIGHT][FWD][0],0,Nbands*2*2*2*NMU*nz*ny*sizeof(real));
  }

  if (NDIM==3){
    memset(y_sbuf[0][UP][RIGHT][FWD][0],0,Nbands*2*2*2*NMU*nz*nx*sizeof(real));
    memset(y_rbuf[0][UP][RIGHT][FWD][0],0,Nbands*2*2*2*NMU*nz*nx*sizeof(real));
    memset(y_oldbuf[0][UP][RIGHT][FWD][0],0,Nbands*2*2*2*NMU*nz*nx*sizeof(real));
  }

  memset(z_sbuf[0][UP][RIGHT][FWD][0],0,Nbands*2*2*2*NMU*nx*ny*sizeof(real));
  memset(z_rbuf[0][UP][RIGHT][FWD][0],0,Nbands*2*2*2*NMU*nx*ny*sizeof(real));
  memset(z_oldbuf[0][UP][RIGHT][FWD][0],0,Nbands*2*2*2*NMU*nx*ny*sizeof(real));

  for(int band = 0; band < Nbands; band++)
    for(int i = 0; i < 2; i++)
      for(int j = 0; j < 2; j++)
        for(int k = 0; k < 2; k++)
          for(int n = 0; n < NMU; n++) {
            if(!ACCH::Present(y_sbuf[band][i][j][k][n], nx*nz*sizeof(real))) std::cout << "Failed" << std::endl;
            if(!ACCH::Present(y_rbuf[band][i][j][k][n], nx*nz*sizeof(real))) std::cout << "Failed" << std::endl;
            if(!ACCH::Present(y_oldbuf[band][i][j][k][n], nx*nz*sizeof(real))) std::cout << "Failed" << std::endl;
            if(!ACCH::Present(x_sbuf[band][i][j][k][n], ny*nz*sizeof(real))) std::cout << "Failed" << std::endl;
            if(!ACCH::Present(x_rbuf[band][i][j][k][n], ny*nz*sizeof(real))) std::cout << "Failed" << std::endl;
            if(!ACCH::Present(x_oldbuf[band][i][j][k][n], ny*nz*sizeof(real))) std::cout << "Failed" << std::endl;
            if(!ACCH::Present(z_sbuf[band][i][j][k][n], nx*ny*sizeof(real))) std::cout << "Failed" << std::endl;
            if(!ACCH::Present(z_rbuf[band][i][j][k][n], nx*ny*sizeof(real))) std::cout << "Failed" << std::endl;
            if(!ACCH::Present(z_oldbuf[band][i][j][k][n], nx*ny*sizeof(real))) std::cout << "Failed" << std::endl;
          }

  memset(numits[0][UP][RIGHT][FWD],0,2*2*2*Nbands*NMU*sizeof(int));

  // init for qrad_tauscale
  int* ranks=new int [cart_sizes[0]];
  MPI_Group MPI_GROUP_WORLD;

  MPI_Comm_group(MPI_COMM_WORLD,&MPI_GROUP_WORLD);
  MPI_Group ** grp_col = new MPI_Group * [cart_sizes[2]];
  comm_col = new MPI_Comm * [cart_sizes[2]];
  
  for(int j=0; j<cart_sizes[2];j++){
    grp_col[j]=new MPI_Group [cart_sizes[1]];
    comm_col[j] = new MPI_Comm [cart_sizes[1]];
  }
  
  for(int j=0; j< cart_sizes[2]; j++){
    for(int k=0; k< cart_sizes[1]; k++){
      for(int i=0; i< cart_sizes[0]; i++)
        ranks[i]=colranks[j][k][cart_sizes[0]-1-i];
      MPI_Group_incl(MPI_GROUP_WORLD,cart_sizes[0],ranks,&(grp_col[j][k]));
      MPI_Comm_create(MPI_COMM_WORLD,grp_col[j][k],&(comm_col[j][k]));
    }
  }
  delete[] ranks;
// all done...?
  del_i3dim(colranks,0,cart_sizes[2]-1,0,cart_sizes[1]-1,0,cart_sizes[0]-1);

  for(int j = 0;j < cart_sizes[2]; j++) delete [] grp_col[j];
  delete [] grp_col;

  /* if (xl or xh) lies in my x-domain and yl or yh lies in my y-domain
   * then save_col, set limits */

  int gxl = Grid.beg[1];
  int gxh = Grid.end[1];
  int gyl = Grid.beg[2];
  int gyh = Grid.end[2];

  int sl_r[4] = {Grid.gbeg[1],Grid.gend[1],Grid.gbeg[2],Grid.gend[2]};

  int x_range = (((gxl >= sl_r[0])&&(gxl <= sl_r[1]))
      ||((gxh >=sl_r[0])&&(gxh <= sl_r[1])));
  int y_range = (((gyl >= sl_r[2])&&(gyl <= sl_r[3]))
      ||((gyh >=sl_r[2])&&(gyh <=sl_r[2])));

  save_col = 0;
  col_bnd[0] =1;
  col_bnd[1] =0;
  col_bnd[2] =1;
  col_bnd[3] =0;

  num_col = (sl_r[1]-sl_r[0]+1)*(sl_r[3]-sl_r[2]+1);
  avg_col = 1.0/((double) num_col);

  if ((x_range)&&(y_range)){
    save_col=1;
    col_bnd[0] =max(gxh-sl_r[0],xl+xo);
    col_bnd[1] =min(gxh-sl_r[1],xh);
    col_bnd[2] =max(gyh-sl_r[2],xl+xo);
    col_bnd[3] =min(gyh-sl_r[3],xh);
  }

  I_band = (double*) ACCH::Malloc(nx*ny*sizeof(double));

  ACCH::UpdateGPU(this, sizeof(RTS));

}
void RTS::load_bins(char* kap_name){
  
  ifstream fp_rt(kap_name, ios::in | ios::binary);

  int pt_rhot;
  int junk4;

  if (fp_rt.is_open()) {

    fp_rt.read((char*)&N5000, sizeof(int));
    fp_rt.read((char*)&NT, sizeof(int));
    fp_rt.read((char*)&Np, sizeof(int));
    fp_rt.read((char*)&Nbands, sizeof(int));
    fp_rt.read((char*)&pt_rhot, sizeof(int));
    fp_rt.read((char*)&fullodf, sizeof(int));
    fp_rt.read((char*)&scatter, sizeof(int));
    fp_rt.read((char*)&junk4, sizeof(int)); 
     
    if (pt_rhot!=0){
      fprintf(stdout,"pt_rhot = %i, but rho-T bins are not currently implemented, aborting",pt_rhot);
      MPI_Abort(MPI_COMM_WORLD,1);
    }

    if (myrank == 0) {
      fprintf(stdout,"Reading in RTS Header, current settings:\n");
      fprintf(stdout,"Full ODF is %d and scatter is %d \n",fullodf,scatter);
      fprintf(stdout,"RT bins are NT %d Np %d Nbands %d \n",NT,Np,Nbands);
      fprintf(stdout,"Reference bin is %d and coronal back heating bin %d \n",N5000,junk4);
    }
    
    if ((scatter==0)&&(eps_const==0)&&(rttype==1)){
      fprintf(stdout, "rt_type =%i, but scattering bins = %i and photon destruction probability is %e. \n",scatter, rttype,eps_const);
      fprintf(stdout, "For scattering either scattering bins, or constant photon destruction probability are required, aborting. \n");
      MPI_Abort(MPI_COMM_WORLD,1);
    }

    tab_T = (double*) ACCH::Malloc(NT*sizeof(double));
    tab_p = (double*) ACCH::Malloc(Np*sizeof(double));

    invT_tab = (double*) ACCH::Malloc(NT*sizeof(double));
    invP_tab = (double*) ACCH::Malloc(Np*sizeof(double));
    
    fp_rt.read((char*)&tab_T[0],NT*sizeof(double));
    fp_rt.read((char*)&tab_p[0],Np*sizeof(double));
    
    for (int i=0;i<NT; i++)
      tab_T[i] = tab_T[i]*TENLOG;
     
    for (int j=0;j<Np; j++)
      tab_p[j] = tab_p[j]*TENLOG;

    for (int l=0; l<=NT-2; l++)
      invT_tab[l]= 1./(tab_T[l+1] - tab_T[l]);
    for (int m=0; m<=Np-2; m++)
      invP_tab[m] = 1./ (tab_p[m+1]- tab_p[m]);

    if (N5000){
      kap_5000_tab = ACCH::Malloc2D<float>(NT, Np);
      B_5000_tab = (float*) ACCH::Malloc(NT*sizeof(float));
       
      fp_rt.read((char*)&kap_5000_tab[0][0],NT*Np*sizeof(float));
      fp_rt.read((char*)&B_5000_tab[0],NT*sizeof(float));
    }
      
    kap_tab = ACCH::Malloc3D<float>(Nbands, NT, Np);
    B_tab = ACCH::Malloc2D<float>(Nbands, NT);
    fp_rt.read((char*)&kap_tab[0][0][0],Nbands*NT*Np*sizeof(float));
    fp_rt.read((char*)&B_tab[0][0],Nbands*NT*sizeof(float));

    if(fullodf){
      nu_tab = f1dim(0,Nlam-1);
      fp_rt.read((char*)&nu_tab[0],Nlam*sizeof(float));

      if (scatter>0){
        acont_pT = f3dim(0,Nlam-1,0,NT-1,0,Np-1);
        kcont_pT = f3dim(0,Nlam-1,0,NT-1,0,Np-1);
        fp_rt.read((char*)&acont_pT[0][0][0],Nlam*NT*Np*sizeof(float));
        fp_rt.read((char*)&kcont_pT[0][0][0],Nlam*NT*Np*sizeof(float));
        
        Npp=1;
        tau_pp_tab = d1dim(0,Npp-1);
        invtau_pp_tab = d1dim(0,Npp-1);

        tau_pp_tab[0] = -99.0;
        invtau_pp_tab[0] = 1.0;

        // If scattering bins and full ODF, but we want to run non-scattering RT,
        // we need to add acont to kappa
        if (rttype==0)
          for (int lam = 0 ;lam < Nlam;lam++)
            for (int bin = 0;bin<Nbin;bin++)
              for (int tt = 0;tt<NT;tt++)
                for (int pp = 0;pp<Np;pp++)
                  kap_tab[Nbin*lam+bin][tt][pp] = log(exp(kap_tab[Nbin*lam+bin][tt][pp])
                      +acont_pT[lam][tt][pp]);
      }
    } else if (scatter>0){
        // output
        Npp = scatter;
        scatter=1;

        sig_tab = f3dim(0,Nbands-1,0,NT-1,0,Np-1);
        abn_tab = f3dim(0,Nbands-1,0,NT-1,0,Np-1);
        fp_rt.read((char*)&abn_tab[0][0][0],Nbands*NT*Np*sizeof(float));
        fp_rt.read((char*)&sig_tab[0][0][0],Nbands*NT*Np*sizeof(float));

        tau_pp_tab = d1dim(0,Npp-1);
        invtau_pp_tab = d1dim(0,Npp-1);

        fp_rt.read((char*)&tau_pp_tab[0],Npp*sizeof(double));

        for (int z=0;z<=Npp-2; z++)
          invtau_pp_tab[z] = 1./(tau_pp_tab[z+1]- tau_pp_tab[z]);
        
        kap_pp_tab = f2dim(0,Nbands-1,0,Npp-1);
        abn_pp_tab = f2dim(0,Nbands-1,0,Npp-1);
        sig_pp_tab = f2dim(0,Nbands-1,0,Npp-1);

        fp_rt.read((char*)&abn_pp_tab[0][0],Nbands*Npp*sizeof(float));
        fp_rt.read((char*)&kap_pp_tab[0][0],Nbands*Npp*sizeof(float));
        fp_rt.read((char*)&sig_pp_tab[0][0],Nbands*Npp*sizeof(float));
    }
        
    if ((scatter==0)&&(rttype==1)){
      Npp=1;
      tau_pp_tab = d1dim(0,Npp-1);
      invtau_pp_tab = d1dim(0,Npp-1);

      tau_pp_tab[0] = -99;
      invtau_pp_tab[0] = 1;
    }

    fp_rt.close();
     
  } else {
    cout << "rt_init: kappa file not found: " << kap_name << endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }
 
}

double RTS::wrapper(int rt_upd,GridData &Grid,RunData &Run,const PhysicsData &Physics){
  //  return 1;
  
  //double * I_band;
  //I_band = (double*) ACCH::Malloc(nx*ny*sizeof(double));        // RT grid

  memset(I_band,0,nx*ny*sizeof(double));
  memset(I_o,0,nx*ny*sizeof(double));

  double DX=Grid.dx[1],DZ=Grid.dx[0],DY=Grid.dx[2];
  int cont_bin = Physics.rt[i_rt_iout];

  const double Temp_TR = Physics.rt[i_rt_tr_tem];
  const double Pres_TR = Physics.rt[i_rt_tr_pre];

  need_I = 0;
  // Do i need to plot? Should I calculate I_o?
  if (Run.NeedsSlice())
    need_I = 1;

  if((Run.iteration%rt_upd)&&(dt_rad>0)){
    // if cont_bin = 2 and need_I = 1 we want bolometric intensity
    if (cont_bin == 2){
    for (int band=Nbands-1;band>=0;--band){
      get_Tau_and_Iout(Grid, Run, Physics,DZ,B_tab[band],kap_tab[band],I_band,need_I);
      ACCH::UpdateGPU(I_o, nx*ny*sizeof(double));
      ACCH::UpdateGPU(I_band, nx*ny*sizeof(double));
#pragma acc parallel loop collapse(2) \
 present(this[:1], I_o[:nx*ny], I_band[:nx*ny])
      for (int y=0;y<ny;y++)
        for (int x=0;x<nx;x++)
          I_o[y*nx+x] +=I_band[y*nx+x];
      }
      ACCH::UpdateCPU(I_o, nx*ny*sizeof(double));
    }
    // If cont_bin = 0 and need_I = 1 we want the output intensity to be the 0th continuum bin
    if (cont_bin==0){
      get_Tau_and_Iout(Grid, Run, Physics,DZ,B_tab[0],kap_tab[0],I_band,need_I);
      ACCH::UpdateGPU(I_o, nx*ny*sizeof(double));
      ACCH::UpdateGPU(I_band, nx*ny*sizeof(double));
#pragma acc parallel loop collapse(2) \
 present(this[:1], I_o[:nx*ny], I_band[:nx*ny])
      for (int y=0;y<ny;y++)
        for (int x=0;x<nx;x++)
          I_o[y*nx+x] +=I_band[y*nx+x];
      ACCH::UpdateCPU(I_o, nx*ny*sizeof(double));
    }
    // if we have the band, we want a tau5000 reference wavelength. If cont_bin=1 and need_I we want
    // output intensity to be 5000A.
    
    if (N5000){
      int I5000_out = 0;
      if ((cont_bin==1)&&(need_I==1))
        I5000_out = 1;
      
      get_Tau_and_Iout(Grid, Run, Physics,DZ,B_5000_tab,kap_5000_tab,I_band,I5000_out);
     
      if ((cont_bin==1)&&(need_I==1)) {
        ACCH::UpdateGPU(I_o, nx*ny*sizeof(double));
        ACCH::UpdateGPU(I_band, nx*ny*sizeof(double));
#pragma acc parallel loop collapse(2) \
 present(this[:1], I_o[:nx*ny], I_band[:nx*ny])
        for (int y=0;y<ny;y++)
          for (int x=0;x<nx;x++)
            I_o[y*nx+x] += I_band[y*nx+x];
        ACCH::UpdateCPU(I_o, nx*ny*sizeof(double));
      }
    }

    calc_Qtot_and_Tau(Grid, Run, Physics);
    //ACCH::Free(I_band, nx*ny*sizeof(double));
    return dt_rad;
  }

// *****************************************************************
// *        interpolate opacity and Planck function (B)            *
// *****************************************************************

  cState *U=Grid.U;

  double N = pow(2,NDIM);

  for(int y=0;y<ny;++y){
    for(int x=0;x<nx;++x){
      int y0 = y+yl;
      int x0 = x+xl;
      int off0 = x0*next[1]+y0*next[2];
      int off1 = x0*next[1]+(y0+yo)*next[2];
      int off2 = (x0+xo)*next[1]+y0*next[2];
      int off3 = (x0+xo)*next[1]+(y0+yo)*next[2];
      int xyoff = y*nx*nz+x*nz;
      for(int z=0;z<nz;++z){
        int z0 = z+zl;
        int inode[]={off0+z0,off0+z0+zo,off2+z0,off2+z0+zo,off1+z0,off1+z0+zo,off3+z0,off3+z0+zo};
        double Tm=0.0,pm=0.0,rm=0.0;
        for(int l=0;l<N;++l){
          Tm+=Grid.temp[inode[l]];
          pm+=Grid.pres[inode[l]];
          rm+=U[inode[l]].d;
        }

        Tm /= N;
        pm /= N;
        rm /= N;

        int ind = xyoff + z;

        //disbale RT if Temp > Temp_TR above the photosphere
        double pswitch = min(max(pm-Pres_TR,0.0),1.0);
        double tswitch = min(max(Temp_TR-Tm,0.0),1.0);
        tr_switch[ind] = (int) max(pswitch,tswitch);
        
        lgTe[ind]=log(Tm);
        lgPe[ind]=log(pm);
        rho[ind] =rm;
    
        lgTe[ind] = min(max(lgTe[ind],(double) tab_T[0]),(double) tab_T[NT-1]);
        lgPe[ind] = min(max(lgPe[ind],(double) tab_p[0]),(double) tab_p[Np-1]);
      }
    
      for(int z=0;z<nz;++z){

        int ind = xyoff + z;

        // Search table once for all RT bands
        int l=0;
        int m=0;
        if(lgTe[ind]<tab_T[0])
          l=0;
        else if(lgTe[ind]>tab_T[NT-1])
          l=NT-2;
        else
          for (l=0; l<=NT-2; l++)
            if ((lgTe[ind] >= tab_T[l]) && (lgTe[ind] <= tab_T[l+1]))
              break;

        if(lgPe[ind]<tab_p[0])
          m=0;
        else if(lgPe[ind]>tab_p[Np-1])
          m=Np-2;
        else
          for (m=0; m<=Np-2; m++)
             if ((lgPe[ind] >= tab_p[m]) && (lgPe[ind] <= tab_p[m+1]))
               break;

        T_ind[ind] = l;
        P_ind[ind] = m;

      }
    }
  }
  
  // Begin RT loop, work band by band,
  memset(St,0,nx*ny*nz*sizeof(double));
  memset(Jt,0,nx*ny*nz*sizeof(double));
  memset(Qt[0][0],0,(nx-xo)*(ny-yo)*(nz-zo)*sizeof(double)); // MHD grid
 
  for(int band=Nbands-1;band>=0;--band){

    if(fullodf){
      nu_ind = band/Nbin;
      bin_ind = band%Nbin;
    }

    if((myrank==0) && (verbose>1)){
      fprintf(stdout,"rt running for band %i of %i \n",band+1,Nbands);
      if (fullodf)
        fprintf(stdout,"nu band = %i, %e and bin %i \n", nu_ind, nu_tab[nu_ind],bin_ind);
    }

    for(int y=0;y<ny;++y){
      for(int x=0;x<nx;++x){
        int xyoff = y*nx*nz + x*nz;
        for(int z=0;z<nz;++z){
          int ind = xyoff + z;
      
          int l = T_ind[ind];
          int m = P_ind[ind];
      
          double xt = (lgTe[ind]-tab_T[l])*invT_tab[l];
          double xp = (lgPe[ind]-tab_p[m])*invP_tab[m];
      
          // Interpolate for kappa and B
          B[ind]=exp(xt*B_tab[band][l+1]+(1.-xt)*B_tab[band][l]);
      
          kap[ind] = 
            exp(xt*(xp*kap_tab[band][l+1][m+1]+(1.-xp)*kap_tab[band][l+1][m])+
            (1.-xt)*(xp*kap_tab[band][l][m+1]+(1.-xp)*kap_tab[band][l][m]));
    
          // TR switch turn off kappa and B
          kap[ind] *= tr_switch[ind];
          B[ind]   *= tr_switch[ind];
        }
      }
    }

// *****************************************************************
// *    update diffusion-approx. boundary condition at bottom      *
// *****************************************************************
  for(int YDIR=FWD;YDIR<=BWD;++YDIR)
    for(int XDIR=RIGHT;XDIR<=LEFT;++XDIR) {
      if(isgbeg[0]==1)
        for(int l=0;l<NMU;++l) {
          for(int y=0;y<ny;++y)
            for(int x=0;x<nx;++x)
              z_rbuf[band][YDIR][XDIR][UP][l][x*ny+y]=B[y*nx*nz+x*nz];
        } // end l
// *****************************************************************
// *    no incoming radiation on the top                           *
// *****************************************************************
      if(isgend[0]==1) {
        memset(z_rbuf[band][YDIR][XDIR][DOWN][0],0,NMU*nx*ny*sizeof(real));
      } // end isgend[0]
    } // end X
// *****************************************************************
// *  solve the transfer                                           *
// *****************************************************************
  driver(DZ,DX,DY,band); 

// *****************************************************************
// *  outgoing radiative flux                                      *
// *****************************************************************

  gFr_mean[band] = 0.0;

  if(isgend[0]==1){
    for(int y=0; y<ny-yo; y++)
      for(int x=0; x<nx-xo; x++)
        gFr_mean[band]+=Fz[y*nx*nz + x*nz + nz-1] +
                        Fz[(y+yo)*nx*nz + x*nz + nz-1] +
                        Fz[y*nx*nz + (x+xo)*nz + nz-1] +
                        Fz[(y+yo)*nx*nz + (x+xo)*nz + nz-1];
    gFr_mean[band]*=0.25;
  }
 
  // If I am saving 5000A intensity then wipe after last bin. if I am saving the continuum bin then wipe the second last bin.

  if (((cont_bin==1)&&(band==0))||((cont_bin==0)&&(band==1))) {
    ACCH::UpdateGPU(I_o, nx*ny*sizeof(double));
    ACCH::UpdateGPU(I_band, nx*ny*sizeof(double));
#pragma acc parallel loop collapse(2) \
 present(this[:1], I_o[:nx*ny], I_band[:nx*ny])
    for(int y=0;y<ny;++y)
      for(int x=0;x<nx;++x)
        I_o[y*nx+x] = 0.0;
    ACCH::UpdateCPU(I_o, nx*ny*sizeof(double));
  }

  }// end loop over bands

  MPI_Allreduce(&gFr_mean[0],&Fr_mean[0],Nbands,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  F_o = 0.0;
 

  for (int band=0;band<Nbands;++band){
    Fr_mean[band]/=(Grid.gsize[1]*Grid.gsize[2]);
    F_o+=Fr_mean[band];
  }
  
  // If I have the band calculate tau5000 reference grid. If I want 5000 A continuum
  // output then switch that on too.
  
  if (N5000){
    int I5000_out = 0;
    if ((cont_bin==1)&&(need_I==1))
      I5000_out = 1;
    
    get_Tau_and_Iout(Grid, Run, Physics,DZ,B_5000_tab,kap_5000_tab,I_band,I5000_out);
   
    if ((cont_bin==1)&&(need_I==1)) {
      for (int y=0;y<ny;y++)
        for (int x=0;x<nx;x++)
          I_o[y*nx+x] +=I_band[y*nx+x];
      PGI_COMPARE(I_o, double, ny*nx, "I_o", "rt.cc", "RTS::wrapper", 14)
    }
  }

  calc_Qtot_and_Tau(Grid, Run, Physics);

  //ACCH::Free(I_band, nx*ny*sizeof(double));
  
  if (Run.NeedsSlice() && Run.RT_HAVG)
    save_1D_avg(Run.path_2D,Run.globiter,Run.time); 

  PGI_COMPARE(&dt_rad, double, 1, "dt_rad", "rt.cc", "RTS::wrapper", 15) 

  return dt_rad;
}

void RTS::calc_Qtot_and_Tau(GridData &Grid, const RunData &Run, const PhysicsData &Physics){

 double tau_min = pow(Physics.rt[i_rt_tau_min],2);
//DP - Tabulated NLTE extension

 cState *U = Grid.U;

// *****************************************************************
// * Limit Qtot if dr_rad<rt_tstep, consider only region with      *
// * rho > rt_rho_min for time step computation. Compute Fout      *
// * from Qrad                                                     *
// *****************************************************************
  dt_rad=0.0;
  double _dt_rad = 0.0;
  double qsum=0.0;

  ACCH::UpdateGPU(Grid.Tau, Grid.bufsize*sizeof(double));
  ACCH::UpdateGPU(Tau, nx*ny*nz*sizeof(double));
  ACCH::UpdateGPU(Grid.Jtot, Grid.bufsize*sizeof(double));
  ACCH::UpdateGPU(Jt, nx*ny*nz*sizeof(double));
  ACCH::UpdateGPU(Grid.Stot, Grid.bufsize*sizeof(double));
  ACCH::UpdateGPU(St, nx*ny*nz*sizeof(double));
  ACCH::UpdateGPU(tr_switch, nx*ny*nz*sizeof(int));
  ACCH::UpdateGPU(U, Grid.bufsize*sizeof(cState));
  ACCH::UpdateGPU3D<double>(Qt, ny-yo, nx-xo, nz-zo);
  ACCH::UpdateGPU(Grid.Qtot, Grid.bufsize*sizeof(double));
#pragma acc parallel loop gang collapse(2) \
  present(this[:1], Grid[:1], Grid.Tau[:Grid.bufsize], Tau[:nx*ny*nz], \
          Grid.Jtot[:Grid.bufsize], Jt[:nx*ny*nz], Grid.Stot[:Grid.bufsize], \
	  St[:nx*ny*nz], tr_switch[:nx*ny*nz], U[:Grid.bufsize], Qt[:ny-yo][:nx-xo][:nz-zo], Grid.Qtot[:Grid.bufsize]) \
  copyin(next[1:2]) reduction(+:qsum) reduction(max:_dt_rad)
  for(int y = yo; y < ny; y++)
    for(int x = xo; x < nx; x++) {
      int off0 = (x+xl)*next[1]+(y+yl)*next[2];
      int xyoff = y*nx*nz + x*nz;
#pragma acc loop vector reduction(+:qsum) reduction(max:dt_rad)
      for(int z = zo; z < nz; z++) {
        int node = off0+z+zl;
	Grid.Tau[node] = tau(z, x, y);
	Grid.Jtot[node] = Jtot(z, x, y);
	Grid.Stot[node] = Stot(z, x, y);
	int ind = xyoff + z;
	double scale = pow(Grid.Tau[node], 2);
	scale = scale/(scale + tau_min)*tr_switch[ind];
	double Qt_step = Qt[y-yo][x-xo][z-zo]*scale;
	double inv_dt = fabs(Qt_step)/U[node].e;
	_dt_rad = max(_dt_rad, inv_dt);
	Grid.Qtot[node] = Qt_step;
	qsum += Qt_step;
      }
    }
  dt_rad = _dt_rad;
  ACCH::UpdateCPU(Grid.Tau, Grid.bufsize*sizeof(double));
  ACCH::UpdateCPU(Grid.Jtot, Grid.bufsize*sizeof(double));
  ACCH::UpdateCPU(Grid.Stot, Grid.bufsize*sizeof(double));
  ACCH::UpdateCPU(Grid.Qtot, Grid.bufsize*sizeof(double));

  exchange_single(Grid,Grid.Tau);

  PGI_COMPARE(Grid.Tau, double, Grid.bufsize, "Tau", "rt.cc", "RTS::calc_Qtot_and_Tau", 21)

  double Fqrad;

  MPI_Allreduce(&qsum,&Fqrad,1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  Fqrad=-Fqrad*Grid.dx[0]/(Grid.gsize[1]*Grid.gsize[2]);

  if(myrank==0) fprintf(stdout,"RT energy flux: %21.15E %21.15E\n",F_o,Fqrad);
  
  F_o=Fqrad;
  dt_rad=Physics.rt[i_rt_cfl]/dt_rad;
  
}

void RTS::integrate(
  double ** coeff, const double c[4], const int stride[2], int ystep, int xstep, int zstep, int yi_i, int yi_f, int xi_i, int xi_f, int zi_i, int zi_f
)
{
  int off[4];
  int i_nu=0;
  double *ii=I_n;
  if(NDIM==3){

    for(int i=0;i<4;i++) off[i]=iystep[i]*stride[0]+ixstep[i]*stride[1]+izstep[i];

    bool ix = (ixstep[0] == 1 && ixstep[1] == 1    && 
               ixstep[2] == 1 && ixstep[3] == 1)   ||
              (ixstep[0] == -1 && ixstep[1] == -1  && 
               ixstep[2] == -1 && ixstep[3] == -1);
    bool iy = (iystep[0] == 1 && iystep[1] == 1    && 
               iystep[2] == 1 && iystep[3] == 1)   ||
              (iystep[0] == -1 && iystep[1] == -1  && 
               iystep[2] == -1 && iystep[3] == -1);
    bool iz = (izstep[0] == 1 && izstep[1] == 1    && 
               izstep[2] == 1 && izstep[3] == 1)   ||
              (izstep[0] == -1 && izstep[1] == -1  && 
               izstep[2] == -1 && izstep[3] == -1);
//    ACCH::UpdateGPU(ii, nx*ny*nz*sizeof(double));

    const int str_inu[] = {(nx-1)*(nz-1), nz-1, 1};
//#pragma acc data copyin(stride[:2], str_inu[:2], off[:4], c[:4])
{
    if(iy) {
      for(int y = 0; y < ny-1; y++) {
        int yi = (yi_i-yl) + y*ystep;
//#pragma acc parallel loop gang async \
// present(this[:1], ii[:nx*ny*nz], coeff[:nx*ny*nz][:2], stride[:2], str_inu[:2], off[:4], c[:4])
        for(int x = 0; x < nx-1; x++) {
          int xi = (xi_i-xl) + x*xstep;
//#pragma acc loop vector
          for(int z = 0; z < nz-1; z++) {
            int zi = (zi_i-zl) + z*zstep;
            int ind = yi*stride[0] + xi*stride[1] + zi;
            int _i_nu = y*str_inu[0] + x*str_inu[1] + z;
            double I_upw=c[0]*ii[ind-off[0]]+c[1]*ii[ind-off[1]]+c[2]*ii[ind-off[2]]+c[3]*ii[ind-off[3]];
            ii[ind]=I_upw*coeff[_i_nu][0]+coeff[_i_nu][1];
          }
        }
      }
    } else if(ix) {
      for(int x = 0; x < nx-1; x++) {
        int xi = (xi_i-xl) + x*xstep;
//#pragma acc parallel loop gang async \
// present(this[:1], ii[:nx*ny*nz], coeff[:nx*ny*nz][:2], stride[:2], str_inu[:2], off[:4], c[:4])
        for(int y = 0; y < ny-1; y++) {
          int yi = (yi_i-yl) + y*ystep;
//#pragma acc loop vector
	  for(int z = 0; z < nz-1; z++) {
            int zi = (zi_i-zl) + z*zstep;
            int ind = yi*stride[0] + xi*stride[1] + zi;
            int _i_nu = y*str_inu[0] + x*str_inu[1] + z;
            double I_upw=c[0]*ii[ind-off[0]]+c[1]*ii[ind-off[1]]+c[2]*ii[ind-off[2]]+c[3]*ii[ind-off[3]];
            ii[ind]=I_upw*coeff[_i_nu][0]+coeff[_i_nu][1];
          }
        }
      }
    } else if(iz) {
      for(int z = 0; z < nz-1; z++) {
        int zi = (zi_i-zl) + z*zstep;
//#pragma acc parallel loop gang async \
// present(this[:1], ii[:nx*ny*nz], coeff[:nx*ny*nz][:2], stride[:2], str_inu[:2], off[:4], c[:4])
        for(int x = 0; x < nx-1; x++) {
          int xi = (xi_i-xl) + x*xstep;
//#pragma acc loop vector
	  for(int y = 0; y < ny-1; y++) {
            int yi = (yi_i-yl) + y*ystep;
            int ind = yi*stride[0] + xi*stride[1] + zi;
            int _i_nu = y*str_inu[0] + x*str_inu[1] + z;
            double I_upw=c[0]*ii[ind-off[0]]+c[1]*ii[ind-off[1]]+c[2]*ii[ind-off[2]]+c[3]*ii[ind-off[3]];
            ii[ind]=I_upw*coeff[_i_nu][0]+coeff[_i_nu][1];
          }
        }
      }
    }
//#pragma acc wait
} // end data
//    ACCH::UpdateCPU(ii, nx*ny*nz*sizeof(double));
    PGI_COMPARE(I_n, double, nx*ny*nz, "I_n", "rt.cc", "RTS::driver", 22)
  }

  if(NDIM==2){
    for(int i=0;i<4;i++) off[i]=ixstep[i]*stride[1]+izstep[i];
    for(int xi=xi_i;xi!=xi_f+xstep;xi=xi+xstep) {
      int xoff=(xi-xl)*stride[1];
      for(int zi=zi_i;zi!=zi_f+zstep;zi=zi+zstep) {
        int ind=xoff+zi;
        double I_upw=c[0]*ii[ind-off[0]]+c[1]*ii[ind-off[1]]+c[2]*ii[ind-off[2]]+c[3]*ii[ind-off[3]];
        ii[ind]=I_upw*coeff[i_nu][0]+coeff[i_nu][1];
        i_nu+=1;
      }
    }
  }

  if(NDIM==1) {
    for(int i=0;i<4;i++) off[i]=izstep[i];
    for(int zi=zi_i;zi!=zi_f+zstep;zi=zi+zstep) {
      double I_upw=c[0]*ii[zi-off[0]]+c[1]*ii[zi-off[1]]+c[2]*ii[zi-off[2]]+c[3]*ii[zi-off[3]];
      ii[zi]=I_upw*coeff[i_nu][0]+coeff[i_nu][1];
      i_nu+=1;
    }
  }
}

void RTS::driver(double DZ, double DX, double DY, int band){
  double etime=0.0,atime=0.0,cmp_time1=0.0,cmp_time2=0.0,buf_time=0.0,err_time=0.0,flx_time=0.0,tau_time=0.0; 
  double ttime=MPI_Wtime();
  
  const int stride[2]={nx*nz,nz}; 
  int stepvec[3][4][3] = { {{1,0,0},{1,0,1},{1,1,0},{1,1,1}},
               {{0,1,0},{1,1,0},{0,1,1},{1,1,1}},
               {{0,0,1},{0,1,1},{1,0,1},{1,1,1}} };
  
  memset(I_n,0,nx*ny*nz*sizeof(double));
  memset(Fz,0,nx*ny*nz*sizeof(double));
  memset(Fx,0,nx*ny*nz*sizeof(double));
  memset(Fy,0,nx*ny*nz*sizeof(double));
  memset(J_band,0,nx*ny*nz*sizeof(double));

  PGI_COMPARE(I_n, double, nx*ny*nz, "I_n", "rt.cc", "RTS::driver", 500);

  double maxerr_up=0.0,maxerr_down=0.0;
// loop over octants & determination of loop direction
  double itavg=0.0;
  double aravg=0.0;

// main loop
  for(int ZDIR=UP;ZDIR<=DOWN;++ZDIR){
    int zi_i=(ZDIR==UP)?zl+1:zh-1,zi_f=(ZDIR==UP)?zh:zl,zstep=(ZDIR==UP)?1:-1;
    for(int XDIR=RIGHT;XDIR<=LEFT;++XDIR){
      int xi_i=(XDIR==RIGHT)?xl+1:xh-1,xi_f=(XDIR==RIGHT)?xh:xl,xstep=(XDIR==RIGHT)?1:-1;
      for(int YDIR=FWD;YDIR<=BWD;++YDIR){
        int yi_i=(YDIR==FWD)?yl+1:yh-1,yi_f=(YDIR==FWD)?yh:yl,ystep=(YDIR==FWD)?1:-1;
          for(int l=0;l<NMU;++l){
            double I_min=max(1.0,threshold*Fr_mean[band]/(NMU*pow(2,NDIM)));
            double c[]={a_00[ibase[l]][l],a_01[ibase[l]][l],a_10[ibase[l]][l],a_11[ibase[l]][l]};
        for(int m=0;m<=3;++m){
          ixstep[m]=stepvec[ibase[l]][m][0]*xstep;
          iystep[m]=stepvec[ibase[l]][m][1]*ystep;
          izstep[m]=stepvec[ibase[l]][m][2]*zstep;
        } 
        double stime=MPI_Wtime();
        interpol(zi_i,zi_f,zstep,xi_i,xi_f,xstep,yi_i,yi_f,ystep,l,coeff,B);
        cmp_time1+=MPI_Wtime()-stime;
        stime=MPI_Wtime(); 
            int rt_iter=0;
            int rt_min_iter=(call_count<2)?0:numits[band][YDIR][XDIR][ZDIR][l]-(!(call_count%3));            
        double gmaxerr=1.0E10;
// Iteration start
        while(gmaxerr>=threshold){
          rt_iter+=1;
          itavg+=1.0;
// read BC from recvbuf
          stime = MPI_Wtime();
          readbuf(band,l,ZDIR,XDIR,YDIR);
          buf_time += MPI_Wtime()-stime;
// loop over the grid points
          stime=MPI_Wtime();
          integrate(coeff, c, stride, ystep, xstep, zstep, yi_i, yi_f, xi_i, xi_f, zi_i, zi_f);
          cmp_time2 += MPI_Wtime()-stime;
          // write new BC, store old BC in oldbuf
          stime = MPI_Wtime();
          writebuf(band,l,ZDIR,XDIR,YDIR); 
          buf_time += MPI_Wtime()-stime;
          stime = MPI_Wtime();
          exchange(band, l, ZDIR, XDIR,YDIR);    
          etime += MPI_Wtime()-stime;
          
          if(rt_iter>=rt_min_iter){
            stime = MPI_Wtime();
            double err=error(band,l,ZDIR,XDIR,YDIR,I_min);
            err_time += MPI_Wtime()-stime;
            stime = MPI_Wtime();
            MPI_Allreduce(&err,&gmaxerr,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
            aravg += 1.0;
            atime += MPI_Wtime()-stime;
          }
        } // end main loop
        if(ZDIR==UP){
          maxerr_up=max(maxerr_up,gmaxerr);
        }else{
          maxerr_down=max(maxerr_down,gmaxerr);
        }
        numits[band][YDIR][XDIR][ZDIR][l]=rt_iter;
// read BC from recvbuf 
            stime=MPI_Wtime();
            readbuf(band,l,ZDIR,XDIR,YDIR);
            buf_time += MPI_Wtime()-stime;
// increment for F and J 
            stime=MPI_Wtime();
            flux(l,ZDIR,XDIR,YDIR);
            flx_time += MPI_Wtime()-stime;
      }// end l
      }// end YDIR
    }// end XDIR
  }// end ZDIR 

// continuum band last -> can use tau in tau_slice as for grey RT 
  double stime=MPI_Wtime();
  tauscale_qrad(band,DX,DY,DZ,B);
  tau_time+=MPI_Wtime()-stime;

  if((myrank==0) && (verbose>1)){
    fprintf(stdout,"rt_driver iter : %f %f \n",aravg/(8.0*NMU),itavg/(8.0*NMU));
    fprintf(stdout,"rt_driver error: %e %e \n",maxerr_up,maxerr_down);
    /*
    int oct=0;
    for(int i1=0;i1<=1;i1++)
      for(int i2=0;i2<=1;i2++)
	for(int i3=0;i3<=1;i3++){
	  oct +=1;
	  for(int l=0;l<NMU;++l)
	    fprintf(stdout,"rt_driver iter : %d %d %d \n",oct,l,numits[0][i1][i2][i3][l]);
	}
    */
  }

  ttime=MPI_Wtime()-ttime;  
  if((myrank==0) && (verbose>2))
    fprintf(stdout,"rt_driver time : %f %f %f %f %f %f %f %f %f %f \n",ttime,cmp_time1,
       cmp_time2,buf_time,err_time,flx_time,tau_time,etime,atime,(etime+atime)/ttime);
   

  call_count+=1;
}

void RTS::interpol(int zi_i,int zi_f,int zstep,int xi_i,int xi_f,int xstep,
           int yi_i,int yi_f,int ystep,int l,double** coeff, double * Ss)
{

  double ds3=ds_upw[l]*inv3,ds6=ds_upw[l]*inv6;
  double c[]={a_00[ibase[l]][l],a_01[ibase[l]][l],a_10[ibase[l]][l],a_11[ibase[l]][l]};

  int zmin=(zi_i<zi_f)?zi_i:zi_f;
  int zmax=(zi_i>zi_f)?zi_i:zi_f;

  double r_upw[zmax+1],k_upw[zmax+1],S_upw[zmax+1],r0[zmax+1],k0[zmax+1],S0[zmax+1];
  double _r_upw, _k_upw, _S_upw, _r0, _k0, _S0;
  int i_nu=0;

  const int stride[] = {nx*nz, nz, 1};
  const int stride_inter[] = {(nx-1)*(nz-1), nz-1, 1};
  const int off[] = {
    iystep[0]*stride[0]+ixstep[0]*stride[1]+izstep[0],
    iystep[1]*stride[0]+ixstep[1]*stride[1]+izstep[1],
    iystep[2]*stride[0]+ixstep[2]*stride[1]+izstep[2],
    iystep[3]*stride[0]+ixstep[3]*stride[1]+izstep[3]
  };

  if(NDIM==3){
    ACCH::UpdateGPU(rho, nx*ny*nz*sizeof(double));
    ACCH::UpdateGPU(kap, nx*ny*nz*sizeof(double));
    ACCH::UpdateGPU(Ss, nx*ny*nz*sizeof(double));
#pragma acc parallel loop gang collapse(2) \
 private(_r_upw, _k_upw, _S_upw, _r0, _k0, _S0) \
 present(this[:1], rho[:nx*ny*nz], kap[:nx*ny*nz], Ss[:nx*ny*nz], coeff[:nx*ny*nz][:2])
    for(int y = 0; y < ny-1; y++) {
      for(int x = 0; x < nx-1; x++) {
        int yi = (yi_i-yl) + y*ystep;
        int xi = (xi_i-xl) + x*xstep;
#pragma acc loop vector \
 private(_r_upw, _k_upw, _S_upw, _r0, _k0, _S0)
        for(int z = 0; z < nz-1; z++) {
          int zi = (zi_i-zl) + z*zstep;
          int ind = yi*stride[0] + xi*stride[1] + zi;
	  int i_nu_ = y*stride_inter[0] + x*stride_inter[1] + z;
          _r_upw=
            c[0]*rho[ind-off[0]]+
            c[1]*rho[ind-off[1]]+
            c[2]*rho[ind-off[2]]+
            c[3]*rho[ind-off[3]];
          
          _k_upw=
            c[0]*kap[ind-off[0]]+
            c[1]*kap[ind-off[1]]+
            c[2]*kap[ind-off[2]]+
            c[3]*kap[ind-off[3]];

          _S_upw=
            c[0]*Ss[ind-off[0]]+
            c[1]*Ss[ind-off[1]]+
            c[2]*Ss[ind-off[2]]+
            c[3]*Ss[ind-off[3]];

          _r0=rho[ind];
          _k0=kap[ind];
          _S0=Ss[ind];

          double dt=ds3*(_k_upw*_r_upw+_k0*_r0)+ds6*(_k0*_r_upw+_k_upw*_r0);
          double expo=exp(-dt);
          double w0,w1;
          if (dt > dtau_min){
            w0=1.0-expo;
            w1=w0-dt*expo;
          }else{
            w0=dt-dt*dt/2.0+dt*dt*dt/6.0;
            w1=dt*dt/2.0-dt*dt*dt/3.0;
          }
          double source=_S0*(w0-w1/dt)+_S_upw*(w1/dt);

          if (dt > dtau_min2){
            coeff[i_nu_][0] = expo;
            coeff[i_nu_][1] = source;
          }else{
            coeff[i_nu_][0] = 1.0; 
            coeff[i_nu_][1] = 0.0;
          }
        }
      }
    }
    ACCH::UpdateCPU2D<double>(coeff, nx*ny*nz, 2);
    //PGI_COMPARE(&coeff[0][0], double, (nx-1)*(ny-1)*(nz-1)*2, "coeff", "rt.cc", "RTS::interpol", 23)

  }

  if(NDIM==1){
    for(int zi=zmin;zi<=zmax;++zi){
      r_upw[zi]=
        c[0]*rho[(zi-izstep[0]-zl)]+
        c[1]*rho[(zi-izstep[1]-zl)]+
        c[2]*rho[(zi-izstep[2]-zl)]+
        c[3]*rho[(zi-izstep[3]-zl)];

      k_upw[zi]=
        c[0]*kap[zi-izstep[0]-zl]+
        c[1]*kap[zi-izstep[1]-zl]+
        c[2]*kap[zi-izstep[2]-zl]+
        c[3]*kap[zi-izstep[3]-zl];

      S_upw[zi]=
        c[0]*Ss[zi-izstep[0]-zl]+
        c[1]*Ss[zi-izstep[1]-zl]+
        c[2]*Ss[zi-izstep[2]-zl]+
        c[3]*Ss[zi-izstep[3]-zl];

      r0[zi]=rho[zi-zl];
      k0[zi]=kap[zi-zl];
      S0[zi]=Ss[zi-zl];
    }

    for(int zi=zi_i;zi!=zi_f+zstep;zi=zi+zstep){
      double dt=ds3*(k_upw[zi]*r_upw[zi]+k0[zi]*r0[zi])+ds6*(k0[zi]*r_upw[zi]+k_upw[zi]*r0[zi]);
      double expo=exp(-dt);
      double w0,w1;
      if (dt > dtau_min){
        w0=1.0-expo;
        w1=w0-dt*expo;
      }else{
        w0=dt-dt*dt/2.0+dt*dt*dt/6.0;
        w1=dt*dt/2.0-dt*dt*dt/3.0;
      }
      double source=S0[zi]*(w0-w1/dt)+S_upw[zi]*(w1/dt);

      if (dt > dtau_min2){
        coeff[i_nu][0] = expo;
        coeff[i_nu][1] = source;
      }else{
        coeff[i_nu][0] = 1.0; 
        coeff[i_nu][1] = 0.0;
      }
    i_nu+=1;
    }
  }

  if(NDIM==2){
    for(int xi=xi_i;xi!=xi_f+xstep;xi=xi+xstep){
      for(int zi=zmin;zi<=zmax;++zi){
        r_upw[zi]=
          c[0]*rho[(xi-ixstep[0]-xl)*nz+(zi-izstep[0]-zl)]+
          c[1]*rho[(xi-ixstep[1]-xl)*nz+(zi-izstep[1]-zl)]+
          c[2]*rho[(xi-ixstep[2]-xl)*nz+(zi-izstep[2]-zl)]+
          c[3]*rho[(xi-ixstep[3]-xl)*nz+(zi-izstep[3]-zl)];

        k_upw[zi]=
          c[0]*kap[(xi-ixstep[0]-xl)*nz+(zi-izstep[0]-zl)]+
          c[1]*kap[(xi-ixstep[1]-xl)*nz+(zi-izstep[1]-zl)]+
          c[2]*kap[(xi-ixstep[2]-xl)*nz+(zi-izstep[2]-zl)]+
          c[3]*kap[(xi-ixstep[3]-xl)*nz+(zi-izstep[3]-zl)];

        S_upw[zi]=
          c[0]*Ss[(xi-ixstep[0]-xl)*nz+(zi-izstep[0]-zl)]+
          c[1]*Ss[(xi-ixstep[1]-xl)*nz+(zi-izstep[1]-zl)]+
          c[2]*Ss[(xi-ixstep[2]-xl)*nz+(zi-izstep[2]-zl)]+
          c[3]*Ss[(xi-ixstep[3]-xl)*nz+(zi-izstep[3]-zl)];

        r0[zi]=rho[(xi-xl)*nz+(zi-zl)];
        k0[zi]=kap[(xi-xl)*nz+(zi-zl)];
        S0[zi]=Ss[(xi-xl)*nz+(zi-zl)];
      }
      for(int zi=zi_i;zi!=zi_f+zstep;zi=zi+zstep){
        double dt=ds3*(k_upw[zi]*r_upw[zi]+k0[zi]*r0[zi])+ds6*(k0[zi]*r_upw[zi]+k_upw[zi]*r0[zi]);
        double expo=exp(-dt);
        double w0,w1;
        if (dt > dtau_min){
          w0=1.0-expo;
          w1=w0-dt*expo;
        }else{
          w0=dt-dt*dt/2.0+dt*dt*dt/6.0;
          w1=dt*dt/2.0-dt*dt*dt/3.0;
        }
        double source=S0[zi]*(w0-w1/dt)+S_upw[zi]*(w1/dt);

        if (dt > dtau_min2){
          coeff[i_nu][0] = expo;
          coeff[i_nu][1] = source;
        }else{
          coeff[i_nu][0] = 1.0; 
          coeff[i_nu][1] = 0.0;
        }
      i_nu+=1;
      }
    }
  }
}

double RTS::error(int band,int l,int ZDIR,int XDIR,int YDIR,double I_min)
{
  int z0=(ZDIR==UP),z1=(ZDIR==DOWN); 
  int x0=(XDIR==RIGHT),x1=(XDIR==LEFT); 
  double err_max=0.0;
//
  if (NDIM==3){
    real *ysb=y_sbuf[band][YDIR][XDIR][ZDIR][l],*yob=y_oldbuf[band][YDIR][XDIR][ZDIR][l];
    PGI_COMPARE(ysb, double, nx*nz, "ysb", "rt.cc", "RTS::error", 100)
    PGI_COMPARE(yob, double, nx*nz, "yob", "rt.cc", "RTS::error", 101)
    ACCH::UpdateGPU(ysb, nx*nz*sizeof(real));
    ACCH::UpdateGPU(yob, nx*nz*sizeof(real));
#pragma acc parallel loop collapse(2) reduction(max:err_max) \
 present(this[:1], ysb[:nx*nz], yob[:nx*nz])
    for(int z=z0;z<nz-z1;++z)
      for(int x=x0;x<nx-x1;++x){
        int ind=z*nx+x;
        err_max=max(err_max,fabs(((double) (ysb[ind]-yob[ind]))/max(I_min,(double) yob[ind])));
     }
  }
//
  if (NDIM>1){
  real *xsb=x_sbuf[band][YDIR][XDIR][ZDIR][l],*xob=x_oldbuf[band][YDIR][XDIR][ZDIR][l];
    PGI_COMPARE(xsb, double, ny*nz, "xsb", "rt.cc", "RTS::error", 102)
    PGI_COMPARE(xob, double, ny*nz, "xob", "rt.cc", "RTS::error", 103)
    ACCH::UpdateGPU(xsb, ny*nz*sizeof(real));
    ACCH::UpdateGPU(xob, ny*nz*sizeof(real));
#pragma acc parallel loop collapse(2) reduction(max:err_max) \
 present(this[:1], xsb[:ny*nz], xob[:ny*nz])
  for(int z=z0;z<nz-z1;++z)
    for(int y=0;y<ny;++y){
      int ind=z*ny+y;
      err_max=max(err_max,fabs(((double) (xsb[ind]-xob[ind]))/max( I_min,(double)xob[ind])));
    }
  }
//
  real *zsb=z_sbuf[band][YDIR][XDIR][ZDIR][l],*zob=z_oldbuf[band][YDIR][XDIR][ZDIR][l];
  PGI_COMPARE(zsb, double, nx*ny, "zsb", "rt.cc", "RTS::error", 104)
  PGI_COMPARE(zob, double, nx*ny, "zob", "rt.cc", "RTS::error", 105)
  ACCH::UpdateGPU(zsb, ny*nx*sizeof(real));
  ACCH::UpdateGPU(zob, ny*nx*sizeof(real));
#pragma acc parallel loop reduction(max:err_max) \
 present(this[:1], zsb[:ny*nx], zob[:ny*nx])
  for(int ind=0;ind<nx*ny;++ind){
    err_max=max(err_max,fabs(((double) (zsb[ind]-zob[ind]))/max(I_min,(double)zob[ind])));
  }
//

  PGI_COMPARE(&err_max, double, 1, "err_max", "rt.cc", "RTS::error", 24)
  return err_max;
}

void RTS::readbuf(int band,int l,int ZDIR,int XDIR,int YDIR)
{//RHC
  ACCH::UpdateGPU(I_n, nx*ny*nz*sizeof(double));  
  if(NDIM==3){
    int y0=(YDIR==FWD)?0:ny-1;
    int y = y0*nx*nz;
    real * ysb = y_sbuf[band][YDIR][XDIR][ZDIR][l];
    real * yrb = y_rbuf[band][YDIR][XDIR][ZDIR][l];
    real * yob = y_oldbuf[band][YDIR][XDIR][ZDIR][l];
    ACCH::UpdateGPU(ysb, nx*nz*sizeof(real));
    ACCH::UpdateGPU(yrb, nx*nz*sizeof(real));
    ACCH::UpdateGPU(yob, nx*nz*sizeof(real));
#pragma acc parallel loop collapse(2) \
 present(this[:1], I_n[:nx*ny*nz], yrb[:nx*nz])
    for(int x=0;x<nx;++x)
      for(int z=0;z<nz;++z)
        I_n[y+x*nz+z]= (double) yrb[z*nx+x];
#pragma acc parallel loop \
 present(this[:1], yob[:nx*nz], ysb[:nx*nz])
    for(int i = 0; i < nx*nz; i++)
      yob[i] = ysb[i];
    ACCH::UpdateCPU(yob, nx*nz*sizeof(real));
  }

  if(NDIM>1){
    int x0=(XDIR==RIGHT)?0:nx-1;
    int x = x0*nz;
    real * xsb = x_sbuf[band][YDIR][XDIR][ZDIR][l];
    real * xrb = x_rbuf[band][YDIR][XDIR][ZDIR][l];
    real * xob = x_oldbuf[band][YDIR][XDIR][ZDIR][l];
    ACCH::UpdateGPU(xsb, ny*nz*sizeof(real));
    ACCH::UpdateGPU(xrb, ny*nz*sizeof(real));
    ACCH::UpdateGPU(xob, ny*nz*sizeof(real));
#pragma acc parallel loop collapse(2) \
 present(this[:1], I_n[:nx*ny*nz], xrb[:ny*nz])
    for(int y=0;y<ny;++y)
      for(int z=0;z<nz;++z)
        I_n[y*nx*nz+x+z]= (double) xrb[z*ny+y];
#pragma acc parallel loop \
 present(this[:1], xob[:ny*nz], xsb[:ny*nz])
    for(int i = 0; i < ny*nz; i++)
      xob[i] = xsb[i];
    ACCH::UpdateCPU(xob, ny*nz*sizeof(real));
  }
  
  int z0=(ZDIR==UP)?0:nz-1;
  real * zsb = z_sbuf[band][YDIR][XDIR][ZDIR][l];
  real * zrb = z_rbuf[band][YDIR][XDIR][ZDIR][l];
  real * zob = z_oldbuf[band][YDIR][XDIR][ZDIR][l];
  ACCH::UpdateGPU(zsb, ny*nx*sizeof(real));
  ACCH::UpdateGPU(zrb, ny*nx*sizeof(real));
  ACCH::UpdateGPU(zob, ny*nx*sizeof(real));
#pragma acc parallel loop collapse(2) \
 present(this[:1], I_n[:nx*ny*nz], zrb[:ny*nx])
  for(int y=0;y<ny;++y)
    for(int x=0;x<nx;++x)
      I_n[y*nx*nz+x*nz+z0]= (double) zrb[x*ny+y];
#pragma acc parallel loop \
 present(this[:1], zob[:ny*nx], zsb[:ny*nx])
  for(int i = 0; i < nx*ny; i++)
    zob[i] = zsb[i];
  ACCH::UpdateCPU(zob, ny*nx*sizeof(real));
  ACCH::UpdateCPU(I_n, nx*ny*nz*sizeof(double)); 
}

void RTS::writebuf(int band, int l,int ZDIR,int XDIR,int YDIR){
  ACCH::UpdateGPU(I_n, nx*ny*nz*sizeof(double));
  if (NDIM==3){
    int y0=(YDIR==FWD)?ny-1:0;
    int y = y0*nx*nz;
    real * ysb = y_sbuf[band][YDIR][XDIR][ZDIR][l];
    ACCH::UpdateGPU(ysb, nx*nz*sizeof(real));
#pragma acc parallel loop collapse(2) \
 present(this[:1], I_n[:nx*ny*nz], ysb[:nx*nz])
    for(int x=0;x<nx;++x)
      for(int z=0;z<nz;++z)
        ysb[z*nx+x]=(real) I_n[y+x*nz+z];
    ACCH::UpdateCPU(ysb, nx*nz*sizeof(real));
  }

  if (NDIM>1){
    int x0=(XDIR==RIGHT)?nx-1:0;
    int x = x0*nz;
    real * xsb = x_sbuf[band][YDIR][XDIR][ZDIR][l];
    ACCH::UpdateGPU(xsb, ny*nz*sizeof(real));
#pragma acc parallel loop collapse(2) \
 present(this[:1], I_n[:nx*ny*nz], xsb[:ny*nz])
    for(int y=0;y<ny;++y)
      for(int z=0;z<nz;++z)
        xsb[z*ny+y]=(real) I_n[y*nx*nz+x+z];
    ACCH::UpdateCPU(xsb, ny*nz*sizeof(real));
  }

  int z0=(ZDIR==UP)?nz-1:0;
  real * zsb = z_sbuf[band][YDIR][XDIR][ZDIR][l];
  ACCH::UpdateGPU(zsb, nx*ny*sizeof(real));
#pragma acc parallel loop collapse(2) \
 present(this[:1], I_n[:nx*ny*nz], zsb[:ny*nx])
  for(int y=0;y<ny;++y)
    for(int x=0;x<nx;++x)
      zsb[x*ny+y]=(real) I_n[y*nx*nz+x*nz+z0];
  ACCH::UpdateCPU(zsb, nx*ny*sizeof(real));

}

void RTS::exchange(int band,int l,int ZDIR,int XDIR,int YDIR)
{
  MPI_Status s1[2],s2[2],s3[2];
  MPI_Request r1[2],r2[2],r3[2];
  int tag1=1,tag2=2,tag3=3;

  // Update v. z_sbuf an den Endpkten, damit der Eckpkt-Wert vom diagonalen Nachbarn uebertragen wird
  int x0=0,y0=0,z0=0;
  int dest_rk;
  int source_rk;

  real * ysb = y_sbuf[band][YDIR][XDIR][ZDIR][l];
  real * yrb = y_rbuf[band][YDIR][XDIR][ZDIR][l];
  real * xsb = x_sbuf[band][YDIR][XDIR][ZDIR][l];
  real * xrb = x_rbuf[band][YDIR][XDIR][ZDIR][l];
  real * zsb = z_sbuf[band][YDIR][XDIR][ZDIR][l];
  real * zrb = z_rbuf[band][YDIR][XDIR][ZDIR][l];

  if (NDIM==3){
    // y-direction
    dest_rk=(YDIR==FWD)?rightr[2]:leftr[2];
    source_rk=(YDIR==FWD)?leftr[2]:rightr[2];
    MPI_Irecv(yrb,nx*nz,REALTYPE,source_rk,tag2,MPI_COMM_WORLD,r2+0);
    MPI_Isend(ysb,nx*nz,REALTYPE,dest_rk,tag2,MPI_COMM_WORLD,r2+1);
    MPI_Waitall(2,r2,s2);
  
    z0=(ZDIR==UP)?nz-1:0;
    y0=(YDIR==FWD)?0:ny-1;
    ACCH::UpdateGPU(zsb, nx*ny*sizeof(real));
    ACCH::UpdateGPU(yrb, nx*nz*sizeof(real));
#pragma acc parallel loop present(zsb[:nx*ny], yrb[:nx*nz])
    for(int x=0;x<nx;++x)
      zsb[x*ny+y0]=yrb[z0*nx+x];
    ACCH::UpdateCPU(zsb, nx*ny*sizeof(real));
  
    x0=(XDIR==RIGHT)?nx-1:0;
    y0=(YDIR==FWD)?0:ny-1;
    ACCH::UpdateGPU(xsb, ny*nz*sizeof(real));
#pragma acc parallel loop present(xsb[:ny*nz], yrb[:nx*nz])
    for(int z=0;z<nz;++z)
      xsb[z*ny+y0]=yrb[z*nx+x0];
    ACCH::UpdateCPU(xsb, ny*nz*sizeof(real));
  }

  if (NDIM >1){
    // x-direction
    dest_rk=(XDIR==RIGHT)?rightr[1]:leftr[1];
    source_rk=(XDIR==RIGHT)?leftr[1]:rightr[1];
    MPI_Irecv(xrb,nz*ny,REALTYPE,source_rk,tag1,MPI_COMM_WORLD,r1+0);
    MPI_Isend(xsb,nz*ny,REALTYPE,dest_rk,tag1,MPI_COMM_WORLD,r1+1);

    MPI_Waitall(2,r1,s1);
 
    x0=(XDIR==RIGHT)?0:nx-1;
    z0=(ZDIR==UP)?nz-1:0;
    ACCH::UpdateGPU(zsb, nx*ny*sizeof(real));
    ACCH::UpdateGPU(xrb, ny*nz*sizeof(real));
#pragma acc parallel loop present(zsb[:nx*ny], xrb[:ny*nz])
    for(int y=0;y<ny;++y)
      zsb[x0*ny+y]=xrb[z0*ny+y];
    ACCH::UpdateCPU(zsb, ny*nx*sizeof(real));
  }
  
// z-direction
  dest_rk  =(ZDIR==UP)?rightr[0]:leftr[0];
  source_rk=(ZDIR==UP)?leftr[0]:rightr[0];
  MPI_Irecv(z_rbuf[band][YDIR][XDIR][ZDIR][l],nx*ny,REALTYPE,source_rk,tag3,MPI_COMM_WORLD,r3+0);
  MPI_Isend(z_sbuf[band][YDIR][XDIR][ZDIR][l],nx*ny,REALTYPE,dest_rk,tag3,MPI_COMM_WORLD,r3+1);
  MPI_Waitall(2,r3,s3);

}

void RTS::flux(int l,int ZDIR,int XDIR,int YDIR)
{
  double zsign=(ZDIR==UP)?1.0:-1.0,xsign=(XDIR==RIGHT)?1.0:-1.0,ysign=(YDIR==FWD)?1.0:-1.0;
  double c_J = 0.125*wmu[l];
  double c_z = 0.5*PI*zsign*wmu[l]*xmu[2][l];
  double c_x = 0.5*PI*xsign*wmu[l]*xmu[0][l];
  double c_y = 0.5*PI*ysign*wmu[l]*xmu[1][l];
  ACCH::UpdateGPU(I_n, nx*ny*nz*sizeof(double));
  ACCH::UpdateGPU(J_band, nx*ny*nz*sizeof(double));
  ACCH::UpdateGPU(Fz, nx*ny*nz*sizeof(double));
  ACCH::UpdateGPU(Fx, nx*ny*nz*sizeof(double));
  ACCH::UpdateGPU(Fy, nx*ny*nz*sizeof(double));
#pragma acc parallel loop collapse(2) gang \
 present(this[:1], I_n[:nx*ny*nz], J_band[:nx*ny*nz], \
         Fz[:nx*ny*nz], Fy[:nx*ny*nz], Fx[:nx*ny*nz])
  for(int y = 0; y < ny; y++) {
    for(int x = 0; x < nx; x++) {
      int yxind = y*nx*nz + x*nz;
#pragma acc loop vector
      for(int z = 0; z < nz; z++) {
        int ind = yxind + z;
        double tmp = I_n[ind];
        J_band[ind] += c_J*tmp;
        Fz[ind]     += c_z*tmp;
        Fx[ind]     += c_x*tmp;
        Fy[ind]     += c_y*tmp;
      }
    }
  }
  ACCH::UpdateCPU(J_band, nx*ny*nz*sizeof(double));
  ACCH::UpdateCPU(Fz, nx*ny*nz*sizeof(double));
  ACCH::UpdateCPU(Fx, nx*ny*nz*sizeof(double));
  ACCH::UpdateCPU(Fy, nx*ny*nz*sizeof(double));
  PGI_COMPARE(J_band, double, nx*ny*nz, "J_band", "rt.cc", "RTS::flux", 29)
  PGI_COMPARE(Fx, double, nx*ny*nz, "Fx", "rt.cc", "RTS::flux", 30)
  PGI_COMPARE(Fy, double, nx*ny*nz, "Fy", "rt.cc", "RTS::flux", 31)
  PGI_COMPARE(Fz, double, nx*ny*nz, "Fz", "rt.cc", "RTS::flux", 32)
}

void RTS::get_Tau_and_Iout(GridData &Grid, const RunData &Run, const PhysicsData &Physics, double DZ, float * B_Iout_tab, float ** kap_Iout_tab, double * I_band, int calc_int){
  
  double ttime=MPI_Wtime(),stime=0.0,atime=0.0;
  
  double ** sbuf = ACCH::Malloc2D<double>(ny, nx);
  double ** rbuf = ACCH::Malloc2D<double>(ny, nx);

  double N = pow(2,Grid.NDIM);

  const double Temp_TR = Physics.rt[i_rt_tr_tem];
  const double Pres_TR = Physics.rt[i_rt_tr_pre];
  ACCH::UpdateGPU(Grid.pres, Grid.bufsize*sizeof(double));
  ACCH::UpdateGPU(Grid.temp, Grid.bufsize*sizeof(double));
  ACCH::UpdateGPU(Grid.U, Grid.bufsize*sizeof(cState));
  ACCH::UpdateGPU(rho, nx*ny*nz*sizeof(double));
  ACCH::UpdateGPU(tab_T, NT*sizeof(double));
  ACCH::UpdateGPU(tab_p, Np*sizeof(double));
  ACCH::UpdateGPU(B, nx*ny*nz*sizeof(double));
  ACCH::UpdateGPU(kap, nx*ny*nz*sizeof(double));
  ACCH::UpdateGPU(Tau, nx*ny*nz*sizeof(double));
  ACCH::UpdateGPU(B_Iout_tab, NT*sizeof(float));
  ACCH::UpdateGPU2D<float>(kap_Iout_tab, NT, Np);
  ACCH::UpdateGPU(I_band, ny*nx*sizeof(double));
  ACCH::UpdateGPU(invT_tab, NT*sizeof(double));
  ACCH::UpdateGPU(invP_tab, Np*sizeof(double));

#pragma acc parallel loop gang collapse(2) \
 present(this[:1], Grid[:1], Grid.pres[:Grid.bufsize], Grid.temp[:Grid.bufsize], Grid.U[:Grid.bufsize], \
         rho[:nx*ny*nz], tab_T[:NT], tab_p[:Np], B_Iout_tab[:NT], kap_Iout_tab[:NT][:Np], \
	 kap[:nx*ny*nz], B[:nx*ny*nz], invT_tab[:NT], invP_tab[:Np]) \
 copyin(next[1:2])
  for(int y = 0; y < ny; y++) {
    for(int x = 0; x < nx; x++) {
      int off0 = (x+xl)*next[1]+(y+yl)*next[2];
      int off1 = (x+xl)*next[1]+(y+yo+yl)*next[2];
      int off2 = (x+xo+xl)*next[1]+(y+yl)*next[2];
      int off3 = (x+xo+yl)*next[1]+(y+yo+yl)*next[2];
#pragma acc loop vector
      for(int z = 0; z < nz; z++) {
        int inode[]={off0+z+zl,off0+z+zo+zl,off2+z+zl,off2+z+zo+zl,off1+z+zl,off1+z+zo+zl,off3+z+zl,off3+z+zo+zl};
        double lgP = 0, lgT = 0, rm = 0;
#pragma acc loop seq
        for(int l=0;l<N;++l){
          lgP+=Grid.pres[inode[l]];
          lgT+=Grid.temp[inode[l]];
          rm +=Grid.U[inode[l]].d;
        }
        lgP /= N;
        lgT /= N;
        rm  /= N;

        //disbale RT if Temp > Temp_TR above the photosphere
        double pswitch  = min(max(lgP-Pres_TR,0.0),1.0);
        double tswitch  = min(max(Temp_TR-lgT,0.0),1.0);
        double trswitch = max(pswitch,tswitch);

        lgP          = log(lgP);
        lgT          = log(lgT);
        rho[y*nx*nz+x*nz+z] = rm;

        lgT = min(max(lgT, (double) tab_T[0]),(double) tab_T[NT-1]);
        lgP = min(max(lgP, (double) tab_p[0]),(double) tab_p[Np-1]);

        int l=0;
        int m=0;
        if(lgT<tab_T[0])
          l=0;
        else if(lgT>tab_T[NT-1])
          l=NT-2;
        else {
#pragma acc loop seq
          for (l=0; l<=NT-2; l++)
            if ((lgT >= tab_T[l]) && (lgT <= tab_T[l+1]))
              break;
        }

        if(lgP<tab_p[0])
          m=0;
        else if(lgP>tab_p[Np-1])
          m=Np-2;
        else {
#pragma acc loop seq
          for (m=0; m<=Np-2; m++)
            if ((lgP >= tab_p[m]) && (lgP <= tab_p[m+1]))
              break;
        }

        double xt = (lgT-tab_T[l])*invT_tab[l];
        double xp = (lgP-tab_p[m])*invP_tab[m];

        // Interpolate for kappa and B
        B[y*nx*nz+x*nz+z]=exp(xt*B_Iout_tab[l+1]+(1.-xt)*B_Iout_tab[l]);

        kap[y*nx*nz+x*nz+z] = exp(xt*(xp*kap_Iout_tab[l+1][m+1]+(1.-xp)*kap_Iout_tab[l+1][m])+
                           (1.-xt)*(xp*kap_Iout_tab[l][m+1]+(1.-xp)*kap_Iout_tab[l][m]));


        // apply tr_switch
        B[y*nx*nz+x*nz+z]   *= trswitch;
        kap[y*nx*nz+x*nz+z] *= trswitch;

      }
    }
  }

#pragma acc parallel loop collapse(2) \
 present(this[:1], Tau[:nx*ny*nz], rho[:nx*ny*nz], kap[:nx*ny*nz])
  for(int y = 0; y < ny; y++)
    for(int x = 0; x < nx; x++) {
      Tau[y*nx*nz+x*nz+(nz-1)]=1.0e-12;
#pragma acc loop seq
      for(int z = nz-2; z >= 0; z--) {
        double k0=kap[y*nx*nz+x*nz+z],r0=rho[y*nx*nz+x*nz+z],k_upw=kap[y*nx*nz+x*nz+(z+1)],r_upw=rho[y*nx*nz+x*nz+(z+1)];
        double delta_tau=DZ*((k0*r0+k_upw*r_upw)*inv3+(k0*r_upw+k_upw*r0)*inv6);
        Tau[y*nx*nz+x*nz+z]=Tau[y*nx*nz+x*nz+(z+1)]+delta_tau;
      }
    }

  ACCH::UpdateCPU(rho, nx*ny*nz*sizeof(double));
  ACCH::UpdateCPU(B, nx*ny*nz*sizeof(double));
  ACCH::UpdateCPU(kap, nx*ny*nz*sizeof(double));

#pragma acc parallel loop collapse(2) \
 present(this[:1], sbuf[:ny][:nx], rbuf[:ny][:nx], Tau[:nx*ny*nz])
  for(int y=0;y<ny;++y){ // loop over RT grid
    for(int x=0;x<nx;++x){
      rbuf[y][x]=0.e0;
      sbuf[y][x]=Tau[y*nx*nz+x*nz];
    }
  }
  double ctime=MPI_Wtime();
  ACCH::UpdateCPU2D<double>(sbuf, ny, nx);
  MPI_Scan(sbuf[0],rbuf[0],nx*ny,MPI_DOUBLE,MPI_SUM,comm_col[lrank[2]][lrank[1]]);
  ACCH::UpdateGPU2D<double>(rbuf, ny, nx);
  stime+=MPI_Wtime()-ctime;
#pragma acc parallel loop gang collapse(2) \
 present(this[:1], rbuf[:ny][:nx], Tau[:nx*ny*nz])
  for(int y=0;y<ny;++y){ // loop over RT grid
    for(int x=0;x<nx;++x){
      rbuf[y][x]-=Tau[y*nx*nz+x*nz];
      for(int z=0;z<nz;++z){
        Tau[y*nx*nz+x*nz+z]+=rbuf[y][x];
      }
    }
  }

  ACCH::UpdateCPU(Tau, nx*ny*nz*sizeof(double));

  if (calc_int){
  //  Outgoing Intensity at top (Long Characteristics)
#pragma acc parallel loop gang collapse(2) \
 present(this[:1], rbuf[:ny][:nx], sbuf[:ny][:nx], B[:nx*ny*nz], Tau[:nx*ny*nz])
  for(int y=0;y<ny;++y){ // loop over RT grid
    for(int x=0;x<nx;++x){
      rbuf[y][x]=0.0;
      sbuf[y][x]=0.0;
      double tmp = 0.0;
#pragma acc loop vector reduction(+:tmp)
      for(int z=1;z<nz;++z){
        double Ss1 = B[y*nx*nz+x*nz+z];
        double Ss2 = B[y*nx*nz+x*nz+z-1];
        double delta_tau=Tau[y*nx*nz+x*nz+z-1]-Tau[y*nx*nz+x*nz+z];
        if(delta_tau>dtau_min){
          double edt=exp(-delta_tau);
          double c1=(1.0-edt)/delta_tau;
          tmp+=(Ss1*(1.0-c1)+Ss2*(c1-edt))*exp(-Tau[y*nx*nz+x*nz+z]);
        }else{
          tmp+=0.5*delta_tau*(Ss1+Ss2);
        }
      }
      sbuf[y][x] = tmp;
    }
  }
  ctime=MPI_Wtime();
  ACCH::UpdateCPU2D<double>(sbuf, ny, nx);
  MPI_Allreduce(sbuf[0],rbuf[0],nx*ny,MPI_DOUBLE,MPI_SUM,comm_col[lrank[2]][lrank[1]]);
  ACCH::UpdateGPU2D<double>(rbuf, ny, nx);
  atime+=MPI_Wtime()-ctime;
  ACCH::UpdateGPU(I_band, nx*ny*sizeof(double));
#pragma acc parallel loop collapse(2) \
 present(this[:1], I_band[:ny*nx], rbuf[:ny][:nx])
  for(int y=0;y<ny;++y){
    for(int x=0;x<nx;++x){
      I_band[y*nx+x]+=rbuf[y][x];
    }
  }
  ACCH::UpdateCPU(I_band, nx*ny*sizeof(double));
  }


  ACCH::Free2D<double>(sbuf, ny, nx);
  ACCH::Free2D<double>(rbuf, ny, nx);
  //
  ttime=MPI_Wtime()-ttime;
  if((myrank==0)&&(verbose>2)) printf("tau5000 time: %f %f %f %f \n",ttime,stime,atime,(stime+atime)/ttime);
  
}    
void RTS::tauscale_qrad(int band, double DX,double DY,double DZ, double * Ss){

  double ttime=MPI_Wtime(),stime=0.0,atime=0.0;
  double idx=1.0/DX,idy=1.0/DY,idz=1.0/DZ;

  if(NDIM==1) {
    idx=0.;
    idy=0.;
  }
  if(NDIM==2) {
    idy=0.0;
  }

  double ** sbuf = ACCH::Malloc2D<double>(ny, nx);
  double ** rbuf = ACCH::Malloc2D<double>(ny, nx);
  double **** Qtemp = ACCH::Malloc4D<double>(ny, nx, nz, 2);

  ACCH::UpdateGPU(rho, nx*ny*nz*sizeof(double));
  ACCH::UpdateGPU(kap, nx*ny*nz*sizeof(double));
#pragma acc parallel loop collapse(2) \
 present(this[:1], Tau[:nx*ny*nz], kap[:nx*ny*nz], rho[:nx*ny*nz])  
  for(int y=0;y<ny;++y){ // loop over RT grid
    for(int x=0;x<nx;++x){
      Tau[y*nx*nz+x*nz+nz-1]=1.0e-12;
#pragma acc loop seq
      for(int z=nz-2;z>=0;--z){
        double k0=kap[y*nx*nz+x*nz+z],r0=rho[y*nx*nz+x*nz+z],k_upw=kap[y*nx*nz+x*nz+z+1],r_upw=rho[y*nx*nz+x*nz+z+1];
        double delta_tau=DZ*((k0*r0+k_upw*r_upw)*inv3+(k0*r_upw+k_upw*r0)*inv6);
        Tau[y*nx*nz+x*nz+z]=Tau[y*nx*nz+x*nz+z+1]+delta_tau;
      }
    }
  }

#pragma acc parallel loop collapse(2) \
 present(this[:1], Tau[:nx*ny*nz], sbuf[:ny][:nx], rbuf[:ny][:nx])
  for(int y=0;y<ny;++y){ // loop over RT grid
    for(int x=0;x<nx;++x){
      rbuf[y][x]=0.e0;
      sbuf[y][x]=Tau[y*nx*nz+x*nz];
    }
  }

  double ctime=MPI_Wtime(); 
  ACCH::UpdateCPU2D<double>(sbuf, ny, nx);
  ACCH::UpdateCPU2D<double>(rbuf, ny, nx);
  MPI_Scan(sbuf[0],rbuf[0],nx*ny,MPI_DOUBLE,MPI_SUM,comm_col[lrank[2]][lrank[1]]);
  ACCH::UpdateGPU2D<double>(rbuf, ny, nx);
  stime+=MPI_Wtime()-ctime;
#pragma acc parallel loop gang collapse(2) \
 present(this[:1], Tau[:nx*ny*nz], rbuf[:ny][:nx])
  for(int y=0;y<ny;++y){ // loop over RT grid
    for(int x=0;x<nx;++x){
      rbuf[y][x]-=Tau[y*nx*nz+x*nz];
#pragma acc loop vector
      for(int z=0;z<nz;++z){
        Tau[y*nx*nz+x*nz+z]+=rbuf[y][x];
      }
    }
  }
  ACCH::UpdateCPU(Tau, nx*ny*nz*sizeof(double));

  if (need_I){
    ACCH::UpdateGPU(Ss, nx*ny*nz*sizeof(double));
    //  Outgoing Intensity at top (Long Characteristics)
#pragma acc parallel loop gang collapse(2) \
 present(this[:1], Tau[:nx*ny*nz], Ss[:nx*ny*nz], \
         rbuf[:ny][:nx], sbuf[:ny][:nx])
    for(int y=0;y<ny;++y){ // loop over RT grid
      for(int x=0;x<nx;++x){
        rbuf[y][x]=0.0;
        sbuf[y][x]=0.0;
#pragma acc loop vector
        for(int z=1;z<nz;++z){
          double Ss1 = Ss[y*nx*nz+x*nz+z];
          double Ss2 = Ss[y*nx*nz+x*nz+z-1];
          double delta_tau=Tau[y*nx*nz+x*nz+z-1]-Tau[y*nx*nz+x*nz+z];
          if(delta_tau>dtau_min){
            double edt=exp(-delta_tau);
            double c1=(1.0-edt)/delta_tau;
            sbuf[y][x]+=(Ss1*(1.0-c1)+Ss2*(c1-edt))*exp(-Tau[y*nx*nz+x*nz+z]);
          }else{
            sbuf[y][x]+=0.5*delta_tau*(Ss1+Ss2);
          }  
        }          
      }
    }
    ctime=MPI_Wtime(); 
    ACCH::UpdateCPU2D<double>(sbuf, ny, nx);
    ACCH::UpdateCPU2D<double>(rbuf, ny, nx);
    MPI_Allreduce(sbuf[0],rbuf[0],nx*ny,MPI_DOUBLE,MPI_SUM,comm_col[lrank[2]][lrank[1]]);
    ACCH::UpdateGPU2D<double>(rbuf, ny, nx);
    ctime+=MPI_Wtime()-ctime;
#pragma acc parallel loop collapse(2) \
 present(this[:1], I_o[:nx*ny], rbuf[:ny][:nx])
    for(int y=0;y<ny;++y){
      for(int x=0;x<nx;++x){
        I_o[y*nx+x]+=rbuf[y][x];
      }
    }
    ACCH::UpdateCPU(I_o, nx*ny*sizeof(double));

  }
//  radiative energy imbalance
  ACCH::UpdateGPU(I_n, nx*ny*nz*sizeof(double));
  ACCH::UpdateGPU(kap, nx*ny*nz*sizeof(double));
  ACCH::UpdateGPU(rho, nx*ny*nz*sizeof(double));
  ACCH::UpdateGPU(J_band, nx*ny*nz*sizeof(double));
  ACCH::UpdateGPU(St, nx*ny*nz*sizeof(double));
  ACCH::UpdateGPU(Jt, nx*ny*nz*sizeof(double));
  ACCH::UpdateGPU(Ss, nx*ny*nz*sizeof(double));
#pragma acc parallel loop gang collapse(2) \
 present(this[:1], I_n[:nx*ny*nz], kap[:nx*ny*nz], rho[:nx*ny*nz], \
         J_band[:nx*ny*nz], St[:nx*ny*nz], Jt[:nx*ny*nz], Ss[:nx*ny*nz])
  for(int y=0;y<ny;++y){
    for(int x=0;x<nx;++x){
#pragma acc loop vector
      for(int z=0;z<nz;++z){
        I_n[y*nx*nz+x*nz+z]=kap[y*nx*nz+x*nz+z]*rho[y*nx*nz+x*nz+z]*(J_band[y*nx*nz+x*nz+z]-Ss[y*nx*nz+x*nz+z]);
        St[y*nx*nz+x*nz+z] +=Ss[y*nx*nz+x*nz+z];
        Jt[y*nx*nz+x*nz+z] +=J_band[y*nx*nz+x*nz+z];
      }
    }
  }
  ACCH::UpdateCPU(I_n, nx*ny*nz*sizeof(double));
  ACCH::UpdateCPU(St, nx*ny*nz*sizeof(double));
  ACCH::UpdateCPU(Jt, nx*ny*nz*sizeof(double));

// 
  double inv_tau_0=1.0e1;
  ACCH::UpdateGPU(Fx, nx*ny*nz*sizeof(double));
  ACCH::UpdateGPU(Fy, nx*ny*nz*sizeof(double));
  ACCH::UpdateGPU(Fz, nx*ny*nz*sizeof(double));
  ACCH::UpdateGPU3D<double>(Qt, ny-yo, nx-xo, nz-zo);

#pragma acc parallel loop gang collapse(2) \
 present(this[:1], Fx[:nx*ny*nz], Fy[:nx*ny*nz], Fz[:nx*ny*nz], \
        I_n[:nx*ny*nz], Qt[:ny-yo][:nx-xo][:nz-zo], Qtemp[:ny][:nx][:nz][:2], \
        Tau[:nx*ny*nz])
  for(int y=0;y<ny-yo;++y){
    for(int x=0;x<nx-xo;++x){
#pragma acc loop vector
      for(int z=0;z<nz-zo;++z){
        double qf1=
          (
            (
              Fz[y*nx*nz + x*nz + z+zo]+
              Fz[y*nx*nz + (x+xo)*nz + (z+zo)]+
              Fz[(y+yo)*nx*nz + x*nz + (z+zo)]+
              Fz[(y+yo)*nx*nz + (x+xo)*nz + (z+zo)]
            )-(
              Fz[y*nx*nz + x*nz + z]+
              Fz[y*nx*nz + (x+xo)*nz + z]+
              Fz[(y+yo)*nx*nz + x*nz + z]+
              Fz[(y+yo)*nx*nz + (x+xo)*nz + z]
            )
          )*idz
          +
          (
            (
              Fx[y*nx*nz + (x+xo)*nz + z]+
              Fx[y*nx*nz + (x+xo)*nz + (z+zo)]+
              Fx[(y+yo)*nx*nz + (x+xo)*nz + z]+
              Fx[(y+yo)*nx*nz + (x+xo)*nz + (z+zo)]
            )-(
              Fx[y*nx*nz + x*nz + z]+
              Fx[y*nx*nz + x*nz + (z+zo)]+
              Fx[(y+yo)*nx*nz + x*nz + z]+
              Fx[(y+yo)*nx*nz + x*nz + (z+zo)]
            )
          )*idx
          +
          (
            (
              Fy[(y+yo)*nx*nz + x*nz + z]+
              Fy[(y+yo)*nx*nz + x*nz + (z+zo)]+
              Fy[(y+yo)*nx*nz + (x+xo)*nz + z]+
              Fy[(y+yo)*nx*nz + (x+xo)*nz + (z+zo)]
            )-(
              Fy[y*nx*nz + x*nz + z]+
              Fy[y*nx*nz + x*nz + (z+zo)]+
              Fy[y*nx*nz + (x+xo)*nz + z]+
              Fy[y*nx*nz + (x+xo)*nz + (z+zo)]
            )
          )*idy;
        qf1*=-0.25e0; 
        double qj1 = I_n[y*nx*nz+x*nz+z] +
                     I_n[y*nx*nz+x*nz+(z+zo)] +
                     I_n[y*nx*nz+(x+xo)*nz+z] +
                     I_n[y*nx*nz+(x+xo)*nz+(z+zo)] +
                     I_n[(y+yo)*nx*nz+x*nz+z] +
                     I_n[(y+yo)*nx*nz+x*nz+(z+zo)] +
                     I_n[(y+yo)*nx*nz+(x+xo)*nz+z] +
                     I_n[(y+yo)*nx*nz+(x+xo)*nz+(z+zo)];
        qj1*=0.5*PI; 
        double tau_local=tau(z+zo,x+xo,y+yo);
        double weight=exp(-tau_local*inv_tau_0);
        Qtemp[y][x][z+zo][0]=qf1;
        Qtemp[y][x][z+zo][1]=qj1;
        Qt[y][x][z-zo+1]+=weight*qj1+(1.0-weight)*qf1;
      }
    }
  }

/*
//#pragma acc parallel loop gang collapse(2) \
 present(this[:1], Fx[:nx*ny*nz], Fy[:nx*ny*nz], Fz[:nx*ny*nz], \
        I_n[:nx*ny*nz], Qt[:ny-yo][:nx-xo][:nz-zo], Qtemp[:ny][:nx][:nz][:2])
  for(int y=0;y<ny-yo;++y){
    for(int x=0;x<nx-xo;++x){
//#pragma acc loop vector
      for(int z=0;z<nz-zo;++z){
        double qf1=
          (
            (
              Fz[y*nx*nz + x*nz + z+zo]+
              Fz[y*nx*nz + (x+xo)*nz + (z+zo)]+
              Fz[(y+yo)*nx*nz + x*nz + (z+zo)]+
              Fz[(y+yo)*nx*nz + (x+xo)*nz + (z+zo)]
            )-(
              Fz[y*nx*nz + x*nz + z]+
              Fz[y*nx*nz + (x+xo)*nz + z]+
              Fz[(y+yo)*nx*nz + x*nz + z]+
              Fz[(y+yo)*nx*nz + (x+xo)*nz + z]
            )
          )*idz
          +
          (
            (
              Fx[y*nx*nz + (x+xo)*nz + z]+
              Fx[y*nx*nz + (x+xo)*nz + (z+zo)]+
              Fx[(y+yo)*nx*nz + (x+xo)*nz + z]+
              Fx[(y+yo)*nx*nz + (x+xo)*nz + (z+zo)]
            )-(
              Fx[y*nx*nz + x*nz + z]+
              Fx[y*nx*nz + x*nz + (z+zo)]+
              Fx[(y+yo)*nx*nz + x*nz + z]+
              Fx[(y+yo)*nx*nz + x*nz + (z+zo)]
            )
          )*idx
          +
          (
            (
              Fy[(y+yo)*nx*nz + x*nz + z]+
              Fy[(y+yo)*nx*nz + x*nz + (z+zo)]+
              Fy[(y+yo)*nx*nz + (x+xo)*nz + z]+
              Fy[(y+yo)*nx*nz + (x+xo)*nz + (z+zo)]
            )-(
              Fy[y*nx*nz + x*nz + z]+
              Fy[y*nx*nz + x*nz + (z+zo)]+
              Fy[y*nx*nz + (x+xo)*nz + z]+
              Fy[y*nx*nz + (x+xo)*nz + (z+zo)]
            )
          )*idy;
        qf1*=-0.25e0; 
        double qj1 = I_n[y*nx*nz+x*nz+z] +
                     I_n[y*nx*nz+x*nz+(z+zo)] +
                     I_n[y*nx*nz+(x+xo)*nz+z] +
                     I_n[y*nx*nz+(x+xo)*nz+(z+zo)] +
                     I_n[(y+yo)*nx*nz+x*nz+z] +
                     I_n[(y+yo)*nx*nz+x*nz+(z+zo)] +
                     I_n[(y+yo)*nx*nz+(x+xo)*nz+z] +
                     I_n[(y+yo)*nx*nz+(x+xo)*nz+(z+zo)];
        qj1*=0.5*PI; 

        double tau_local=tau(z+zo,x+xo,y+yo);
        double weight=exp(-tau_local*inv_tau_0);
        Qtemp[y][x][z+zo][0]=qf1;
        Qtemp[y][x][z+zo][1]=qj1;

        Qt[y][x][z-zo+1]+=weight*qj1+(1.0-weight)*qf1;
      }
    }
  }
*/
  ACCH::UpdateCPU3D<double>(Qt, ny-yo, nx-xo, nz-zo);
  ACCH::UpdateCPU4D<double>(Qtemp, ny, nx, nz, 2);

   if (save_col){
     ACCH::UpdateGPU3D<double>(Col_out, Nbands, col_nz, col_nvar);
     ACCH::UpdateGPU(J_band, nx*ny*nz*sizeof(double));
     ACCH::UpdateGPU(Ss, nx*ny*nz*sizeof(double));
     ACCH::UpdateGPU(kap, nx*ny*nz*sizeof(double));
     ACCH::UpdateGPU(abn, nx*ny*nz*sizeof(double));
     ACCH::UpdateGPU(sig, nx*ny*nz*sizeof(double));
     ACCH::UpdateGPU(B, nx*ny*nz*sizeof(double));
     ACCH::UpdateGPU(Tau, nx*ny*nz*sizeof(double));
     ACCH::UpdateGPU4D<double>(Qtemp, ny, nx, nz, 2);
     int col_bnd2 = col_bnd[2];
     int col_bnd3 = col_bnd[3];
     int col_bnd0 = col_bnd[0];
     int col_bnd1 = col_bnd[1];
#pragma acc parallel loop gang collapse(2) \
 present(this[:1], Col_out[:Nbands][:col_nz][:col_nvar], \
         J_band[:nx*ny*nz], Ss[:nx*ny*nz], kap[:nx*ny*nz], abn[:nx*ny*nz], \
         sig[:nx*ny*nz], B[:nx*ny*nz], Tau[:nx*ny*nz], Qtemp[:ny][:nx][:nz][:2])
     for (int y=col_bnd2;y<=col_bnd3;++y){
       for (int x=col_bnd0;x<=col_bnd1;++x){
#pragma acc loop vector
         for (int z=zl+zo;z<=zh;++z){
           int ind = (y-yl)*nx*nz + (x-xl)*nz + (z-zl);
           Col_out[band][z-2*zo][0] += J_band[ind]*avg_col;
           Col_out[band][z-2*zo][1] += Ss[ind]*avg_col;
           Col_out[band][z-2*zo][2] += kap[ind]*avg_col;
           Col_out[band][z-2*zo][3] += abn[ind]*avg_col;
           Col_out[band][z-2*zo][4] += sig[ind]*avg_col;
           Col_out[band][z-2*zo][5] += B[ind]*avg_col;
           Col_out[band][z-2*zo][6] += Tau[ind]*avg_col;
           Col_out[band][z-2*zo][7] += Qtemp[y][x][z][0]*avg_col;
           Col_out[band][z-2*zo][8] += Qtemp[y][x][z][1]*avg_col;
         }
       }
     }
     ACCH::UpdateCPU3D<double>(Col_out, Nbands, col_nz, col_nvar);
   }

  ACCH::Free4D<double>(Qtemp, ny, nx, nz, 2);
  ACCH::Free2D<double>(sbuf, ny, nx);
  ACCH::Free2D<double>(rbuf, ny, nx);
  //
  ttime=MPI_Wtime()-ttime;

  if((myrank==0)&&(verbose>2)) printf("tauscale time: %f %f %f %f \n",ttime,stime,atime,(stime+atime)/ttime);
}

void RTS::save_1D_avg(char * path, int iter,double time){

  // Save 1D avg
  
  static int ini_flag = 1;

  static MPI_Datatype x_subarray;

  if (ini_flag) {
    int array_of_sizes[1];
    int array_of_subsizes[1];
    int array_of_starts[1];

    array_of_sizes[0]=col_nzt;
    array_of_subsizes[0]=col_nz;
    array_of_starts[0]=col_offz;

    MPI_Type_create_subarray(1,array_of_sizes,array_of_subsizes,
                 array_of_starts,MPI_ORDER_FORTRAN,
                 MPI_FLOAT,&x_subarray);
    MPI_Type_commit(&x_subarray);

    ini_flag = 0;
  }

  double *** Col_out_glo = d3dim(0,Nbands-1,0,col_nz-1,0,col_nvar-1);
  memset(Col_out_glo[0][0],0,col_nz*Nbands*col_nvar*sizeof(double));

  MPI_Reduce(Col_out[0][0],Col_out_glo[0][0],col_nvar*Nbands*col_nz,MPI_DOUBLE,MPI_SUM,0,YZ_COMM);

  if(yz_rank==0){ // MPI_Reduce results are only meaningful on rank 0!

    int bufsize = Nbands*col_nz*col_nvar;
    float * iobuf = new float[bufsize];
 
    for(int band=0;band<Nbands;band++)
      for(int ind=0;ind<col_nz;ind++) 
        for(int var=0;var<col_nvar;var++)
          iobuf[band*(col_nz*col_nvar)+var*col_nz+ind] = (float) Col_out_glo[band][ind][var];

    char filename[128];
    sprintf(filename,"%s%s.%06d",path,"RT_mean1D",iter);

    MPI_File fhandle_mpi;

    MPI_File_open(XCOL_COMM,filename,MPI_MODE_CREATE | MPI_MODE_WRONLY,
          MPI_INFO_NULL,&fhandle_mpi);

    if( xcol_rank == 0 ){
      float header[8];            
      header[0] = (float) Nbands;
      header[1] = (float) col_nvar;
      header[2] = (float) col_nzt;      
      header[3] = (float) time;
      header[4] = (float) col_bnd[0];
      header[5] = (float) col_bnd[1];
      header[6] = (float) col_bnd[2];      
      header[7] = (float) col_bnd[3];
      MPI_File_write(fhandle_mpi,&header[0],8,MPI_FLOAT,MPI_STATUS_IGNORE);
    }
    
    int offset = 8*sizeof(float);
      
    MPI_File_set_view(fhandle_mpi,offset,MPI_FLOAT,x_subarray,(char*) "native", MPI_INFO_NULL);
    MPI_File_write_all(fhandle_mpi,&iobuf[0],bufsize,MPI_FLOAT, MPI_STATUS_IGNORE);
    MPI_File_close(&fhandle_mpi);  
    delete [] iobuf;
  }
  
  del_d3dim(Col_out_glo,0,Nbands-1,0,col_nz-1,0,col_nvar-1);

} 
