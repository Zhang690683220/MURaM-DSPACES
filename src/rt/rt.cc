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


RTS::~RTS(void)
{
  for(int j = 0;j < cart_sizes[2]; j++) delete [] comm_col[j];
  delete [] comm_col;

///*
  del_d2dim(I_o,yl,yh,xl,xh); // MHD Grid
  del_i3dim(tr_switch,yl,yh,xl,xh,zl,zh); // switch for the lowest point at which we have reached the transition region
  del_d3dim(lgTe,yl,yh,xl,xh,zl,zh); // temperature
  del_d3dim(lgPe,yl,yh,xl,xh,zl,zh); // temperature
  del_i3dim(T_ind,yl,yh,xl,xh,zl,zh); // temperature
  del_i3dim(P_ind,yl,yh,xl,xh,zl,zh); // temperature
  del_d3dim(rho,yl,yh,xl,xh,zl,zh);  // RT grid

  del_d3dim(Tau,yl,yh,xl,xh,zl,zh);  // RT grid
  del_d3dim(Qt,yl+yo,yh,xl+xo,xh,zl+zo,zh);   // MHD grid

  del_d3dim(Jt,yl,yh,xl,xh,zl,zh);   // MHD grid
  del_d3dim(St,yl,yh,xl,xh,zl,zh);   // MHD grid

// frequency dependent quantities
  del_d3dim(B,yl,yh,xl,xh,zl,zh);
  del_d3dim(kap,yl,yh,xl,xh,zl,zh);
  
  if (rttype ==1){
    del_d3dim(sig,yl,yh,xl,xh,zl,zh);
    del_d3dim(abn,yl,yh,xl,xh,zl,zh);
  }
  
  del_d3dim(J_band,yl,yh,xl,xh,zl,zh);
  del_d3dim(I_n,yl,yh,xl,xh,zl,zh);

  del_d3dim(Fx,yl,yh,xl,xh,zl,zh);
  del_d3dim(Fy,yl,yh,xl,xh,zl,zh);
  del_d3dim(Fz,yl,yh,xl,xh,zl,zh);

  del_d3dim(Col_out,0,Nbands-1,0,col_nz-1,0,col_nvar-1);
   
  del_i5dim(numits,0,Nbands-1,FWD,BWD,RIGHT,LEFT,UP,DOWN,0,NMU-1);
  
  if (NDIM>1){
    del_r6dim(x_sbuf,0,Nbands-1,FWD,BWD,RIGHT,LEFT,UP,DOWN,0,NMU-1,0,nz*ny-1);
    del_r6dim(x_rbuf,0,Nbands-1,FWD,BWD,RIGHT,LEFT,UP,DOWN,0,NMU-1,0,nz*ny-1);
    del_r6dim(x_oldbuf,0,Nbands-1,FWD,BWD,RIGHT,LEFT,UP,DOWN,0,NMU-1,0,nz*ny-1);
  }

  if (NDIM==3){
    del_r6dim(y_sbuf,0,Nbands-1,FWD,BWD,RIGHT,LEFT,UP,DOWN,0,NMU-1,0,nz*nx-1);
    del_r6dim(y_rbuf,0,Nbands-1,FWD,BWD,RIGHT,LEFT,UP,DOWN,0,NMU-1,0,nz*nx-1);
    del_r6dim(y_oldbuf,0,Nbands-1,FWD,BWD,RIGHT,LEFT,UP,DOWN,0,NMU-1,0,nz*nx-1);
  }

  del_r6dim(z_sbuf,0,Nbands-1,FWD,BWD,RIGHT,LEFT,UP,DOWN,0,NMU-1,0,nx*ny-1);
  del_r6dim(z_rbuf,0,Nbands-1,FWD,BWD,RIGHT,LEFT,UP,DOWN,0,NMU-1,0,nx*ny-1);
  del_r6dim(z_oldbuf,0,Nbands-1,FWD,BWD,RIGHT,LEFT,UP,DOWN,0,NMU-1,0,nx*ny-1);
//*/
  del_d1dim(tab_T,0,NT-1);
  del_d1dim(tab_p,0,Np-1);
  del_d1dim(invT_tab,0,NT-1);
  del_d1dim(invP_tab,0,Np-1);
  
  if(fullodf)
    del_f1dim(nu_tab,0,Nlam-1);
    
  if (N5000){
    del_f1dim(B_5000_tab,0,NT-1); 
    del_f2dim(kap_5000_tab,0,NT-1,0,Np-1); 
  }

  del_f3dim(kap_tab,0,Nbands-1,0,NT-1,0,Np-1);
  del_f2dim(B_tab,0,Nbands-1,0,NT-1);

}
 
double RTS::tau(int z,int x,int y){
  double Tau_local = Tau[y][x][z]+Tau[y-yo][x][z]+
  Tau[y][x-xo][z]+Tau[y-yo][x-xo][z]+
  Tau[y][x][z-zo]+Tau[y-yo][x][z-zo]+
  Tau[y][x-xo][z-zo]+Tau[y-yo][x-xo][z-zo];
  Tau_local *= 0.125;

  return Tau_local;
}

double RTS::Qtot(int z,int x,int y)
{  
  return Qt[y][x][z];
}


double RTS::Jtot(int z,int x,int y){
  double J_local = Jt[y][x][z]+Jt[y-yo][x][z]+
    Jt[y][x-xo][z]+Jt[y-yo][x-xo][z]+
    Jt[y][x][z-zo]+Jt[y-yo][x][z-zo]+
    Jt[y][x-xo][z-zo]+Jt[y-yo][x-xo][z-zo];
  J_local *= 0.125;

  return J_local;
}

double RTS::Stot(int z,int x,int y){
  double S_local = St[y][x][z]+St[y-yo][x][z]+
    St[y][x-xo][z]+St[y-yo][x-xo][z]+
    St[y][x][z-zo]+St[y-yo][x][z-zo]+
    St[y][x-xo][z-zo]+St[y-yo][x-xo][z-zo];
  S_local *= 0.125;

  return S_local;
}


double RTS::Iout(int x,int y)
{
  double io = I_o[y][x] + I_o[y-yo][x] + I_o[y][x-xo] + I_o[y-yo][x-xo];
  io *= 0.25;
  return io;
}

RTS::RTS(GridData&Grid,RunData &Run,PhysicsData &Physics){
  
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

  load_bins(Run.kap_name);

// output intensity
  I_o=d2dim(yl,yh,xl,xh);        // RT grid

  Fr_mean=d1dim(0,Nbands-1);
  gFr_mean=d1dim(0,Nbands-1);
  memset(Fr_mean,0,Nbands*sizeof(double));
  memset(gFr_mean,0,Nbands*sizeof(double));

  lgTe=d3dim(yl,yh,xl,xh,zl,zh); // temperature
  lgPe=d3dim(yl,yh,xl,xh,zl,zh); // pressure
  T_ind=i3dim(yl,yh,xl,xh,zl,zh); // pressure
  P_ind=i3dim(yl,yh,xl,xh,zl,zh); // pressure
  rho=d3dim(yl,yh,xl,xh,zl,zh);  // RT grid
  tr_switch =i3dim(yl,yh,xl,xh,zl,zh);

  Tau=d3dim(yl,yh,xl,xh,zl,zh);  // RT grid

  Qt=d3dim(yl+yo,yh,xl+xo,xh,zl+zo,zh);   // MHD grid
  Jt=d3dim(yl,yh,xl,xh,zl,zh);   // MHD grid
  St=d3dim(yl,yh,xl,xh,zl,zh);   // MHD grid

// frequency dependent quantities
  B=d3dim(yl,yh,xl,xh,zl,zh);
  kap=d3dim(yl,yh,xl,xh,zl,zh);
  I_n=d3dim(yl,yh,xl,xh,zl,zh);

  if (rttype==0){
    sig=kap;
    abn=kap;
  } else {
    sig=d3dim(yl,yh,xl,xh,zl,zh);
    abn=d3dim(yl,yh,xl,xh,zl,zh);
  }

  J_band=d3dim(yl,yh,xl,xh,zl,zh);

  Fx=d3dim(yl,yh,xl,xh,zl,zh);
  Fy=d3dim(yl,yh,xl,xh,zl,zh);
  Fz=d3dim(yl,yh,xl,xh,zl,zh);
 
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

  Col_out = d3dim(0,Nbands-1,0,col_nz-1,0,col_nvar-1);

  memset(Col_out[0][0],0,col_nz*Nbands*col_nvar*sizeof(double));

  numits = i5dim(0,Nbands-1,FWD,BWD,RIGHT,LEFT,UP,DOWN,0,NMU-1);

  if (NDIM==3){
    y_sbuf=r6dim(0,Nbands-1,FWD,BWD,RIGHT,LEFT,UP,DOWN,0,NMU-1,0,nz*nx-1);
    y_rbuf=r6dim(0,Nbands-1,FWD,BWD,RIGHT,LEFT,UP,DOWN,0,NMU-1,0,nz*nx-1);
    y_oldbuf=r6dim(0,Nbands-1,FWD,BWD,RIGHT,LEFT,UP,DOWN,0,NMU-1,0,nz*nx-1);
  }

  if (NDIM>1){
    x_sbuf=r6dim(0,Nbands-1,FWD,BWD,RIGHT,LEFT,UP,DOWN,0,NMU-1,0,nz*ny-1);
    x_rbuf=r6dim(0,Nbands-1,FWD,BWD,RIGHT,LEFT,UP,DOWN,0,NMU-1,0,nz*ny-1);
    x_oldbuf=r6dim(0,Nbands-1,FWD,BWD,RIGHT,LEFT,UP,DOWN,0,NMU-1,0,nz*ny-1);
  }
  z_sbuf=r6dim(0,Nbands-1,FWD,BWD,RIGHT,LEFT,UP,DOWN,0,NMU-1,0,nx*ny-1);
  z_rbuf=r6dim(0,Nbands-1,FWD,BWD,RIGHT,LEFT,UP,DOWN,0,NMU-1,0,nx*ny-1);
  z_oldbuf=r6dim(0,Nbands-1,FWD,BWD,RIGHT,LEFT,UP,DOWN,0,NMU-1,0,nx*ny-1);
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

    if(((aM[2]<=aM[0])&&((aM[2]<=aM[1])||(NDIM==2)))||(NDIM==1)){
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
      for(int i=0; i< cart_sizes[0]; i++) ranks[i]=colranks[j][k][cart_sizes[0]-1-i];
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

    tab_T = d1dim(0,NT-1);
    tab_p = d1dim(0,Np-1);
     
    invT_tab = d1dim(0,NT-1);
    invP_tab = d1dim(0,Np-1);
    
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
      kap_5000_tab = f2dim(0,NT-1,0,Np-1); 
      B_5000_tab = f1dim(0,NT-1);
       
      fp_rt.read((char*)&kap_5000_tab[0][0],NT*Np*sizeof(float));
      fp_rt.read((char*)&B_5000_tab[0],NT*sizeof(float));
    }
      
    kap_tab = f3dim(0,Nbands-1,0,NT-1,0,Np-1);
    B_tab   = f2dim(0,Nbands-1,0,NT-1);
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
  
  double ** I_band;
  I_band=d2dim(yl,yh,xl,xh);        // RT grid

  memset(I_band[yl]+xl,0,nx*ny*sizeof(double));
  memset(I_o[yl]+xl,0,nx*ny*sizeof(double));

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
      for (int y=yl;y<=yh;y++)
        for (int x=xl;x<=xh;x++)
          I_o[y][x] +=I_band[y][x];
      }
    }
    // If cont_bin = 0 and need_I = 1 we want the output intensity to be the 0th continuum bin
    if (cont_bin==0){
      get_Tau_and_Iout(Grid, Run, Physics,DZ,B_tab[0],kap_tab[0],I_band,need_I);
      for (int y=yl;y<=yh;y++)
        for (int x=xl;x<=xh;x++)
          I_o[y][x] +=I_band[y][x];
    }
    // if we have the band, we want a tau5000 reference wavelength. If cont_bin=1 and need_I we want
    // output intensity to be 5000A.
    
    if (N5000){
      int I5000_out = 0;
      if ((cont_bin==1)&&(need_I==1))
        I5000_out = 1;
      
      get_Tau_and_Iout(Grid, Run, Physics,DZ,B_5000_tab,kap_5000_tab,I_band,I5000_out);
     
      if ((cont_bin==1)&&(need_I==1))
        for (int y=yl;y<=yh;y++)
          for (int x=xl;x<=xh;x++)
            I_o[y][x] += I_band[y][x];
    }

    calc_Qtot_and_Tau(Grid, Run, Physics);
    del_d2dim(I_band,yl,yh,xl,xh);
    return dt_rad;
  }

// *****************************************************************
// *        interpolate opacity and Planck function (B)            *
// *****************************************************************

  cState *U=Grid.U;

  double N = pow(2,NDIM);

  for(int y=yl;y<=yh;++y){
    for(int x=xl;x<=xh;++x){
      int off0 = x*next[1]+y*next[2];
      int off1 = x*next[1]+(y+yo)*next[2];
      int off2 = (x+xo)*next[1]+y*next[2];
      int off3 = (x+xo)*next[1]+(y+yo)*next[2];
      for(int z=zl;z<=zh;++z){
        int inode[]={off0+z,off0+z+zo,off2+z,off2+z+zo,off1+z,off1+z+zo,off3+z,off3+z+zo};
        double Tm=0.0,pm=0.0,rm=0.0;
        for(int l=0;l<N;++l){
          Tm+=Grid.temp[inode[l]];
          pm+=Grid.pres[inode[l]];
          rm+=U[inode[l]].d;
        }

        Tm /= N;
        pm /= N;
        rm /= N;

        //disbale RT if Temp > 2e4 above the photosphere
        double pswitch = min(max(pm-Pres_TR,0.0),1.0);
        double tswitch = min(max(Temp_TR-Tm,0.0),1.0);
        tr_switch[y][x][z] = (int) max(pswitch,tswitch);
        
        lgTe[y][x][z]=log(Tm);
        lgPe[y][x][z]=log(pm);
        rho[y][x][z] =rm;
    
        lgTe[y][x][z] = min(max(lgTe[y][x][z],(double) tab_T[0]),(double) tab_T[NT-1]);
        lgPe[y][x][z] = min(max(lgPe[y][x][z],(double) tab_p[0]),(double) tab_p[Np-1]);
      }
    
      for(int z=zl;z<=zh;++z){

    // Search table once for all RT bands
        int l=0;
        int m=0;
        if(lgTe[y][x][z]<tab_T[0])
          l=0;
        else if(lgTe[y][x][z]>tab_T[NT-1])
          l=NT-2;
        else
          for (l=0; l<=NT-2; l++)
            if ((lgTe[y][x][z] >= tab_T[l]) && (lgTe[y][x][z] <= tab_T[l+1]))
              break;

        if(lgPe[y][x][z]<tab_p[0])
          m=0;
        else if(lgPe[y][x][z]>tab_p[Np-1])
          m=Np-2;
        else
          for (m=0; m<=Np-2; m++)
             if ((lgPe[y][x][z] >= tab_p[m]) && (lgPe[y][x][z] <= tab_p[m+1]))
               break;

        T_ind[y][x][z] = l;
        P_ind[y][x][z] = m;

      }
    }
  }
  
  // Begin RT loop, work band by band,
  memset(St[yl][xl]+zl,0,nx*ny*nz*sizeof(double));
  memset(Jt[yl][xl]+zl,0,nx*ny*nz*sizeof(double));
  memset(Qt[yl+yo][xl+xo]+zl+zo,0,(nx-xo)*(ny-yo)*(nz-zo)*sizeof(double)); // MHD grid
 
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

    for(int y=yl;y<=yh;++y){
      for(int x=xl;x<=xh;++x){
        for(int z=zl;z<=zh;++z){
      
          int l = T_ind[y][x][z];
          int m = P_ind[y][x][z];
      
          double xt = (lgTe[y][x][z]-tab_T[l])*invT_tab[l];
          double xp = (lgPe[y][x][z]-tab_p[m])*invP_tab[m];
      
          // Interpolate for kappa and B
          B[y][x][z]=exp(xt*B_tab[band][l+1]+(1.-xt)*B_tab[band][l]);
      
          kap[y][x][z] = 
          exp(xt*(xp*kap_tab[band][l+1][m+1]+(1.-xp)*kap_tab[band][l+1][m])+
          (1.-xt)*(xp*kap_tab[band][l][m+1]+(1.-xp)*kap_tab[band][l][m]));
    
          // TR switch turn off kappa and B
          kap[y][x][z] *= tr_switch[y][x][z];
          B[y][x][z]   *= tr_switch[y][x][z];
        }
      }
    }

// *****************************************************************
// *    update diffusion-approx. boundary condition at bottom      *
// *****************************************************************
  for(int YDIR=FWD;YDIR<=BWD;++YDIR)
    for(int XDIR=RIGHT;XDIR<=LEFT;++XDIR){
      if(isgbeg[0]==1)
        for(int l=0;l<NMU;++l)
        for(int y=0;y<ny;++y)
          for(int x=0;x<nx;++x)
            z_rbuf[band][YDIR][XDIR][UP][l][x*ny+y]=B[yl+y][xl+x][zl];
          
// *****************************************************************
// *    no incoming radiation on the top                           *
// *****************************************************************
          if(isgend[0]==1) memset(z_rbuf[band][YDIR][XDIR][DOWN][0],0,NMU*nx*ny*sizeof(real));
          }
// *****************************************************************
// *  solve the transfer                                           *
// *****************************************************************
  driver(DZ,DX,DY,band); 

// *****************************************************************
// *  outgoing radiative flux                                      *
// *****************************************************************

  gFr_mean[band] = 0.0;

  if(isgend[0]==1){
    for(int y=yl;y<=yh-yo;++y)
      for(int x=xl;x<=xh-xo;++x)
        gFr_mean[band]+=Fz[y][x][zh]+Fz[y+yo][x][zh]+
          Fz[y][x+xo][zh]+Fz[y+yo][x+xo][zh];
    gFr_mean[band]*=0.25;
  }
 
  // If I am saving 5000A intensity then wipe after last bin. if I am saving the continuum bin then wipe the second last bin.

  if (((cont_bin==1)&&(band==0))||((cont_bin==0)&&(band==1)))
    for(int y=yl;y<=yh;++y)
      for(int x=xl;x<=xh;++x)
        I_o[y][x] = 0.0;

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
   
    if ((cont_bin==1)&&(need_I==1))
      for (int y=yl;y<=yh;y++)
        for (int x=xl;x<=xh;x++)
          I_o[y][x] +=I_band[y][x];
  }

  calc_Qtot_and_Tau(Grid, Run, Physics);
  
  del_d2dim(I_band,yl,yh,xl,xh);
  
  if (Run.NeedsSlice() && Run.RT_HAVG)
    save_1D_avg(Run.path_2D,Run.globiter,Run.time);  

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

  double qsum=0.0;

  for(int y=yl+yo;y<=yh;++y) // loop over MHD grid
    for(int x=xl+xo;x<=xh;++x){
      int off0 = x*next[1]+y*next[2];
      for(int z=zl+zo;z<=zh;++z){
        int node = off0+z;

        Grid.Tau[node]=tau(z,x,y);
        Grid.Jtot[node]=Jtot(z,x,y);
        Grid.Stot[node]=Stot(z,x,y);
        
        double scale = pow(Grid.Tau[node],2)*tr_switch[y][x][z];
        scale = scale/(scale + tau_min);

        double Qt_step=Qt[y][x][z]*scale;
        double inv_dt=fabs(Qt_step)/U[node].e;
        
        dt_rad=(dt_rad>inv_dt)?dt_rad:inv_dt;

        Grid.Qtot[node]=Qt_step;

        qsum += Qt_step;
    }
  }

  exchange_single(Grid,Grid.Tau);

  double Fqrad;

  MPI_Allreduce(&qsum,&Fqrad,1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  Fqrad=-Fqrad*Grid.dx[0]/(Grid.gsize[1]*Grid.gsize[2]);

  if(myrank==0) fprintf(stdout,"RT energy flux: %21.15E %21.15E\n",F_o,Fqrad);
  
  F_o=Fqrad;
  dt_rad=Physics.rt[i_rt_cfl]/dt_rad;
  
}

void RTS::driver(double DZ, double DX, double DY, int band){
  double etime=0.0,atime=0.0,cmp_time1=0.0,cmp_time2=0.0,buf_time=0.0,err_time=0.0,flx_time=0.0,tau_time=0.0; 
  double ttime=MPI_Wtime();
  
  const int stride[2]={nx*nz,nz}; 
  int stepvec[3][4][3] = { {{1,0,0},{1,0,1},{1,1,0},{1,1,1}},
               {{0,1,0},{1,1,0},{0,1,1},{1,1,1}},
               {{0,0,1},{0,1,1},{1,0,1},{1,1,1}} };
  double ** coeff=d2dim(0,1,0,nx*ny*nz-1);
  
  memset(I_n[yl][xl]+zl,0,nx*ny*nz*sizeof(double)); // is there a better way?
  memset(Fz[yl][xl]+zl,0,nx*ny*nz*sizeof(double));
  memset(Fx[yl][xl]+zl,0,nx*ny*nz*sizeof(double));
  memset(Fy[yl][xl]+zl,0,nx*ny*nz*sizeof(double));
  memset(J_band[yl][xl]+zl,0,nx*ny*nz*sizeof(double));

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
              int off[4];
          int i_nu=0;
          double *ii=I_n[yl][xl];
          if(NDIM==3){
            for(int i=0;i<4;i++) off[i]=iystep[i]*stride[0]+ixstep[i]*stride[1]+izstep[i];
            for(int yi=yi_i;yi!=yi_f+ystep;yi=yi+ystep)
              for(int xi=xi_i;xi!=xi_f+xstep;xi=xi+xstep){
                int xyoff=(yi-yl)*stride[0]+(xi-xl)*stride[1];
                for(int zi=zi_i;zi!=zi_f+zstep;zi=zi+zstep){
                  int ind=xyoff+zi;
                  double I_upw=c[0]*ii[ind-off[0]]+c[1]*ii[ind-off[1]]+c[2]*ii[ind-off[2]]+c[3]*ii[ind-off[3]];
                  ii[ind]=I_upw*coeff[0][i_nu]+coeff[1][i_nu];
                  i_nu+=1;
                }
              }
          }

          if(NDIM==2){
            for(int i=0;i<4;i++) off[i]=ixstep[i]*stride[1]+izstep[i];
            for(int xi=xi_i;xi!=xi_f+xstep;xi=xi+xstep){
              int xoff=(xi-xl)*stride[1];
              for(int zi=zi_i;zi!=zi_f+zstep;zi=zi+zstep){
                int ind=xoff+zi;
                double I_upw=c[0]*ii[ind-off[0]]+c[1]*ii[ind-off[1]]+c[2]*ii[ind-off[2]]+c[3]*ii[ind-off[3]];
                ii[ind]=I_upw*coeff[0][i_nu]+coeff[1][i_nu];
                i_nu+=1;
              }
            }
          }

          if(NDIM==1){
            for(int i=0;i<4;i++) off[i]=izstep[i];
            for(int zi=zi_i;zi!=zi_f+zstep;zi=zi+zstep){
              double I_upw=c[0]*ii[zi-off[0]]+c[1]*ii[zi-off[1]]+c[2]*ii[zi-off[2]]+c[3]*ii[zi-off[3]];
              ii[zi]=I_upw*coeff[0][i_nu]+coeff[1][i_nu];
              i_nu+=1;
            }
          }

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
  }

  ttime=MPI_Wtime()-ttime;  
  if((myrank==0) && (verbose>2))
    fprintf(stdout,"rt_driver time : %f %f %f %f %f %f %f %f %f %f \n",ttime,cmp_time1,
       cmp_time2,buf_time,err_time,flx_time,tau_time,etime,atime,(etime+atime)/ttime);
   

  call_count+=1;
  del_d2dim(coeff,0,1,0,nx*ny*nz-1);
}

void RTS::interpol(int zi_i,int zi_f,int zstep,int xi_i,int xi_f,int xstep,
           int yi_i,int yi_f,int ystep,int l,double** coeff, double *** Ss)
{

  double ds3=ds_upw[l]*inv3,ds6=ds_upw[l]*inv6;
  double c[]={a_00[ibase[l]][l],a_01[ibase[l]][l],a_10[ibase[l]][l],a_11[ibase[l]][l]};

  int zmin=(zi_i<zi_f)?zi_i:zi_f;
  int zmax=(zi_i>zi_f)?zi_i:zi_f;

  double source[zmax+1],expo[zmax+1],dt[zmax+1],r_upw[zmax+1],k_upw[zmax+1],S_upw[zmax+1],r0[zmax+1],k0[zmax+1],S0[zmax+1];

  int i_nu=0;

  if(NDIM==3){
    for(int yi=yi_i;yi!=yi_f+ystep;yi=yi+ystep){
      for(int xi=xi_i;xi!=xi_f+xstep;xi=xi+xstep){
        for(int zi=zmin;zi<=zmax;++zi){
          r_upw[zi]=c[0]*rho[yi-iystep[0]][xi-ixstep[0]][zi-izstep[0]]+
          c[1]*rho[yi-iystep[1]][xi-ixstep[1]][zi-izstep[1]]+  
          c[2]*rho[yi-iystep[2]][xi-ixstep[2]][zi-izstep[2]]+ 
          c[3]*rho[yi-iystep[3]][xi-ixstep[3]][zi-izstep[3]];
          
          k_upw[zi]=c[0]*kap[yi-iystep[0]][xi-ixstep[0]][zi-izstep[0]]+
          c[1]*kap[yi-iystep[1]][xi-ixstep[1]][zi-izstep[1]]+
          c[2]*kap[yi-iystep[2]][xi-ixstep[2]][zi-izstep[2]]+
          c[3]*kap[yi-iystep[3]][xi-ixstep[3]][zi-izstep[3]];

          S_upw[zi]=c[0]*Ss[yi-iystep[0]][xi-ixstep[0]][zi-izstep[0]]+
          c[1]*Ss[yi-iystep[1]][xi-ixstep[1]][zi-izstep[1]]+
          c[2]*Ss[yi-iystep[2]][xi-ixstep[2]][zi-izstep[2]]+
          c[3]*Ss[yi-iystep[3]][xi-ixstep[3]][zi-izstep[3]];

          r0[zi]=rho[yi][xi][zi];
          k0[zi]=kap[yi][xi][zi];
          S0[zi]=Ss[yi][xi][zi];
        }
        for(int zi=zmin;zi<=zmax;++zi){
          dt[zi]=ds3*(k_upw[zi]*r_upw[zi]+k0[zi]*r0[zi])+ds6*(k0[zi]*r_upw[zi]+k_upw[zi]*r0[zi]);
          expo[zi]=exp(-dt[zi]);
          double w0=(1.0-expo[zi])/dt[zi];
          source[zi]=S0[zi]*(1.0-w0)+S_upw[zi]*(w0-expo[zi]);
        }   

        for(int zi=zi_i;zi!=zi_f+zstep;zi=zi+zstep){
          coeff[0][i_nu]=expo[zi];
          if (dt[zi] > dtau_min)
            coeff[1][i_nu]= source[zi];
          else
            coeff[1][i_nu] = 0.5*dt[zi]*(S0[zi]+S_upw[zi]);
            
          i_nu+=1;
        }
      }
    }
  }

  if(NDIM==1){
    for(int zi=zmin;zi<=zmax;++zi){
      r_upw[zi]=c[0]*rho[0][0][zi-izstep[0]]+
        c[1]*rho[0][0][zi-izstep[1]]+
        c[2]*rho[0][0][zi-izstep[2]]+
        c[3]*rho[0][0][zi-izstep[3]];

      k_upw[zi]=c[0]*kap[0][0][zi-izstep[0]]+
        c[1]*kap[0][0][zi-izstep[1]]+
        c[2]*kap[0][0][zi-izstep[2]]+
        c[3]*kap[0][0][zi-izstep[3]];

      S_upw[zi]=c[0]*Ss[0][0][zi-izstep[0]]+
        c[1]*Ss[0][0][zi-izstep[1]]+
        c[2]*Ss[0][0][zi-izstep[2]]+
        c[3]*Ss[0][0][zi-izstep[3]];

      r0[zi]=rho[0][0][zi];
      k0[zi]=kap[0][0][zi];
      S0[zi]=Ss[0][0][zi];
    }

    for(int zi=zmin;zi<=zmax;++zi){
      dt[zi]=ds3*(k_upw[zi]*r_upw[zi]+k0[zi]*r0[zi])+ds6*(k0[zi]*r_upw[zi]+k_upw[zi]*r0[zi]);
      expo[zi]=exp(-dt[zi]);
      double w0=(1.0-expo[zi])/dt[zi];
      source[zi]=S0[zi]*(1.0-w0)+S_upw[zi]*(w0-expo[zi]);
    }

    for(int zi=zi_i;zi!=zi_f+zstep;zi=zi+zstep){
      coeff[0][i_nu]=expo[zi];
      if (dt[zi] > dtau_min)
        coeff[1][i_nu]=source[zi];
      else
        coeff[1][i_nu] = 0.5*dt[zi]*(S0[zi]+S_upw[zi]);

      i_nu+=1;
    }
  }

  if(NDIM==2){
    for(int xi=xi_i;xi!=xi_f+xstep;xi=xi+xstep){
      for(int zi=zmin;zi<=zmax;++zi){
        r_upw[zi]=c[0]*rho[0][xi-ixstep[0]][zi-izstep[0]]+
          c[1]*rho[0][xi-ixstep[1]][zi-izstep[1]]+
          c[2]*rho[0][xi-ixstep[2]][zi-izstep[2]]+
          c[3]*rho[0][xi-ixstep[3]][zi-izstep[3]];

        k_upw[zi]=c[0]*kap[0][xi-ixstep[0]][zi-izstep[0]]+
          c[1]*kap[0][xi-ixstep[1]][zi-izstep[1]]+
          c[2]*kap[0][xi-ixstep[2]][zi-izstep[2]]+
          c[3]*kap[0][xi-ixstep[3]][zi-izstep[3]];

        S_upw[zi]=c[0]*Ss[0][xi-ixstep[0]][zi-izstep[0]]+
          c[1]*Ss[0][xi-ixstep[0]][zi-izstep[1]]+
          c[2]*Ss[0][xi-ixstep[0]][zi-izstep[2]]+
          c[3]*Ss[0][xi-ixstep[0]][zi-izstep[3]];

        r0[zi]=rho[0][xi][zi];
        k0[zi]=kap[0][xi][zi];
        S0[zi]=Ss[0][xi][zi];
      }

      for(int zi=zmin;zi<=zmax;++zi){
        dt[zi]=ds3*(k_upw[zi]*r_upw[zi]+k0[zi]*r0[zi])+ds6*(k0[zi]*r_upw[zi]+k_upw[zi]*r0[zi]);
        expo[zi]=exp(-dt[zi]);
        double w0=(1.0-expo[zi])/dt[zi];
        source[zi]=S0[zi]*(1.0-w0)+S_upw[zi]*(w0-expo[zi]);
     }

      for(int zi=zi_i;zi!=zi_f+zstep;zi=zi+zstep){
        coeff[0][i_nu]=expo[zi];
        coeff[1][i_nu]=source[zi];
      
        if (dt[zi] > dtau_min)
          coeff[1][i_nu]=source[zi];
        else
          coeff[1][i_nu] = 0.5*dt[zi]*(S0[zi]+S_upw[zi]);

        i_nu+=1;
      }
    }
  }
}

double RTS::error(int band,int l,int ZDIR,int XDIR,int YDIR,double I_min)
{
  int z0=(ZDIR==UP),z1=(ZDIR==DOWN); 
  int y0=(YDIR==FWD),y1=(YDIR==BWD); 
  double err_max=0.0;
//
  if (NDIM==3){
    real *ysb=y_sbuf[band][YDIR][XDIR][ZDIR][l],*yob=y_oldbuf[band][YDIR][XDIR][ZDIR][l];
    for(int z=0+z0;z<nz-z1;++z)
      for(int x=0;x<nx;++x){
        int ind=z*nx+x;
        err_max=max(err_max,fabs(((double) (ysb[ind]-yob[ind]))/max(I_min,(double) yob[ind])));
     }
  }
//
  if (NDIM>1){
  real *xsb=x_sbuf[band][YDIR][XDIR][ZDIR][l],*xob=x_oldbuf[band][YDIR][XDIR][ZDIR][l];
  for(int z=0+z0;z<nz-z1;++z)
    for(int y=0+y0;y<ny-y1;++y){
      int ind=z*ny+y;
      err_max=max(err_max,fabs(((double) (xsb[ind]-xob[ind]))/max( I_min,(double)xob[ind])));
    }
  }
//
  real *zsb=z_sbuf[band][YDIR][XDIR][ZDIR][l],*zob=z_oldbuf[band][YDIR][XDIR][ZDIR][l];
  for(int ind=0;ind<nx*ny;++ind) err_max=max(err_max,fabs(((double) (zsb[ind]-zob[ind]))/max(I_min,(double)zob[ind])));
//
  return err_max;
}

void RTS::readbuf(int band,int l,int ZDIR,int XDIR,int YDIR)
{//RHC
  if(NDIM==3){
    int y0=(YDIR==FWD)?0:ny-1;
    for(int x=0;x<nx;++x) for(int z=0;z<nz;++z) I_n[yl+y0][xl+x ][zl+z ]= (double) y_rbuf[band][YDIR][XDIR][ZDIR][l][z*nx+x];
  memcpy(y_oldbuf[band][YDIR][XDIR][ZDIR][l],y_sbuf[band][YDIR][XDIR][ZDIR][l],nz*nx*sizeof(real));
  }

  if(NDIM>1){
    int x0=(XDIR==RIGHT)?0:nx-1;
    for(int y=0;y<ny;++y) for(int z=0;z<nz;++z) I_n[yl+y ][xl+x0][zl+z ]= (double) x_rbuf[band][YDIR][XDIR][ZDIR][l][z*ny+y];
    memcpy(x_oldbuf[band][YDIR][XDIR][ZDIR][l],x_sbuf[band][YDIR][XDIR][ZDIR][l],ny*nz*sizeof(real));  
  }
  
  int z0=(ZDIR==UP)?0:nz-1;
  for(int y=0;y<ny;++y) for(int x=0;x<nx;++x) I_n[yl+y ][xl+x ][zl+z0]= (double) z_rbuf[band][YDIR][XDIR][ZDIR][l][x*ny+y];
  memcpy(z_oldbuf[band][YDIR][XDIR][ZDIR][l],z_sbuf[band][YDIR][XDIR][ZDIR][l],nx*ny*sizeof(real));  
}

void RTS::writebuf(int band, int l,int ZDIR,int XDIR,int YDIR)
{
  if (NDIM==3){
    int y0=(YDIR==FWD)?ny-1:0;
    for(int x=0;x<nx;++x) for(int z=0;z<nz;++z) y_sbuf[band][YDIR][XDIR][ZDIR][l][z*nx+x]=(real) I_n[yl+y0][xl+x ][zl+z ];
  }

  if (NDIM>1){
    int x0=(XDIR==RIGHT)?nx-1:0;
    for(int y=0;y<ny;++y) for(int z=0;z<nz;++z) x_sbuf[band][YDIR][XDIR][ZDIR][l][z*ny+y]=(real) I_n[yl+y ][xl+x0][zl+z ];
  }

  int z0=(ZDIR==UP)?nz-1:0;
  for(int y=0;y<ny;++y) for(int x=0;x<nx;++x) z_sbuf[band][YDIR][XDIR][ZDIR][l][x*ny+y]=(real) I_n[yl+y ][xl+x ][zl+z0];
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

  if (NDIM==3){
	// y-direction
	dest_rk=(YDIR==FWD)?rightr[2]:leftr[2];
	source_rk=(YDIR==FWD)?leftr[2]:rightr[2];
	MPI_Irecv(y_rbuf[band][YDIR][XDIR][ZDIR][l],nx*nz,REALTYPE,source_rk,tag2,MPI_COMM_WORLD,r2+0);
	MPI_Isend(y_sbuf[band][YDIR][XDIR][ZDIR][l],nx*nz,REALTYPE,dest_rk,tag2,MPI_COMM_WORLD,r2+1);
	MPI_Waitall(2,r2,s2);
  
    z0=(ZDIR==UP)?nz-1:0;
    y0=(YDIR==FWD)?0:ny-1;
    for(int x=0;x<nx;++x)
      z_sbuf[band][YDIR][XDIR][ZDIR][l][x*ny+y0]=y_rbuf[band][YDIR][XDIR][ZDIR][l][z0*nx+x];
  
    x0=(XDIR==RIGHT)?nx-1:0;
    y0=(YDIR==FWD)?0:ny-1;
    for(int z=0;z<nz;++z)
        x_sbuf[band][YDIR][XDIR][ZDIR][l][z*ny+y0]=y_rbuf[band][YDIR][XDIR][ZDIR][l][z*nx+x0];
  }

  if (NDIM >1){
    // x-direction
    dest_rk=(XDIR==RIGHT)?rightr[1]:leftr[1];
    source_rk=(XDIR==RIGHT)?leftr[1]:rightr[1];
    MPI_Irecv(x_rbuf[band][YDIR][XDIR][ZDIR][l],nz*ny,REALTYPE,source_rk,tag1,MPI_COMM_WORLD,r1+0);
    MPI_Isend(x_sbuf[band][YDIR][XDIR][ZDIR][l],nz*ny,REALTYPE,dest_rk,tag1,MPI_COMM_WORLD,r1+1);

    MPI_Waitall(2,r1,s1);
 
    x0=(XDIR==RIGHT)?0:nx-1;
    z0=(ZDIR==UP)?nz-1:0;
    for(int y=0;y<ny;++y)
      z_sbuf[band][YDIR][XDIR][ZDIR][l][x0*ny+y]=x_rbuf[band][YDIR][XDIR][ZDIR][l][z0*ny+y];
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
  for(int y=yl;y<=yh;++y){
    for(int x=xl;x<=xh;++x){
      for(int z=zl;z<=zh;++z){
        J_band[y][x][z]+=c_J*I_n[y][x][z];
        Fz[y][x][z]    +=c_z*I_n[y][x][z];
        Fx[y][x][z]    +=c_x*I_n[y][x][z];
        Fy[y][x][z]    +=c_y*I_n[y][x][z];
      }
    }
  }
}

void RTS::get_Tau_and_Iout(GridData &Grid, const RunData &Run, const PhysicsData &Physics, double DZ, float * B_Iout_tab, float ** kap_Iout_tab, double ** I_band, int calc_int){
  
  double ttime=MPI_Wtime(),stime=0.0,atime=0.0;
  
  double **sbuf=d2dim(yl,yh,xl,xh),**rbuf=d2dim(yl,yh,xl,xh);
  
  double N = pow(2,Grid.NDIM);

  for(int y=yl;y<=yh;++y){ // loop over RT grid
    for(int x=xl;x<=xh;++x){
      int off0 = x*next[1]+y*next[2];
      int off1 = x*next[1]+(y+yo)*next[2];
      int off2 = (x+xo)*next[1]+y*next[2];
      int off3 = (x+xo)*next[1]+(y+yo)*next[2];
      for(int z=zl;z<=zh;++z){
        int inode[]={off0+z,off0+z+zo,off2+z,off2+z+zo,off1+z,off1+z+zo,off3+z,off3+z+zo};
        double lgP = 0, lgT = 0, rm = 0;
        
        for(int l=0;l<N;++l){
          lgP+=Grid.pres[inode[l]];
          lgT+=Grid.temp[inode[l]];
          rm +=Grid.U[inode[l]].d;
        }
    
    lgP /= N;
    lgT /= N;
    rm  /= N;

    lgP          = log(lgP);
    lgT          = log(lgT);
    rho[y][x][z] = rm;

    lgT = min(max(lgT, (double) tab_T[0]),(double) tab_T[NT-1]);
    lgP = min(max(lgP, (double) tab_p[0]),(double) tab_p[Np-1]);

    int l=0;
    int m=0;
    if(lgT<tab_T[0])
      l=0;
    else if(lgT>tab_T[NT-1])
      l=NT-2;
    else
      for (l=0; l<=NT-2; l++)
        if ((lgT >= tab_T[l]) && (lgT <= tab_T[l+1]))
          break;

    if(lgP<tab_p[0])
      m=0;
    else if(lgP>tab_p[Np-1])
      m=Np-2;
    else
      for (m=0; m<=Np-2; m++)
        if ((lgP >= tab_p[m]) && (lgP <= tab_p[m+1]))
          break;

    double xt = (lgT-tab_T[l])*invT_tab[l];
    double xp = (lgP-tab_p[m])*invP_tab[m];

    // Interpolate for kappa and B
    B[y][x][z]=exp(xt*B_Iout_tab[l+1]+(1.-xt)*B_Iout_tab[l]);
    
    kap[y][x][z] = exp(xt*(xp*kap_Iout_tab[l+1][m+1]+(1.-xp)*kap_Iout_tab[l+1][m])+
        (1.-xt)*(xp*kap_Iout_tab[l][m+1]+(1.-xp)*kap_Iout_tab[l][m]));
      }
      Tau[y][x][zh]=1.0e-12;
      for(int z=zh-1;z>=zl;--z){
        double k0=kap[y][x][z],r0=rho[y][x][z],k_upw=kap[y][x][z+1],r_upw=rho[y][x][z+1];
        double delta_tau=DZ*((k0*r0+k_upw*r_upw)*inv3+(k0*r_upw+k_upw*r0)*inv6);
        Tau[y][x][z]=Tau[y][x][z+1]+delta_tau;
      }
    }
  }

  for(int y=yl;y<=yh;++y){ // loop over MHD grid
    for(int x=xl;x<=xh;++x){
      rbuf[y][x]=0.e0;
      sbuf[y][x]=Tau[y][x][zl];
    }
  }

  double ctime=MPI_Wtime(); 
  MPI_Scan(sbuf[yl]+xl,rbuf[yl]+xl,nx*ny,MPI_DOUBLE,MPI_SUM,comm_col[lrank[2]][lrank[1]]);
  stime+=MPI_Wtime()-ctime;
  for(int y=yl;y<=yh;++y){ // loop over MHD grid
    for(int x=xl;x<=xh;++x){
      rbuf[y][x]-=Tau[y][x][zl];
      for(int z=zl;z<=zh;++z){
        Tau[y][x][z]+=rbuf[y][x];
      }
    }
  }

  if (calc_int){
  //  Outgoing Intensity at top (Long Characteristics)
  for(int y=yl;y<=yh;++y){ // loop over MHD grid
    for(int x=xl;x<=xh;++x){
      rbuf[y][x]=0.0;
      sbuf[y][x]=0.0;
      for(int z=zl+1;z<=zh;++z){
    double Ss1 = B[y][x][z];
    double Ss2 = B[y][x][z-1];
    double delta_tau=Tau[y][x][z-1]-Tau[y][x][z];
    if(delta_tau>dtau_min){
      double edt=exp(-delta_tau);
      double c1=(1.0-edt)/delta_tau;
      sbuf[y][x]+=(Ss1*(1.0-c1)+Ss2*(c1-edt))*exp(-Tau[y][x][z]);
    }else{
      sbuf[y][x]+=0.5*delta_tau*(Ss1+Ss2);
    }  
      }          
    }
  }
  ctime=MPI_Wtime(); 
  MPI_Allreduce(sbuf[yl]+xl,rbuf[yl]+xl,nx*ny,MPI_DOUBLE,MPI_SUM,comm_col[lrank[2]][lrank[1]]);
  atime+=MPI_Wtime()-ctime;
  for(int y=yl;y<=yh;++y){
    for(int x=xl;x<=xh;++x){
      I_band[y][x]+=rbuf[y][x];
    }
  }
  }
  del_d2dim(sbuf,yl,yh,xl,xh);
  del_d2dim(rbuf,yl,yh,xl,xh);
  //
  ttime=MPI_Wtime()-ttime;
  if((myrank==0)&&(verbose>2)) printf("tau5000 time: %f %f %f %f \n",ttime,stime,atime,(stime+atime)/ttime);
  
}    
void RTS::tauscale_qrad(int band, double DX,double DY,double DZ, double *** Ss){

  double ttime=MPI_Wtime(),stime=0.0,atime=0.0;
  double idx=1.0/DX,idy=1.0/DY,idz=1.0/DZ;

  if(NDIM==1) {
    idx=0.;
    idy=0.;
  }
  if(NDIM==2) {
    idy=0.0;
  }

  double **sbuf=d2dim(yl,yh,xl,xh),**rbuf=d2dim(yl,yh,xl,xh);
  double ****Qtemp=d4dim(yl,yh,xl,xh,zl,zh,0,1);
  
  for(int y=yl;y<=yh;++y){ // loop over RT grid
    for(int x=xl;x<=xh;++x){
      Tau[y][x][zh]=1.0e-12;
      for(int z=zh-1;z>=zl;--z){
        double k0=kap[y][x][z],r0=rho[y][x][z],k_upw=kap[y][x][z+1],r_upw=rho[y][x][z+1];
        double delta_tau=DZ*((k0*r0+k_upw*r_upw)*inv3+(k0*r_upw+k_upw*r0)*inv6);
        Tau[y][x][z]=Tau[y][x][z+1]+delta_tau;
      }
    }
  }

  for(int y=yl;y<=yh;++y){ // loop over MHD grid
    for(int x=xl;x<=xh;++x){
      rbuf[y][x]=0.e0;
      sbuf[y][x]=Tau[y][x][zl];
    }
  }

  double ctime=MPI_Wtime(); 
  MPI_Scan(sbuf[yl]+xl,rbuf[yl]+xl,nx*ny,MPI_DOUBLE,MPI_SUM,comm_col[lrank[2]][lrank[1]]);
  stime+=MPI_Wtime()-ctime;
  for(int y=yl;y<=yh;++y){ // loop over MHD grid
    for(int x=xl;x<=xh;++x){
      rbuf[y][x]-=Tau[y][x][zl];
      for(int z=zl;z<=zh;++z){
        Tau[y][x][z]+=rbuf[y][x];
      }
    }
  }
  if (need_I){
    //  Outgoing Intensity at top (Long Characteristics)
    for(int y=yl;y<=yh;++y){ // loop over MHD grid
      for(int x=xl;x<=xh;++x){
        rbuf[y][x]=0.0;
        sbuf[y][x]=0.0;
        for(int z=zl+1;z<=zh;++z){
          double Ss1 = Ss[y][x][z];
          double Ss2 = Ss[y][x][z-1];
          double delta_tau=Tau[y][x][z-1]-Tau[y][x][z];
          if(delta_tau>dtau_min){
            double edt=exp(-delta_tau);
            double c1=(1.0-edt)/delta_tau;
            sbuf[y][x]+=(Ss1*(1.0-c1)+Ss2*(c1-edt))*exp(-Tau[y][x][z]);
          }else{
            sbuf[y][x]+=0.5*delta_tau*(Ss1+Ss2);
          }  
        }          
      }
    }

    ctime=MPI_Wtime(); 
    MPI_Allreduce(sbuf[yl]+xl,rbuf[yl]+xl,nx*ny,MPI_DOUBLE,MPI_SUM,comm_col[lrank[2]][lrank[1]]);
    atime+=MPI_Wtime()-ctime;
    for(int y=yl;y<=yh;++y){
      for(int x=xl;x<=xh;++x){
        I_o[y][x]+=rbuf[y][x];
      }
    }
  }
//  radiative energy imbalance
  for(int y=yl;y<=yh;++y){
    for(int x=xl;x<=xh;++x){
      for(int z=zl;z<=zh;++z){
        I_n[y][x][z]=kap[y][x][z]*rho[y][x][z]*(J_band[y][x][z]-Ss[y][x][z]);
        St[y][x][z] +=Ss[y][x][z];
        Jt[y][x][z] +=J_band[y][x][z];
      }
    }
  }
// 
  double inv_tau_0=1.0e1;
  for(int y=yl;y<=yh-yo;++y){
    for(int x=xl;x<=xh-xo;++x){
      for(int z=zl;z<=zh-zo;++z){
        double qf1=((Fz[y   ][x   ][z+zo]+Fz[y   ][x+xo][z+zo]+Fz[y+yo][x   ][z+zo]+Fz[y+yo][x+xo][z+zo])-
        (Fz[y   ][x   ][z   ]+Fz[y   ][x+xo][z   ]+Fz[y+yo][x   ][z   ]+Fz[y+yo][x+xo][z   ]))*idz+
        ((Fx[y   ][x+xo][z   ]+Fx[y   ][x+xo][z+zo]+Fx[y+yo][x+xo][z   ]+Fx[y+yo][x+xo][z+zo])-
        (Fx[y   ][x   ][z   ]+Fx[y   ][x   ][z+zo]+Fx[y+yo][x   ][z   ]+Fx[y+yo][x   ][z+zo]))*idx+
        ((Fy[y+yo][x   ][z   ]+Fy[y+yo][x   ][z+zo]+Fy[y+yo][x+xo][z   ]+Fy[y+yo][x+xo][z+zo])-
        (Fy[y   ][x   ][z   ]+Fy[y   ][x   ][z+zo]+Fy[y   ][x+xo][z   ]+Fy[y   ][x+xo][z+zo]))*idy;
        qf1*=-0.25e0; 
        double qj1 = I_n[y   ][x][z]+I_n[y   ][x][z+zo]+I_n[y   ][x+xo][z]+I_n[y   ][x+xo][z+zo]
        +I_n[y+yo][x][z]+I_n[y+yo][x][z+zo]+I_n[y+yo][x+xo][z]+I_n[y+yo][x+xo][z+zo];
        qj1*=0.5*PI; 

        double tau_local=tau(z+zo,x+xo,y+yo);
        double weight=exp(-tau_local*inv_tau_0);
        Qtemp[y][x][z+zo][0]=qf1;
        Qtemp[y][x][z+zo][1]=qj1;

        Qt[y+yo][x+xo][z+1]+=weight*qj1+(1.0-weight)*qf1;
      }
    }
  }

   if (save_col){
     for (int y=col_bnd[2];y<=col_bnd[3];++y){
       for (int x=col_bnd[0];x<=col_bnd[1];++x){
         for (int z=zl+zo;z<=zh;++z){
           Col_out[band][z-2*zo][0] += J_band[y][x][z]*avg_col;
           Col_out[band][z-2*zo][1] += Ss[y][x][z]*avg_col;
           Col_out[band][z-2*zo][2] += kap[y][x][z]*avg_col;
           Col_out[band][z-2*zo][3] += abn[y][x][z]*avg_col;
           Col_out[band][z-2*zo][4] += sig[y][x][z]*avg_col;
           Col_out[band][z-2*zo][5] += B[y][x][z]*avg_col;
           Col_out[band][z-2*zo][6] += Tau[y][x][z]*avg_col;
           Col_out[band][z-2*zo][7] += Qtemp[y][x][z][0]*avg_col;
           Col_out[band][z-2*zo][8] += Qtemp[y][x][z][1]*avg_col;
         }
       }
     }
   }

  del_d4dim(Qtemp,yl,yh,xl,xh,zl,zh,0,1);
  del_d2dim(sbuf,yl,yh,xl,xh);
  del_d2dim(rbuf,yl,yh,xl,xh);
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
