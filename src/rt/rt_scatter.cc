#include <string.h>
#include <cmath>
#include <fstream>
#include <errno.h>
#include "mem.h"
#include "rt.h"
#include "rt_scatter.h"
#include "comm_split.H"
#include "ACCH.h"


RTS_SCATTER::~RTS_SCATTER(void)
{
  // source function for each band
  if (!fullodf){
    ACCH::Free2D<double>(S, Nbands, nx*ny*nz);
  } else {
    ACCH::Free2D<double>(S, 1, nx*ny*nz);
  }

  del_d2dim(J_havg,0,Nbands-1,zl,zh);
  
  // lambda_star for NLTE scattering problem
  del_d3dim(lambda_star,yl,yh,xl,xh,zl,zh);
  
  del_i3dim(Tau_ind,yl,yh,xl,xh,zl,zh);
  del_d3dim(lgTau,yl,yh,xl,xh,zl,zh);

  // and the nel
  ACCH::Free(ne, nx*ny*nz*sizeof(double));

  if ((fullodf)&&(scatter==0)){
    del_f3dim(acont_pT,0,Nlam-1,0,NT-1,0,Np-1);
    del_f3dim(kcont_pT,0,Nlam-1,0,NT-1,0,Np-1);
  }
  if ((!fullodf)&&(scatter==0)){
    del_f3dim(sig_tab,0,Nbands-1,0,NT-1,0,Np-1);
    del_f3dim(abn_tab,0,Nbands-1,0,NT-1,0,Np-1);

    del_f2dim(abn_pp_tab,0,Nbands-1,0,Npp-1);
    del_f2dim(kap_pp_tab,0,Nbands-1,0,Npp-1);
    del_f2dim(sig_pp_tab,0,Nbands-1,0,Npp-1);
  }

  del_d1dim(tau_pp_tab,0,Npp-1);
  del_d1dim(invtau_pp_tab,0,Npp-1);
}

RTS_SCATTER::RTS_SCATTER(GridData & Grid,RunData & Run,PhysicsData & Physics):RTS(Grid,Run,Physics){
  
  // electron number
  ne = (double*) ACCH::Malloc(nx*ny*nz*sizeof(double));

  // The lambda_star
  lambda_star = d3dim(yl,yh,xl,xh,zl,zh);
  memset(lambda_star[yl][xl]+zl,0,nx*ny*nz*sizeof(double));
  
  Tau_ind = i3dim(yl,yh,xl,xh,zl,zh);

  lgTau = d3dim(yl,yh,xl,xh,zl,zh);

 /*Source function for each band:
   * if full ODF, then we use 1 band, and recalculate source function from a horizontally
   * averaged J */

  if (!fullodf){
    S = ACCH::Malloc2D<double>(Nbands, nx*ny*nz);
    memset(S[0],0,nx*ny*nz*Nbands*sizeof(double)); // MHD grid
  } else {
    S = ACCH::Malloc2D<double>(1, nx*ny*nz);
    memset(S[0],0,nx*ny*nz*sizeof(double)); // MHD grid
  }

  J_havg=d2dim(0,Nbands-1,zl,zh);
  memset(J_havg[0]+zl,0,nz*Nbands*sizeof(double));  

  for (int band=0;band<Nbands;band++)
    for (int z = zl;z<=zh;z++)
      J_havg[band][z] = 1.0e5;

  ini_flag = i1dim(0,Nbands-1);

  for (int band=0;band<Nbands;band++)
    ini_flag[band] = 0;
}

double RTS_SCATTER::wrapper(int rt_upd,GridData &Grid,RunData &Run,const PhysicsData &Physics){

  // Van regemorter function P(deltaE/kt)
  
  double VR_deltaE[14] ={0.0025,0.005,0.01,0.02,0.04,0.1,0.2,0.4,1.0,2.0,4.0,10.0,15.0,20.0};
  double VR_P[14] = {1.4925147328502524,1.3014386981388302, 1.160,0.956,0.758,0.493,
                    0.331,0.209,0.100,0.063,0.040,0.023,0.017041126723312636,0.014758048651498613};

  memset(I_band,0,nx*ny*sizeof(double));
  memset(I_o,0,nx*ny*sizeof(double));

  double DX=Grid.dx[1],DZ=Grid.dx[0],DY=Grid.dx[2];
  int cont_bin = Physics.rt[i_rt_iout];
  const double Temp_TR = Physics.rt[i_rt_tr_tem];
  const double Pres_TR = Physics.rt[i_rt_tr_pre];
  
  fail_count = 0 ;

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
          I_o[(y-yl)*nx+(x-xl)] +=I_band[(y-yl)*nx+(x-xl)];
      }
    }
    // If cont_bin = 0 and need_I = 1 we want the output intensity to be the 0th continuum bin
    if (cont_bin==0){
      get_Tau_and_Iout(Grid, Run, Physics,DZ,B_tab[0],kap_tab[0],I_band,need_I);
      for (int y=yl;y<=yh;y++)
        for (int x=xl;x<=xh;x++)
          I_o[(y-yl)*nx+(x-xl)] +=I_band[(y-yl)*nx+(x-xl)];
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
            I_o[(y-yl)*nx+(x-xl)] +=I_band[(y-yl)*nx+(x-xl)];
    }

    calc_Qtot_and_Tau(Grid, Run, Physics);
    return dt_rad;
  }
  
  memset(St,0,nx*ny*nz*sizeof(double)); // MHD grid
  memset(Jt,0,nx*ny*nz*sizeof(double)); // MHD grid
  memset(Qt[0][0],0,(nx-xo)*(ny-yo)*(nz-zo)*sizeof(double)); // MHD grid
    
// *****************************************************************
// *        interpolate opacity and Planck function (B)            *
// *****************************************************************
  cState *U=Grid.U;

  double N = pow(2,NDIM);
  get_Tau_and_Iout(Grid, Run, Physics,DZ,B_5000_tab,kap_5000_tab,I_band,0);

  for(int y=yl;y<=yh;++y){
    for(int x=xl;x<=xh;++x){
      int off0 = x*next[1]+y*next[2];
      int off1 = x*next[1]+(y+yo)*next[2];
      int off2 = (x+xo)*next[1]+y*next[2];
      int off3 = (x+xo)*next[1]+(y+yo)*next[2];
      for(int z=zl;z<=zh;++z){
        int inode[]={off0+z,off0+z+zo,off2+z,off2+z+zo,off1+z,off1+z+zo,off3+z,off3+z+zo};
        double Tm=0.0,pm=0.0,rm=0.0,nm=0.0;

        for(int l=0;l<N;++l){
          Tm+=Grid.temp[inode[l]];
          pm+=Grid.pres[inode[l]];
          rm+=U[inode[l]].d;
          nm+=Grid.ne[inode[l]];
        }
 
        Tm /= N;
        pm /= N;
        rm /= N;
        nm /= N;

        //disbale RT if Temp > 2e4 above the photosphere
        double pswitch = min(max(pm-Pres_TR,0.0),1.0);
        double tswitch = min(max(Temp_TR-Tm,0.0),1.0);
        tr_switch[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)] = (int) max(pswitch,tswitch);

        ne[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)] = nm;
        lgTe[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)]=log(Tm);
        lgPe[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)]=log(pm);
        rho[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)] = rm;

        lgTe[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)] = min(max(lgTe[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)],(double) tab_T[0]),(double) tab_T[NT-1]);
        lgPe[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)] = min(max(lgPe[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)],(double) tab_p[0]),(double) tab_p[Np-1]);
        
        // Search table once for all RT bands
        int l=0;
        int m=0;
        if(lgTe[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)]<tab_T[0])
          l=0;
        else if(lgTe[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)]>tab_T[NT-1])
          l=NT-2;
        else
          for (l=0; l<=NT-2; l++)
            if ((lgTe[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)] >= tab_T[l]) && (lgTe[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)] <= tab_T[l+1]))
              break;

        if(lgPe[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)]<tab_p[0])
          m=0;
        else if(lgPe[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)]>tab_p[Np-1])
          m=Np-2;
        else
          for (m=0; m<=Np-2; m++)
             if ((lgPe[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)] >= tab_p[m]) && (lgPe[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)] <= tab_p[m+1]))
               break;

        T_ind[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)] = l;
        P_ind[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)] = m;

        lgTau[y][x][z] = min(max(log(Tau[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)]),(double) tau_pp_tab[0]),(double) tau_pp_tab[Npp-1]);

        int t=0;
        if(lgTau[y][x][z] < tau_pp_tab[0])
          t=0;
        else if (lgTau[y][x][z] > tau_pp_tab[Npp-1])
          t=Npp-2;
        else
          for (t=0;t<=Npp-2;t++)
            if ((lgTau[y][x][z] >= tau_pp_tab[t]) && (lgTau[y][x][z] <= tau_pp_tab[t+1]))
              break;

        Tau_ind[y][x][z] = t;
      }
    }
  }
  
  // Begin RT loop, work band by band,
  
  total_iter = 0.0;
  for(int band=Nbands-1;band>=0;--band){
  
    // If full ODF only 1-band S function
    int band_S = band*(1-fullodf);
    
    nu_ind = 0;
    bin_ind = 0;

    if(fullodf){
      nu_ind = band/Nbin;
      bin_ind = band%Nbin;
    }

    if((myrank==0) && (verbose>1)){
      fprintf(stdout,"|------------------------------------------------------------|\n");
      fprintf(stdout,"rt running for band %i of %i \n",band+1,Nbands);
      if (fullodf) 
        fprintf(stdout,"nu band = %i, %e and bin %i \n", nu_ind, nu_tab[nu_ind],bin_ind);
    }
   
    /* Four Cases:
     * 1) Binned opacities + constant photon destruction probability
     * 2) Binned opacities, including seperate absorption and scattering (skartlien-like)
     * 3) full odf + constant photon destruction probability
     * 4) full odf w/ seperate absorption opacity
     */
  

    if (!fullodf){
      if(eps_const>0){
        for(int y=yl;y<=yh;++y){
          for(int x=xl;x<=xh;++x){
            for(int z=zl;z<=zh;++z){
              int l = T_ind[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)];
              int m = P_ind[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)];
        
              double xt = (lgTe[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)]-tab_T[l])*invT_tab[l];
              double xp = (lgPe[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)]-tab_p[m])*invP_tab[m];
        
              // Interpolate for kappa and B
              B[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)]=exp(xt*B_tab[band][l+1]+(1.-xt)*B_tab[band][l]);
        
              double kap1 = exp(xt*(xp*kap_tab[band][l+1][m+1]+(1.-xp)*kap_tab[band][l+1][m])+
              (1.-xt)*(xp*kap_tab[band][l][m+1]+(1.-xp)*kap_tab[band][l][m]));
              
              kap[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)] = kap1;
              abn[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)] = eps_const;
              sig[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)] = (1.0-eps_const);

              /* TR switch turn off kappa and B */
              B[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)] *= tr_switch[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)];
              sig[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)] *= tr_switch[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)];
              abn[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)] *= tr_switch[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)];
              kap[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)] *= tr_switch[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)];
            }
          }
        }
      } else {
        for(int y=yl;y<=yh;++y){
          for(int x=xl;x<=xh;++x){
            for(int z=zl;z<=zh;++z){
              int l = T_ind[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)];
              int m = P_ind[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)];

              int t = Tau_ind[y][x][z];

              double xt   = (lgTe[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)]-tab_T[l])*invT_tab[l];
              double xp   = (lgPe[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)]-tab_p[m])*invP_tab[m];
             
              double xtau = (lgTau[y][x][z] - tau_pp_tab[t])*invtau_pp_tab[t];

              // Interpolate for kappa and B
              B[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)]=exp(xt*B_tab[band][l+1]+(1.-xt)*B_tab[band][l]);
        
              double kap1 = exp(xt*(xp*kap_tab[band][l+1][m+1]+(1.0-xp)*kap_tab[band][l+1][m])+
              (1.-xt)*(xp*kap_tab[band][l][m+1]+(1.-xp)*kap_tab[band][l][m]));

              double sig1 = exp(xt*(xp*sig_tab[band][l+1][m+1]+(1.0-xp)*sig_tab[band][l+1][m])+
              (1.-xt)*(xp*sig_tab[band][l][m+1]+(1.-xp)*sig_tab[band][l][m]));

              double abn1 = exp(xt*(xp*abn_tab[band][l+1][m+1]+(1.0-xp)*abn_tab[band][l+1][m])+
              (1.-xt)*(xp*abn_tab[band][l][m+1]+(1.-xp)*abn_tab[band][l][m]));
             
              double sig2 = xtau*sig_pp_tab[band][t+1] + (1.0-xtau)*sig_pp_tab[band][t];
              double kap2 = xtau*kap_pp_tab[band][t+1] + (1.0-xtau)*kap_pp_tab[band][t];
              //double abn2 = xtau*abn_pp_tab[band][t+1] + (1.0-xtau)*abn_pp_tab[band][t];
              
              double weight=exp(-kap1*30.0);

              abn1 = abn1; //*(1.0-weight) + abn2*weight;
              kap1 = kap1*(1.0-weight) + kap2*weight;
              sig1 = sig1*(1.0-weight) + sig2*weight;

              kap[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)] = abn1+sig1;
              abn[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)] = abn1/(abn1+sig1);
              sig[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)] = sig1/(abn1+sig1);

              /* TR switch turn off kappa and B */

              B[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)] *= tr_switch[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)];
              sig[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)] *= tr_switch[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)];
              abn[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)] *= tr_switch[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)];
              kap[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)] *= tr_switch[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)];
            }
          }
        }
      }
      if (ini_flag[band]){
        for(int y=yl;y<=yh;++y)
          for(int x=xl;x<=xh;++x){
            S[band_S][(y-yl)*nx*nz+(x-xl)*nz] = B[(y-yl)*nx*nz+(x-xl)*nz];
            S[band_S][(y-yl)*nx*nz+(x-xl)*nz+(nz-1)] = 0.0;
          }
        ini_flag = 0;
      }
    } else {
      if(eps_const>0){
        for(int y=yl;y<=yh;++y){
          for(int x=xl;x<=xh;++x){
            for(int z=zl;z<=zh;++z){
            
              int l = T_ind[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)];
              int m = P_ind[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)];
      
              double xt = (lgTe[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)]-tab_T[l])*invT_tab[l];
              double xp = (lgPe[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)]-tab_p[m])*invP_tab[m];
      
              // Interpolate for kappa and B
              B[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)]=exp(xt*B_tab[band][l+1]+(1.-xt)*B_tab[band][l]);
      
              double kap1 = exp(xt*(xp*kap_tab[band][l+1][m+1]+(1.-xp)*kap_tab[band][l+1][m])+
              (1.-xt)*(xp*kap_tab[band][l][m+1]+(1.-xp)*kap_tab[band][l][m]));
              
              double acont = xt*(xp*acont_pT[nu_ind][l+1][m+1]+(1.-xp)*acont_pT[nu_ind][l+1][m])+
                (1.-xt)*(xp*acont_pT[nu_ind][l][m+1]+(1.-xp)*acont_pT[nu_ind][l][m]);

              kap1 = kap1+acont;
              kap[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)] = kap1;
              abn[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)] = eps_const;
              sig[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)] = (1-eps_const);
              
              /* TR switch turn off kappa and B */
              
              B[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)] *= tr_switch[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)];
              sig[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)] *= tr_switch[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)];
              abn[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)] *= tr_switch[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)];
              kap[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)] *= tr_switch[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)];
              
              // For a full ODF only one band in S to save memory. Instead, use an averaged J as a first guess.

              S[band_S][(y-yl)*nx*nz+(x-xl)*nz+(z-zl)]=B[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)]*abn[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)]+sig[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)]*J_havg[band][z];
              J_havg[band][z] = 0.0;
            }
          }
        }    
      } else {
        for(int y=yl;y<=yh;++y){
          for(int x=xl;x<=xh;++x){
            for(int z=zl;z<=zh;++z){
        
              double nel = ne[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)];
              double TT = exp(lgTe[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)]);
              
              int l = T_ind[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)];
              int m = P_ind[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)];
      
              double xt = (lgTe[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)]-tab_T[l])*invT_tab[l];
              double xp = (lgPe[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)]-tab_p[m])*invP_tab[m];
      
              // Interpolate for kappa and B
              B[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)]=exp(xt*B_tab[band][l+1]+(1.-xt)*B_tab[band][l]);
      
              double kap1 = exp(xt*(xp*kap_tab[band][l+1][m+1]+(1.-xp)*kap_tab[band][l+1][m])+
              (1.-xt)*(xp*kap_tab[band][l][m+1]+(1.-xp)*kap_tab[band][l][m]));
              
              double acont = xt*(xp*acont_pT[nu_ind][l+1][m+1]+(1.-xp)*acont_pT[nu_ind][l+1][m])+
                (1.-xt)*(xp*acont_pT[nu_ind][l][m+1]+(1.-xp)*acont_pT[nu_ind][l][m]);

              double kcont = xt*(xp*kcont_pT[nu_ind][l+1][m+1]+(1.-xp)*kcont_pT[nu_ind][l+1][m])+
              (1.-xt)*(xp*kcont_pT[nu_ind][l][m+1]+(1.-xp)*kcont_pT[nu_ind][l][m]);

              double deltaE = H*nu_tab[nu_ind]/k_B/TT;
                
              double interp = 0.0;
              int ee = 0;
              if (deltaE <= VR_deltaE[0]){
                interp = pow(3.0,0.5)/(2.0*PI)*(-0.57722-log(deltaE));
              } else if (deltaE >= VR_deltaE[13]){
                interp = 0.066*pow(deltaE,-0.5);
              } else {
                for (ee=0;ee<12;ee++)
                  if (deltaE <= VR_deltaE[ee+1])
                    break;
                interp = VR_P[ee]+ (VR_P[ee+1]-VR_P[ee])/(VR_deltaE[ee+1]-VR_deltaE[ee])*(deltaE-VR_deltaE[ee]);
              }

              double ConA=20.6*pow(Cvac/nu_tab[nu_ind],3)*nel*pow(TT,-0.5)*interp;

              double VR =  ConA/(ConA+1.0);
              
              double abs1 = kap1*VR+kcont;
              
              kap1 = kap1+acont;

              kap[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)] = kap1;
              abn[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)] = abs1/kap1;
              sig[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)] = (kap1-abs1)/kap1; // True for full ODF

              /* TR switch turn off kappa and B */
              B[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)] *= tr_switch[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)];
              sig[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)] *= tr_switch[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)];
              abn[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)] *= tr_switch[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)];
              kap[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)] *= tr_switch[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)];
              
              S[band_S][(y-yl)*nx*nz+(x-xl)*nz+(z-zl)]=B[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)]*abn[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)]+sig[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)]*J_havg[band][z];
              J_havg[band][z] = 0.0;
            }
          }
        }
      }
    }

// *****************************************************************
// *    update diffusion-approx. boundary condition at bottom      *
// *****************************************************************
  for(int YDIR=FWD;YDIR<=BWD;++YDIR)
    for(int XDIR=RIGHT;XDIR<=LEFT;++XDIR)
      for(int l=0;l<NMU;++l){
        if(isgbeg[0]==1)
          for(int y=0;y<ny;++y)
            for(int x=0;x<nx;++x)
              z_rbuf[band][YDIR][XDIR][UP][l][x*ny+y]=B[y*nx*nz+x*nz];

        if(isgend[0]==1)
          for(int y=0;y<ny;++y)
            for(int x=0;x<nx;++x)
              z_rbuf[band][YDIR][XDIR][DOWN][l][x*ny+y]=0.0e0;
      }
// *****************************************************************
// *  get lambda operator                                          *
// *****************************************************************
  
  get_lambdastar();
  
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
        gFr_mean[band]+=Fz[y*nx*nz+x*nz+nz-1] +
                        Fz[(y+yo)*nx*nz+x*nz+nz-1] +
                        Fz[y*nx*nz+(x+xo)*nz+nz-1] +
                        Fz[(y+yo)*nx*nz+(x+xo)*nz+nz-1];
    gFr_mean[band]*=0.25;
  }

  // If I am saving 5000A intensity then wipe after last bin. if I am saving the continuum bin then wipe the second last bin.

  if (((cont_bin==1)&&(band==0))||((cont_bin==0)&&(band==1)))
    for(int y=yl;y<=yh;++y)
      for(int x=xl;x<=xh;++x)
        I_o[(y-yl)*nx+(x-xl)] = 0.0;

  }// end loop over bands

  MPI_Allreduce(&gFr_mean[0],&Fr_mean[0],Nbands,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  F_o = 0.0;

  for (int band=0;band<Nbands;++band){
    Fr_mean[band]/=(Grid.gsize[1]*Grid.gsize[2]);
    F_o+=Fr_mean[band];
  }
  
  if (N5000){
    int I5000_out = 0;
    if ((cont_bin==1)&&(need_I==1))
      I5000_out = 1;
    
    get_Tau_and_Iout(Grid, Run, Physics,DZ,B_5000_tab,kap_5000_tab,I_band,I5000_out);
   
    if ((cont_bin==1)&&(need_I==1))
      for (int y=yl;y<=yh;y++)
        for (int x=xl;x<=xh;x++)
          I_o[(y-yl)*nx+(x-xl)] +=I_band[(y-yl)*nx+(x-xl)];
  }
  
  calc_Qtot_and_Tau(Grid, Run, Physics);

  if (myrank==0)
    fprintf(stdout,"NLTE_scattering Iter %i bands failed to converge, total iterations over all bands %i \n",fail_count,total_iter);
  
  if (Run.NeedsSlice())
    save_1D_avg(Run.path_2D,Run.globiter,Run.time);  

  return dt_rad;
}

void RTS_SCATTER::driver(double DZ, double DX, double DY, int band) {
  double etime=0.0,atime=0.0,cmp_time1=0.0,cmp_time2=0.0,buf_time=0.0,err_time=0.0,flx_time=0.0,tau_time=0.0,ttime = 0.0;
 
  const int stride[2]={nx*nz,nz}; 
  int stepvec[3][4][3] = { {{1,0,0},{1,0,1},{1,1,0},{1,1,1}},
               {{0,1,0},{1,1,0},{0,1,1},{1,1,1}},
               {{0,0,1},{0,1,1},{1,0,1},{1,1,1}} };

  int nlte_iter = 0;
  int max_nlte_iter = 100000;

  double S_tol = 1.0e-3;

  // J= sum of rays, to get under S_tol rt_tol should ideally be
  // S_tol/Nrays = 8*NMU
  
  double rt_tol = min(threshold,S_tol);
  double mean_S = 1.0e7;
  double mean_delta_S = 0.0;
  
  double min_lambda_star = 1.0;
  double max_lambda_star = 0.0;
  
  double S_all_err = 1.0e12;
  double S_err = 1.0e8;

  int band_S = band*(1-fullodf);
  // While source error is greater than tolerance and iterations are less than maxiterations then continue
  // Do a couple of iterations on the first iter.
  memset(I_n,0,nx*ny*nz*sizeof(double)); // is there a better way?
  
            
  while((S_all_err > S_tol)&&(nlte_iter<=max_nlte_iter)){ 

    mean_S = 0.0;
    mean_delta_S = 0.0;
    
    min_lambda_star = 1.0;
    max_lambda_star = 0.0;

    memset(Fz,0,nx*ny*nz*sizeof(double));
    memset(Fx,0,nx*ny*nz*sizeof(double));
    memset(J_band,0,nx*ny*nz*sizeof(double));
    memset(Fy,0,nx*ny*nz*sizeof(double));

    if (nlte_iter == max_nlte_iter){
      fail_count+=1;
      break;
    }

    // Minimum source function to speed convergence
    double S_min = 1.0e-40;
    double maxerr_up=0.0,maxerr_down=0.0;
  // loop over octants & determination of loop direction
    double itavg=0.0;
    double aravg=0.0;
    ttime = MPI_Wtime();
  // main loop
    for(int ZDIR=UP;ZDIR<=DOWN;++ZDIR){
      int zi_i=(ZDIR==UP)?zl+1:zh-1,zi_f=(ZDIR==UP)?zh:zl,zstep=(ZDIR==UP)?1:-1;
      for(int YDIR=FWD;YDIR<=BWD;++YDIR){
        int yi_i=(YDIR==FWD)?yl+1:yh-1,yi_f=(YDIR==FWD)?yh:yl,ystep=(YDIR==FWD)?1:-1;
        for(int XDIR=RIGHT;XDIR<=LEFT;++XDIR){
          int xi_i=(XDIR==RIGHT)?xl+1:xh-1,xi_f=(XDIR==RIGHT)?xh:xl,xstep=(XDIR==RIGHT)?1:-1;
          for(int l=0;l<NMU;++l){
            double I_min = max(1.0,rt_tol*Fr_mean[band]/(NMU*pow(2,NDIM)));
            double c[]={a_00[ibase[l]][l],a_01[ibase[l]][l],a_10[ibase[l]][l],a_11[ibase[l]][l]};
            for(int m=0;m<=3;++m){
              ixstep[m]=stepvec[ibase[l]][m][0]*xstep;
              iystep[m]=stepvec[ibase[l]][m][1]*ystep;
              izstep[m]=stepvec[ibase[l]][m][2]*zstep;
            } 
            double stime=MPI_Wtime();
            interpol(zi_i,zi_f,zstep,xi_i,xi_f,xstep,yi_i,yi_f,ystep,l,coeff,S[band_S]);
            cmp_time1+=MPI_Wtime()-stime;
            stime=MPI_Wtime(); 
            int rt_iter=0;
            double gmaxerr=1.0E10;
            // Iteration start
            while(gmaxerr>=rt_tol){
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
              double *ii=I_n;
              if(NDIM==3){
                for(int i=0;i<4;i++) 
                  off[i]=iystep[i]*stride[0]+ixstep[i]*stride[1]+izstep[i];
                for(int yi=yi_i;yi!=yi_f+ystep;yi=yi+ystep){
                  for(int xi=xi_i;xi!=xi_f+xstep;xi=xi+xstep){
                    int xyoff=(yi-yl)*stride[0]+(xi-xl)*stride[1];
                    for(int zi=zi_i;zi!=zi_f+zstep;zi=zi+zstep){
                      int ind=xyoff+zi-zl;
                      double I_upw=c[0]*ii[ind-off[0]]+c[1]*ii[ind-off[1]]+c[2]*ii[ind-off[2]]+c[3]*ii[ind-off[3]];
                      ii[ind]=I_upw*coeff[i_nu]+coeff[nx*ny*nz+i_nu];
                      i_nu+=1;
                    }
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
                    ii[ind]=I_upw*coeff[i_nu]+coeff[i_nu+nx*ny*nz];
                    i_nu+=1;
                  }
                }
              }

              if(NDIM==1){
                for(int i=0;i<4;i++) off[i]=izstep[i];
                for(int zi=zi_i;zi!=zi_f+zstep;zi=zi+zstep){
                  double I_upw=c[0]*ii[zi-off[0]]+c[1]*ii[zi-off[1]]+c[2]*ii[zi-off[2]]+c[3]*ii[zi-off[3]];
                  ii[zi]=I_upw*coeff[i_nu]+coeff[i_nu+nx*ny*nz];
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

              stime = MPI_Wtime();
              double err=error(band,l,ZDIR,XDIR,YDIR,I_min);
              err_time += MPI_Wtime()-stime;
              stime = MPI_Wtime();
              MPI_Allreduce(&err,&gmaxerr,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
              aravg += 1.0;
              atime += MPI_Wtime()-stime;
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

            double sstime=MPI_Wtime();
            flux(l,ZDIR,XDIR,YDIR);
            flx_time += MPI_Wtime()-sstime;
          }// end l
        }// end YDIR
      }// end XDIR
    }// end ZDIR 
    
    if((myrank==0) && (verbose>2)){
      fprintf(stdout,"rt_driver iter : %f %f \n",aravg/(8.0*NMU),itavg/(8.0*NMU));
      fprintf(stdout,"rt_driver error: %e %e \n",maxerr_up,maxerr_down);
    }

    ttime=MPI_Wtime()-ttime;  
    if((myrank==0) && (verbose>2))
      fprintf(stdout,"rt_driver time : %f %f %f %f %f %f %f %f %f %f \n",ttime,cmp_time1,
         cmp_time2,buf_time,err_time,flx_time,tau_time,etime,atime,(etime+atime)/ttime);
    
    // increment for S and J 

    S_err = S_tol*S_tol;
    
    for(int y=yl;y<=yh;y++){
      for(int x=xl;x<=xh;x++){
        for(int z=zl;z<=zh;++z){
          double ls = max(min(lambda_star[y][x][z],0.99999),0.0);
          double S_old = S[band_S][(y-yl)*nx*nz+(x-xl)*nz+(z-zl)];
          double tr = tr_switch[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)];
          min_lambda_star = min(ls,min_lambda_star);
          max_lambda_star = max(ls,max_lambda_star);

          double dS_n = tr*(sig[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)]*J_band[(y-yl)*nx*nz + (x-xl)*nz + (z-zl)] + abn[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)]*B[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)] - S_old)/(1.0-sig[(y-yl)*nx*nz+(x-xl)*nz+(z-zl)]*ls);

          J_havg[band][z]+=J_band[(y-yl)*nx*nz + (x-xl)*nz + (z-zl)]/(nx*ny);

          S[band_S][(y-yl)*nx*nz+(x-xl)*nz+(z-zl)] = S_old+dS_n; //Add correction
          mean_S += S[band_S][(y-yl)*nx*nz+(x-xl)*nz+(z-zl)];

          dS_n = abs(dS_n)/max(S_old,S_min);
          S_err = max(dS_n,S_err);  // For convergence we will want delta_S to be under some small number
          mean_delta_S += dS_n;
        } //end loop
      }
    }

    //Reduce for Error.
    MPI_Allreduce(&S_err,&S_all_err,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);

    if (verbose > 2){
      mean_delta_S /= (nx*ny*nz);
      mean_S /= (nx*ny*nz);

      double all_mean_delta_S = 0;
      double mean_S_all = 0;

      MPI_Allreduce(&mean_S,&mean_S_all,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      MPI_Allreduce(&mean_delta_S,&all_mean_delta_S,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

      mean_S_all/=(cart_sizes[1]*cart_sizes[2]*cart_sizes[0]);
      mean_delta_S/=(cart_sizes[1]*cart_sizes[2]*cart_sizes[0]);
      
      if ((myrank==0)){
        
        fprintf(stdout,"Source correction iteration complete, number of iterations: %i \n",nlte_iter);
        fprintf(stdout,"Final delta / S: %21.15E, mean source %21.15E and maximum error: %21.15E\n", all_mean_delta_S,mean_S_all,S_all_err);
        fprintf(stdout,"min/max lambda_star: %21.15E %21.15E\n", min_lambda_star, max_lambda_star);
      }
    }

    nlte_iter +=1;
    total_iter+=1;
    
  } // end while
  
  if (verbose <= 2){
    mean_delta_S /= (nx*ny*nz);
    mean_S /= (nx*ny*nz);

    // Reduce for Error.
    double all_mean_delta_S = 0;
    double mean_S_all = 0;

    MPI_Allreduce(&mean_S,&mean_S_all,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&mean_delta_S,&all_mean_delta_S,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

    mean_S_all/=(cart_sizes[1]*cart_sizes[2]*cart_sizes[0]);
    mean_delta_S/=(cart_sizes[1]*cart_sizes[2]*cart_sizes[0]);

    if ((myrank==0)){
        
      fprintf(stdout,"Source correction iteration complete, number of iterations: %i \n",nlte_iter);
      fprintf(stdout,"Final delta / S: %21.15E, mean source %21.15E and maximum error: %21.15E\n", all_mean_delta_S,mean_S_all,S_all_err);
      fprintf(stdout,"min/max lambda_star: %21.15E %21.15E\n", min_lambda_star, max_lambda_star);
    }
  }
    
  double stime=MPI_Wtime();
  tauscale_qrad(band,DX,DY,DZ,S[band_S]);
  tau_time+=MPI_Wtime()-stime;

  call_count+=1;
}

void RTS_SCATTER::get_lambdastar(){

  memset(lambda_star[yl][xl]+zl,0,nx*ny*nz*sizeof(double));

  int stepvec[3][4][3] = { {{1,0,0},{1,0,1},{1,1,0},{1,1,1}},
               {{0,1,0},{1,1,0},{0,1,1},{1,1,1}},
               {{0,0,1},{0,1,1},{1,0,1},{1,1,1}} };

  for(int ZDIR=UP;ZDIR<=DOWN;++ZDIR){
    int zi_i=(ZDIR==UP)?zl+1:zh-1,zi_f=(ZDIR==UP)?zh:zl,zstep=(ZDIR==FWD)?1:-1;
    for(int YDIR=FWD;YDIR<=BWD;++YDIR){
      int yi_i=(YDIR==FWD)?yl+1:yh-1,yi_f=(YDIR==FWD)?yh:yl,ystep=(YDIR==FWD)?1:-1;
      for(int XDIR=RIGHT;XDIR<=LEFT;++XDIR){
        int xi_i=(XDIR==RIGHT)?xl+1:xh-1,xi_f=(XDIR==RIGHT)?xh:xl,xstep=(XDIR==RIGHT)?1:-1;
        for(int l=0;l<NMU;++l){
  
          int zmin=(zi_i<zi_f)?zi_i:zi_f;
          int zmax=(zi_i>zi_f)?zi_i:zi_f;
            
          for(int m=0;m<=3;++m){
            ixstep[m]=stepvec[ibase[l]][m][0]*xstep;
            iystep[m]=stepvec[ibase[l]][m][1]*ystep;
            izstep[m]=stepvec[ibase[l]][m][2]*zstep;
          } 
          double expo[zmax+1],dt[zmax+1],r_upw[zmax+1],k_upw[zmax+1],r0[zmax+1],k0[zmax+1];
          double ds3=ds_upw[l]*inv3,ds6=ds_upw[l]*inv6;
          double c[]={a_00[ibase[l]][l],a_01[ibase[l]][l],a_10[ibase[l]][l],a_11[ibase[l]][l]};
          double c_J = 0.125*wmu[l];

          if(NDIM==3){
            for(int yi=yi_i;yi!=yi_f+ystep;yi=yi+ystep){
              for(int xi=xi_i;xi!=xi_f+xstep;xi=xi+xstep){
                for(int zi=zmin;zi<=zmax;++zi){
                  r_upw[zi]=
                    c[0]*rho[(yi-iystep[0]-yl)*nx*nz+(xi-ixstep[0]-xl)*nz+(zi-izstep[0]-zl)]+
                    c[1]*rho[(yi-iystep[1]-yl)*nx*nz+(xi-ixstep[1]-xl)*nz+(zi-izstep[1]-zl)]+
                    c[2]*rho[(yi-iystep[2]-yl)*nx*nz+(xi-ixstep[2]-xl)*nz+(zi-izstep[2]-zl)]+
                    c[3]*rho[(yi-iystep[3]-yl)*nx*nz+(xi-ixstep[3]-xl)*nz+(zi-izstep[3]-zl)];

                    k_upw[zi]=
                    c[0]*kap[(yi-iystep[0]-yl)*nx*nz+(xi-ixstep[0]-xl)*nz+(zi-izstep[0]-zl)]+
                    c[1]*kap[(yi-iystep[1]-yl)*nx*nz+(xi-ixstep[1]-xl)*nz+(zi-izstep[1]-zl)]+
                    c[2]*kap[(yi-iystep[2]-yl)*nx*nz+(xi-ixstep[2]-xl)*nz+(zi-izstep[2]-zl)]+
                    c[3]*kap[(yi-iystep[3]-yl)*nx*nz+(xi-ixstep[3]-xl)*nz+(zi-izstep[3]-zl)];

                  r0[zi]=rho[(yi-yl)*nx*nz+(xi-xl)*nz+(zi-zl)];
                  k0[zi]=kap[(yi-yl)*nx*nz+(xi-xl)*nz+(zi-zl)];
                }
                for(int zi=zmin;zi<=zmax;++zi){
                  dt[zi]=ds3*(k_upw[zi]*r_upw[zi]+k0[zi]*r0[zi])+ds6*(k0[zi]*r_upw[zi]+k_upw[zi]*r0[zi]);
                  expo[zi]=exp(-dt[zi]);
                  double w0=(1.0-expo[zi])/dt[zi];

                  if (dt[zi] > dtau_min)
                    lambda_star[yi][xi][zi] += c_J*(1.0-w0);
                  else
                    lambda_star[yi][xi][zi] += c_J*(0.5*dt[zi]);
                }   
              }
            }
          }

          if(NDIM==1){
             for(int zi=zmin;zi<=zmax;++zi){
               r_upw[zi]=
                 c[0]*rho[(zi-izstep[0]-zl)]+
                 c[1]*rho[(zi-izstep[1]-zl)]+
                 c[2]*rho[(zi-izstep[2]-zl)]+
                 c[3]*rho[(zi-izstep[3]-zl)];

              k_upw[zi]=
              c[0]*kap[(zi-izstep[0]-zl)]+
              c[1]*kap[(zi-izstep[1]-zl)]+
              c[2]*kap[(zi-izstep[2]-zl)]+
              c[3]*kap[(zi-izstep[3]-zl)];

              r0[zi]=rho[zi-zl];
              k0[zi]=kap[zi-zl];
            }

            for(int zi=zmin;zi<=zmax;++zi){
              dt[zi]=ds3*(k_upw[zi]*r_upw[zi]+k0[zi]*r0[zi])+ds6*(k0[zi]*r_upw[zi]+k_upw[zi]*r0[zi]);
              expo[zi]=exp(-dt[zi]);
              double w0=(1.0-expo[zi])/dt[zi];

              if (dt[zi] > dtau_min)
                lambda_star[0][0][zi] += c_J*(1.0-w0);
              else
                lambda_star[0][0][zi] += c_J*(0.5*dt[zi]);
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

                r0[zi]=rho[(xi-xl)*nz+(zi-zl)];
                k0[zi]=kap[(xi-xl)*nz+(zi-zl)];
              }
              for(int zi=zmin;zi<=zmax;++zi){
                dt[zi]=ds3*(k_upw[zi]*r_upw[zi]+k0[zi]*r0[zi])+ds6*(k0[zi]*r_upw[zi]+k_upw[zi]*r0[zi]);
                expo[zi]=exp(-dt[zi]);
                double w0=(1.0-expo[zi])/dt[zi];

                if (dt[zi] > dtau_min)
                  lambda_star[0][xi][zi] += c_J*(1.0-w0);
                else
                  lambda_star[0][xi][zi] += c_J*(0.5*dt[zi]);
              }
            }
          }
        }
      }
    }
  }
}
