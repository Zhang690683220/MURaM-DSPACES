#ifndef WEAK_ENCAPSULATION
#define WEAK_ENCAPSULATION // For performance
#endif

#include <euler3d.H>

typedef PRECISION real; // PRECISION is defined outside

const real TOLERANCE = sizeof(real)==sizeof(float) ? 1.e-5 : 1.e-12;
const real GAMMA_PLUS_1 = EULER::GAMMA+1.;

inline real Dstar(const real& d, const real& p, const real& pstar) {
  return (pstar < p) ? d*pow(pstar/p,EULER::INV_GAMMA) :
    d*(EULER::GAMMA_1*p+GAMMA_PLUS_1*pstar)/
    (GAMMA_PLUS_1*p+EULER::GAMMA_1*pstar);
}

inline real MassFlux(const real& d, const real& p, const real& pstar) {
  return ( pstar > p*(1.-TOLERANCE) ) ?
    d*sqrt(EULER::GAMMA*(p/d)*(1.+0.5*EULER::INV_GAMMA*GAMMA_PLUS_1*(pstar/p-1.))) :
    d*sqrt(EULER::GAMMA*p/d)*0.5*EULER::INV_GAMMA*EULER::GAMMA_1*(1.-pstar/p)/
    (1.-pow(pstar/p,0.5*EULER::INV_GAMMA*EULER::GAMMA_1));
}

inline real PstarOffset(const real& dl, const real& ul, const real& pl,
			const real& dr, const real& ur, const real& pr,
			const real& pstar) {
  real ml = MassFlux(dl,pl,pstar);
  real mr = MassFlux(dr,pr,pstar);
  return 1.-(ml*pr+mr*pl-ml*mr*(ur-ul))/(pstar*(ml+mr));
}

real Pstar(const real& dl, const real& ul, real cl,
	   const real& dr, const real& ur, real cr) {
  using namespace EULER;

  static real pl,pr,pmin,pmax,pstar,Fmin,Fmax,Fstar;

  pl = dl*cl*cl*INV_GAMMA;
  pr = dr*cr*cr*INV_GAMMA;

  pmax = pow((cl+cr-0.5*GAMMA_1*(ur-ul))/
	     (cl/pow(pl,0.5*INV_GAMMA*GAMMA_1)+
	      cr/pow(pr,0.5*INV_GAMMA*GAMMA_1)),
	     2.*GAMMA_INV_GAMMA_1);

  Fmax  = PstarOffset(dl,ul,pl,dr,ur,pr,pmax);

  if( (pmax<pl && pmax<pr) || fabs(Fmax)<TOLERANCE )
    return pmax;

  cl *= dl; // cl is no longer the speed of sound
  cr *= dr; // cr -------#-------#-------#-------

  pmin = (cl*pr+cr*pl-cl*cr*(ur-ul))/(cl+cr);
  Fmin  = PstarOffset(dl,ul,pl,dr,ur,pr,pmin);

  pstar = 0.5*(pmin+pmax);
  Fstar = PstarOffset(dl,ul,pl,dr,ur,pr,pstar);

  while( (pmax-pmin) > TOLERANCE*pstar && fabs(Fstar) > TOLERANCE ) {

#ifndef NDEBUG
    if( Fmin*Fmax>0. ) {
      cerr << "GODUNOV SOLVER: CONVERGENCE PROBLEM" << endl;
      return -1.;
    }
#endif

    if( Fstar*Fmin>0. ) {
      pmin = pstar;
      Fmin = Fstar;
    }
    else {
      pmax = pstar;
      Fmax = Fstar;
    }

    pstar = (Fmax*pmin-Fmin*pmax)/(Fmax-Fmin);

#ifndef NDEBUG
    if( pstar<0. ) {
      cerr << "GODUNOV SOLVER: NEGATIVE PRESSURE" << endl;
      return pstar;
    }
#endif

    Fstar = PstarOffset(dl,ul,pl,dr,ur,pr,pstar);
  }
  return pstar;
}


real EULER::GodunovFlux(EULER::cState<real>& F,
			const EULER::pState<real>& WL,
			const EULER::pState<real>& WR,
			const Vector3D<real>& normal) {
  static pState<real> Wstar;
  static real unl,unr,cl,cr,cstar,dstar,pstar,ustar;

  unl = WL.V*normal;     unr = WR.V*normal;
  cl  = WL.SoundSpeed(); cr  = WR.SoundSpeed();

#ifndef NDEBUG
  if( GAMMA_1*unl+2.*cl < GAMMA_1*unr-2.*cr ) {
    cerr << "GODUNOV SOLVER: VACUUM DETECTED" << endl;

    real ustar;
    real unlstar = unl+2.*INV_GAMMA_1*cl;
    real unrstar = unr+2.*INV_GAMMA_1*cr;

    if( unlstar>0. )
      if( unl-cl>0. )
	WL.SetFlux(F,normal);
      else {
	ustar = unlstar*GAMMA_1/(GAMMA+1.);
	Wstar.d = WL.d*pow(ustar/cl,2.*INV_GAMMA_1);
	Wstar.V = WL.V+normal*(ustar-unl);
	Wstar.p = WL.p*pow(ustar/cl,2.*GAMMA_INV_GAMMA_1);
	Wstar.SetFlux(F,normal);
      }
    else if( unrstar<0. )
      if( unr+cr<0. )
	WR.SetFlux(F,normal);
      else {
	ustar = unrstar*GAMMA_1/(GAMMA+1.);
	Wstar.d = WR.d*pow(-ustar/cr,2.*INV_GAMMA_1);
	Wstar.V = WR.V+normal*(ustar-unr);
	Wstar.p = WR.p*pow(-ustar/cr,2.*GAMMA_INV_GAMMA_1);
	Wstar.SetFlux(F,normal);
      }
    else
      F.zero();
  }
  else
#endif
    {
      pstar = Pstar(WL.d,unl,cl,WR.d,unr,cr);
      ustar = unr+(pstar-WR.p)/MassFlux(WR.d,WR.p,pstar);

      if( ustar>0. )
	if( pstar>WL.p )
	  if( unl-MassFlux(WL.d,WL.p,pstar)/WL.d>0. )
	    WL.SetFlux(F,normal);
	  else {
	    Wstar.d = Dstar(WL.d,WL.p,pstar);
	    Wstar.V = WL.V+normal*(ustar-unl);
	    Wstar.p = pstar;
	  }
	else
	  if( unl-cl>0. )
	    WL.SetFlux(F,normal);
	  else {
	    dstar = Dstar(WL.d,WL.p,pstar);
	    cstar = sqrt(GAMMA*pstar/dstar);
	    if( ustar-cstar>0. ) {
	      ustar = (ustar*cl-unl*cstar)/(ustar-cstar-unl+cl);
	      dstar = WL.d*pow(ustar/cl,2.*INV_GAMMA_1);
	      pstar = WL.p*pow(dstar/WL.d,GAMMA);
	    }
	    Wstar.d = dstar;
	    Wstar.V = WL.V+normal*(ustar-unl);
	    Wstar.p = pstar;
	  }
      else
	if( pstar>WR.p )
	  if( unr+MassFlux(WR.d,WR.p,pstar)/WR.d<0. )
	    WR.SetFlux(F,normal);
	  else {
	    Wstar.d = Dstar(WR.d,WR.p,pstar);
	    Wstar.V = WR.V+normal*(ustar-unr);
	    Wstar.p = pstar;
	  }
	else
	  if( unr+cr<0. )
	    WR.SetFlux(F,normal);
	  else {
	    dstar = Dstar(WR.d,WR.p,pstar);
	    cstar = sqrt(GAMMA*pstar/dstar);
	    if( ustar+cstar<0. ) {
	      ustar = (ustar*cr-unr*cstar)/(unr+cr-ustar-cstar);
	      dstar = WR.d*pow(-ustar/cr,2.*INV_GAMMA_1);
	      pstar = WR.p*pow(dstar/WR.d,GAMMA);
	    }
	    Wstar.d = dstar;
	    Wstar.V = WR.V+normal*(ustar-unr);
	    Wstar.p = pstar;
	  }

      Wstar.SetFlux(F,normal);
    }

  cl += fabs(unl); cr += fabs(unr);
  return cl > cr ? cl : cr;
}

real EULER::GodunovFlux(cState<real>& F,
			const cState<real>& UL,
			const cState<real>& UR,
			const Vector3D<real>& normal) {
  return GodunovFlux(F,UL.Prim(),UR.Prim(),normal);
}
