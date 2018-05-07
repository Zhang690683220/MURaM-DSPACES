#ifndef WEAK_ENCAPSULATION
#define WEAK_ENCAPSULATION // For performance
#endif

#include <euler3d.H>

typedef PRECISION real; // PRECISION is defined outside

real EULER::HLLELFlux(cState<real>& F,
		      const pState<real>& WL,
		      const pState<real>& WR,
		      const Vector3D<real>& normal) {
  static real a,s,c,sl,sr,u,ul,ur,el,er;
  static cState<real> FL,FR;

  sl = WL.p/WL.d;   sr = WR.p/WR.d;
  ul = WL.V*normal; ur = WR.V*normal;
  el = 0.5*WL.V.sqr()+INV_GAMMA_1*sl;
  er = 0.5*WR.V.sqr()+INV_GAMMA_1*sr;

  s = 1./(1.+sqrt(WR.d/WL.d)); u = 1.-s;
  c = sqrt(GAMMA_1*((el+sl)*s+(er+sr)*u-0.5*(WL.V*s+WR.V*u).sqr()));
  u  = ul*s+ur*u;

  s = u-c; sl = ul-sqrt(GAMMA*sl); sl = sl<s ? sl:s;
  s = u+c; sr = ur+sqrt(GAMMA*sr); sr = sr>s ? sr:s;
  s = sr+sl>0. ? sr : -sl;

  if( sl>0. ) {
    F.d = WL.d*ul;
    F.M = WL.V*F.d; F.M += normal*WL.p;
    F.e = F.d*el+ul*WL.p;
  }
  else if( sr<0. ) {
    F.d = WR.d*ur;
    F.M = WR.V*F.d; F.M += normal*WR.p;
    F.e = F.d*er+ur*WR.p;
  }
  else {
    FL.d = WL.d; FL.M = WL.V*WL.d; FL.e = WL.d*el;
    FR.d = WR.d; FR.M = WR.V*WR.d; FR.e = WR.d*er;
    F = FR-FL;

    FL *= ul; FL.M += normal*WL.p; FL.e += ul*WL.p;
    FR *= ur; FR.M += normal*WR.p; FR.e += ur*WR.p;

    a = c*(c*fabs(F.d)+fabs(F.M.x)+fabs(F.M.y)+fabs(F.M.z))+fabs(F.e);
    a = a==0. ? 0.:
      1.-(c*(c*fabs(FR.d-FL.d-u*F.d)+
	     fabs(FR.M.x-FL.M.x-u*F.M.x)+
	     fabs(FR.M.y-FL.M.y-u*F.M.y)+
	     fabs(FR.M.z-FL.M.z-u*F.M.z))+
	  fabs(FR.e-FL.e-u*F.e))/(c*a);
    a = a>0. ? a : 0.;

    if( u>0. ) {
      c = sl/(sr-sl);
      F *= c*((1.-a)*sr+a*u); F += FL*(1.+c); F -= FR*c;
    }
    else {
      c = sr/(sr-sl);
      F *= c*((1.-a)*sl+a*u); F += FL*c; F += FR*(1.-c);
    }
  }

  return s;
}

real EULER::HLLELFlux(cState<real>& F,
		      const cState<real>& UL,
		      const cState<real>& UR,
		      const Vector3D<real>& normal) {
  static real a,c,s,sl,sr,u,ul,ur,pl,pr;
  static cState<real> FL,FR;

  ul = UL.M*normal;   ur = UR.M*normal;
  pl = UL.Pressure(); pr = UR.Pressure();

  s = 1./(1.+sqrt(UR.d/UL.d)); u = (1.-s)/UR.d; s /= UL.d;
  c = sqrt(GAMMA_1*((UL.e+pl)*s+(UR.e+pr)*u-0.5*(UL.M*s+UR.M*u).sqr()));
  u = ul*s+ur*u; ul /= UL.d; ur /= UR.d;

  s = u-c; sl = ul-sqrt(GAMMA*pl/UL.d); sl = sl<s ? sl:s;
  s = u+c; sr = ur+sqrt(GAMMA*pr/UR.d); sr = sr>s ? sr:s;
  s = sr+sl>0. ? sr : -sl;

  if( sl>0.) {
    F = UL*ul; F.M += normal*pl; F.e += pl*ul;
  }
  else if( sr<0. ) {
    F = UR*ur; F.M += normal*pr; F.e += pr*ur;
  }
  else {
    F  = UR-UL;
    FL = UL*ul; FL.M += normal*pl; FL.e += pl*ul;
    FR = UR*ur; FR.M += normal*pr; FR.e += pr*ur;

    a = c*(c*fabs(F.d)+fabs(F.M.x)+fabs(F.M.y)+fabs(F.M.z))+fabs(F.e);

    a = a==0. ? 0.:
      1.-(c*(c*fabs(FR.d-FL.d-u*F.d)+
	     fabs(FR.M.x-FL.M.x-u*F.M.x)+
	     fabs(FR.M.y-FL.M.y-u*F.M.y)+
	     fabs(FR.M.z-FL.M.z-u*F.M.z))+
	  fabs(FR.e-FL.e-u*F.e))/(c*a);
    a = a>0. ? a : 0.;

    if( u>0. ) {
      c = sl/(sr-sl);
      F *= c*((1.-a)*sr+a*u); F += FL*(1.+c); F -= FR*c;
    }
    else {
      c = sr/(sr-sl);
      F *= c*((1.-a)*sl+a*u); F += FL*c; F += FR*(1.-c);
    }
  }

  return s;
}
