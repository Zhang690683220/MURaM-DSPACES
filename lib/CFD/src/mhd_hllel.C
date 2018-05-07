#ifndef WEAK_ENCAPSULATION
#define WEAK_ENCAPSULATION // For performance
#endif

#include <mhd3d.H>

typedef PRECISION real; // PRECISION is defined outside

real MHD::HLLELFlux(cState<real>& F,
		    const pState<real>& WL,
		    const pState<real>& WR,
		    const Vector3D<real>& normal) {
  static real a,d,c,s,sl,sr,u,ul,ur,b,bl,br,pl,pr,el,er;
  static cState<real> FL,FR;

  ul = WL.V*normal; ur = WR.V*normal;
  bl = WL.B*normal; br = WR.B*normal;
  pl = WL.p+0.5*WL.B.sqr(); pr = WR.p+0.5*WR.B.sqr();
  el = 0.5*WL.V.sqr()+(pl-GAMMA_2*INV_GAMMA_1*WL.p)/WL.d;
  er = 0.5*WR.V.sqr()+(pr-GAMMA_2*INV_GAMMA_1*WR.p)/WR.d;

  s   = sqrt(WR.d/WL.d); d = s*WL.d;
  s   = 1./(1.+s); u = 1.-s;
  F.B = WR.B*s+WL.B*u; b = F.B.sqr();
  c   = GAMMA_1*((el+pl/WL.d)*s+(er+pr/WR.d)*u-0.5*(WL.V*s+WR.V*u).sqr()-b/d);
  sl  = (1.-2*s)*(F.B*(WR.B-WL.B)); sl = sl>0. ? sl:0.; sl = 0.5*(c+(b-0.5*sl)/d);
  b   = br*s+bl*u; c = sqrt(sl+sqrt(fabs(sl*sl-c*b*b/d))); u = ul*s+ur*u;

  sl = pl+0.5*GAMMA_2*WL.p; sl = ul-sqrt((sl+sqrt(fabs(sl*sl-GAMMA*WL.p*bl*bl)))/WL.d);
  sr = pr+0.5*GAMMA_2*WR.p; sr = ur+sqrt((sr+sqrt(fabs(sr*sr-GAMMA*WR.p*br*br)))/WR.d);

  s = u-c; sl = sl<s ? sl:s;
  s = u+c; sr = sr>s ? sr:s;

  if( sl>0. ) {
    F.d = WL.d*ul;
    F.M = WL.V*F.d; F.M += normal*pl; F.M -= WL.B*bl;
    F.e = ul*(WL.d*el+pl)-bl*(WL.V*WL.B);
    F.B = WL.B*ul; F.B -= WL.V*bl;
  }
  else if( sr<0. ) {
    F.d = WR.d*ur;
    F.M = WR.V*F.d; F.M += normal*pr; F.M -= WR.B*br;
    F.e = ur*(WR.d*er+pr)-br*(WR.V*WR.B);
    F.B = WR.B*ur; F.B -= WR.V*br;
  }
  else {
    FL.d = WL.d; FL.M = WL.V*WL.d; FL.e = WL.d*el; FL.B = WL.B;
    FR.d = WR.d; FR.M = WR.V*WR.d; FR.e = WR.d*er; FR.B = WR.B;
    F = FR-FL;

    FL *= ul; FL.M += normal*pl; FL.M -= WL.B*bl; FL.e += ul*pl-bl*(WL.V*WL.B); FL.B -= WL.V*bl;
    FR *= ur; FR.M += normal*pr; FR.M -= WR.B*br; FR.e += ur*pr-br*(WR.V*WR.B); FR.B -= WR.V*br;

    s = sqrt(d);
    a = c*(c*fabs(F.d)+fabs(F.M.x)+fabs(F.M.y)+fabs(F.M.z)+
	   s*(fabs(F.B.x)+fabs(F.B.y)+fabs(F.B.z)))+fabs(F.e);
    a = a==0. ? 0. :
      1.-(c*(c*fabs(FR.d-FL.d-u*F.d)+
	     fabs(FR.M.x-FL.M.x-u*F.M.x)+
	     fabs(FR.M.y-FL.M.y-u*F.M.y)+
	     fabs(FR.M.z-FL.M.z-u*F.M.z)+
	     s*(fabs(FR.B.x-FL.B.x-u*F.B.x)+
		fabs(FR.B.y-FL.B.y-u*F.B.y)+
		fabs(FR.B.z-FL.B.z-u*F.B.z)))
	  +fabs(FR.e-FL.e-u*F.e))/(c*a);
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

  s = sr+sl>0. ? sr : -sl;

  u = 0.5*fabs(s)*(br-bl);
  F.e -= u*b; F.B -= normal*u;

  return s;
}

real MHD::HLLELFlux(cState<real>& F,
		    const cState<real>& UL,
		    const cState<real>& UR,
		    const Vector3D<real>& normal) {
  static real a,d,c,s,sl,sr,u,ul,ur,b,bl,br,pl,pr;
  static cState<real> FL,FR;

  bl = UL.B*normal; br = UR.B*normal;
  pl = 0.5*UL.B.sqr(); pr = 0.5*UR.B.sqr();
  ul = GAMMA_1*(UL.e-0.5*UL.M.sqr()/UL.d-pl);
  ur = GAMMA_1*(UR.e-0.5*UR.M.sqr()/UR.d-pr);
  pl += ul; pr += ur;

  s   = sqrt(UR.d/UL.d); d = s*UL.d;
  s   = 1./(1.+s); u = 1.-s; sl = s/UL.d; sr = u/UR.d;
  F.B = UR.B*s+UL.B*u; b = F.B.sqr();
  c   = GAMMA_1*((UL.e+pl)*sl+(UR.e+pr)*sr-0.5*(UL.M*sl+UR.M*sr).sqr()-b/d);
  sl  = (1.-2.*s)*(F.B*(UR.B-UL.B)); sl = sl>0. ? sl:0.; sl = 0.5*(c+(b-0.5*sl)/d);
  b   = br*s+bl*u; c  = sqrt(sl+sqrt(fabs(sl*sl-c*b*b/d)));

  sl = pl+0.5*GAMMA_2*ul; sl = sqrt((sl+sqrt(fabs(sl*sl-GAMMA*ul*bl*bl)))/UL.d);
  sr = pr+0.5*GAMMA_2*ur; sr = sqrt((sr+sqrt(fabs(sr*sr-GAMMA*ur*br*br)))/UR.d);
  ul = (UL.M*normal)/UL.d; ur = (UR.M*normal)/UR.d; u = ul*s+ur*u;

  s = u-c; sl = sl<s ? sl:s;
  s = u+c; sr = sr>s ? sr:s;

  if( sl>0. ) {
    F = UL*ul; s = bl/UL.d;
    F.M += normal*pl; F.M -= UL.B*bl;
    F.e += ul*pl-s*(UL.M*UL.B);
    F.B -= UL.M*s;
  }
  else if( sr<0. ) {
    F = UR*ur; s = br/UR.d;
    F.M += normal*pr; F.M -= UR.B*br;
    F.e += ur*pr-s*(UR.M*UR.B);
    F.B -= UR.M*s;
  }
  else {
    F = UR-UL;

    FL = UL*ul; s = bl/UL.d;
    FL.M += normal*pl; FL.M -= UL.B*bl;
    FL.e += ul*pl-s*(UL.M*UL.B);
    FL.B -= UL.M*s;

    FR = UR*ur; s = br/UR.d;
    FR.M += normal*pr; FR.M -= UR.B*br;
    FR.e += ur*pr-s*(UR.M*UR.B);
    FR.B -= UR.M*s;

    s = sqrt(d);
    a = c*(c*fabs(F.d)+fabs(F.M.x)+fabs(F.M.y)+fabs(F.M.z)+
	   s*(fabs(F.B.x)+fabs(F.B.y)+fabs(F.B.z)))+fabs(F.e);
    a = a==0. ? 0. :
      1.-(c*(c*fabs(FR.d-FL.d-u*F.d)+
	     fabs(FR.M.x-FL.M.x-u*F.M.x)+
	     fabs(FR.M.y-FL.M.y-u*F.M.y)+
	     fabs(FR.M.z-FL.M.z-u*F.M.z)+
	     s*(fabs(FR.B.x-FL.B.x-u*F.B.x)+
		fabs(FR.B.y-FL.B.y-u*F.B.y)+
		fabs(FR.B.z-FL.B.z-u*F.B.z)))
	  +fabs(FR.e-FL.e-u*F.e))/(c*a);
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

  s = sr+sl>0. ? sr : -sl;

  u = 0.5*fabs(s)*(br-bl);
  F.e -= u*b; F.B -= normal*u;
  
  return s;
}
