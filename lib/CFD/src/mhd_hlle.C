#ifndef WEAK_ENCAPSULATION
#define WEAK_ENCAPSULATION // For performance
#endif

#include <mhd3d.H>

typedef PRECISION real; // PRECISION is defined outside

real MHD::HLLEFlux(cState<real>& F,
		   const pState<real>& WL,
		   const pState<real>& WR,
		   const Vector3D<real>& normal) {
  static real d,c,s,sl,sr,u,ul,ur,b,bl,br,pl,pr,el,er;

  ul = WL.V*normal; ur = WR.V*normal;
  bl = WL.B*normal; br = WR.B*normal;
  pl = WL.p+0.5*WL.B.sqr(); pr = WR.p+0.5*WR.B.sqr();
  el = 0.5*WL.V.sqr()+(pl-GAMMA_2*INV_GAMMA_1*WL.p)/WL.d;
  er = 0.5*WR.V.sqr()+(pr-GAMMA_2*INV_GAMMA_1*WR.p)/WR.d;

  s    = sqrt(WR.d/WL.d); d = s*WL.d;
  s    = 1./(1.+s); u = 1.-s;
  F.B  = WR.B*s+WL.B*u; b = F.B.sqr();
  c    = GAMMA_1*((el+pl/WL.d)*s+(er+pr/WR.d)*u-0.5*(WL.V*s+WR.V*u).sqr()-b/d);
  sl   = (1.-2*s)*(F.B*(WR.B-WL.B)); sl = sl>0. ? sl:0.; sl = 0.5*(c+(b-0.5*sl)/d);
  b    = br*s+bl*u; c = sqrt(sl+sqrt(fabs(sl*sl-c*b*b/d))); u = ul*s+ur*u;

  sl = pl+0.5*GAMMA_2*WL.p; sl = ul-sqrt((sl+sqrt(fabs(sl*sl-GAMMA*WL.p*bl*bl)))/WL.d);
  sr = pr+0.5*GAMMA_2*WR.p; sr = ur+sqrt((sr+sqrt(fabs(sr*sr-GAMMA*WR.p*br*br)))/WR.d);

  s = u-c; sl = sl<s ? sl:s; sl = sl<0. ? sl:0.;
  s = u+c; sr = sr>s ? sr:s; sr = sr>0. ? sr:0.;
  s = sr+sl>0. ? sr : -sl;

  c = sr-sl; sr /= c; sl /= c;
  d = sr*(ul-sl*c)*WL.d;
  c = sl*(ur-sr*c)*WR.d;
  u = 0.5*fabs(u)*(br-bl);

  F.d  = d-c;
  F.M  = WL.V*d;
  F.M -= WR.V*c;
  F.M += normal*(sr*pl-sl*pr);
  F.M -= WL.B*(sr*bl);
  F.M += WR.B*(sl*br);
  F.e  = d*el-c*er+sr*(ul*pl-bl*(WL.B*WL.V))-sl*(ur*pr-br*(WR.B*WR.V))-u*b;
  F.B  = WL.B*(d/WL.d);
  F.B -= WR.B*(c/WR.d);
  F.B -= WL.V*(sr*bl);
  F.B += WR.V*(sl*br);
  F.B -= normal*u;

  return s;
}

real MHD::HLLEFlux(cState<real>& F,
		   const cState<real>& UL,
		   const cState<real>& UR,
		   const Vector3D<real>& normal) {
  static real d,c,s,sl,sr,u,ul,ur,b,bl,br,pl,pr;

  bl = UL.B*normal; br = UR.B*normal;
  pl = 0.5*UL.B.sqr(); pr = 0.5*UR.B.sqr();
  ul = GAMMA_1*(UL.e-0.5*UL.M.sqr()/UL.d-pl);
  ur = GAMMA_1*(UR.e-0.5*UR.M.sqr()/UR.d-pr);
  pl += ul; pr += ur;

  s    = sqrt(UR.d/UL.d); d = s*UL.d;
  s    = 1./(1.+s); u = 1.-s; sl = s/UL.d; sr = u/UR.d;
  F.B  = UR.B*s+UL.B*u; b = F.B.sqr();
  c    = GAMMA_1*((UL.e+pl)*sl+(UR.e+pr)*sr-0.5*(UL.M*sl+UR.M*sr).sqr()-b/d);
  sl   = (1.-2.*s)*(F.B*(UR.B-UL.B)); sl = sl>0. ? sl:0.; sl = 0.5*(c+(b-0.5*sl)/d);
  b    = br*s+bl*u; c  = sqrt(sl+sqrt(fabs(sl*sl-c*b*b/d)));

  sl = pl+0.5*GAMMA_2*ul; sl = sqrt((sl+sqrt(fabs(sl*sl-GAMMA*ul*bl*bl)))/UL.d);
  sr = pr+0.5*GAMMA_2*ur; sr = sqrt((sr+sqrt(fabs(sr*sr-GAMMA*ur*br*br)))/UR.d);
  ul = (UL.M*normal)/UL.d; ur = (UR.M*normal)/UR.d; u = ul*s+ur*u;

  s = u-c; sl = sl<s ? sl:s; sl = sl<0. ? sl:0.;
  s = u+c; sr = sr>s ? sr:s; sr = sr>0. ? sr:0.;
  s = sr+sl>0. ? sr : -sl;

  c = sr-sl; sr /= c; sl /= c;
  d = sr*(ul-sl*c);
  c = sl*(ur-sr*c);
  u = 0.5*fabs(u)*(br-bl);
  bl /= UL.d; br /= UR.d;

  F.d  = d*UL.d-c*UR.d;
  F.M  = UL.M*d;
  F.M -= UR.M*c;
  F.M += normal*(sr*pl-sl*pr);
  F.M -= UL.B*(sr*bl*UL.d);
  F.M += UR.B*(sl*br*UR.d);
  F.e  = d*UL.e-c*UR.e+sr*(ul*pl-bl*(UL.B*UL.M))-sl*(ur*pr-br*(UR.B*UR.M))-u*b;
  F.B  = UL.B*d;
  F.B -= UR.B*c;
  F.B -= UL.M*(sr*bl);
  F.B += UR.M*(sl*br);
  F.B -= normal*u;

  return s;
}
