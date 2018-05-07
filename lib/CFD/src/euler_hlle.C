#ifndef WEAK_ENCAPSULATION
#define WEAK_ENCAPSULATION // For performance
#endif

#include <euler3d.H>

typedef PRECISION real; // PRECISION is defined outside

real EULER::HLLEFlux(cState<real>& F,
		     const pState<real>& WL,
		     const pState<real>& WR,
		     const Vector3D<real>& normal) {
  static real s,c,sl,sr,u,ul,ur,el,er;

  sl = WL.p/WL.d;   sr = WR.p/WR.d;
  ul = WL.V*normal; ur = WR.V*normal;
  el = 0.5*WL.V.sqr()+INV_GAMMA_1*sl;
  er = 0.5*WR.V.sqr()+INV_GAMMA_1*sr;

  s = 1./(1.+sqrt(WR.d/WL.d)); u = 1.-s;
  c = sqrt(GAMMA_1*((el+sl)*s+(er+sr)*u-0.5*(WL.V*s+WR.V*u).sqr()));
  u  = ul*s+ur*u;

  s = u-c; sl = ul-sqrt(GAMMA*sl); sl = sl<s ? sl:s; sl = sl<0. ? sl:0.;
  s = u+c; sr = ur+sqrt(GAMMA*sr); sr = sr>s ? sr:s; sr = sr>0. ? sr:0.;
  s = sr+sl>0. ? sr : -sl;

  c = sr-sl; sr /= c; sl /= c;
  u = sr*(ul-sl*c)*WL.d;
  c = sl*(ur-sr*c)*WR.d;

  F.d = u-c;
  F.M = WL.V*u; F.M -= WR.V*c; F.M += normal*(sr*WL.p-sl*WR.p);
  F.e = u*el-c*er+sr*ul*WL.p-sl*ur*WR.p;

  return s;
}

real EULER::HLLEFlux(cState<real>& F,
		     const cState<real>& UL,
		     const cState<real>& UR,
		     const Vector3D<real>& normal) {
  static real c,s,sl,sr,u,ul,ur,pl,pr;

  ul = UL.M*normal;   ur = UR.M*normal;
  pl = UL.Pressure(); pr = UR.Pressure();

  s = 1./(1.+sqrt(UR.d/UL.d)); u = (1.-s)/UR.d; s /= UL.d;
  c = sqrt(GAMMA_1*((UL.e+pl)*s+(UR.e+pr)*u-0.5*(UL.M*s+UR.M*u).sqr()));
  u = ul*s+ur*u; ul /= UL.d; ur /= UR.d;

  s = u-c; sl = ul-sqrt(GAMMA*pl/UL.d); sl = sl<s ? sl:s; sl = sl<0. ? sl:0.;
  s = u+c; sr = ur+sqrt(GAMMA*pr/UR.d); sr = sr>s ? sr:s; sr = sr>0. ? sr:0.;
  s = sr+sl>0. ? sr : -sl;

  c = sr-sl; sr /= c; sl /= c;
  u = sr*(ul-sl*c);
  c = sl*(ur-sr*c);

  F.d = u*UL.d-c*UR.d;
  F.M = UL.M*u; F.M -= UR.M*c; F.M += normal*(sr*pl-sl*pr);
  F.e = u*UL.e-c*UR.e+sr*ul*pl-sl*ur*pr;

  return s;
}
