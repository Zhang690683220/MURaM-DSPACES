#ifndef WEAK_ENCAPSULATION
#define WEAK_ENCAPSULATION // For performance
#endif

#include <euler3d.H>

typedef PRECISION real; // PRECISION is defined outside

real EULER::LaxFriedrichsFlux(cState<real>& F,
			      const pState<real>& WL,
			      const pState<real>& WR,
			      const Vector3D<real>& normal) {
  static real s,sl,sr,ul,ur;

  ul = WL.V*normal; ur = WR.V*normal;

  sl = fabs(ul)+WL.SoundSpeed();
  sr = fabs(ur)+WR.SoundSpeed();
  s  = sl > sr ? sl : sr;

  sl = 0.5*WL.d*(ul+s);
  sr = 0.5*WR.d*(ur-s);

  F.d  = sl+sr;
  F.M  = WL.V*sl; F.M += WR.V*sr; F.M += normal*(0.5*(WL.p+WR.p));
  F.e  = sl*WL.SpecificEnergy()+sr*WR.SpecificEnergy()+0.5*(ul*WL.p+ur*WR.p);

  return s;
}

real EULER::LaxFriedrichsFlux(cState<real>& F,
			      const cState<real>& UL,
			      const cState<real>& UR,
			      const Vector3D<real>& normal) {
  static real s,sl,sr,ul,ur,pl,pr;

  ul = (UL.M*normal)/UL.d; ur = (UR.M*normal)/UR.d;
  pl = UL.Pressure(); pr = UR.Pressure();

  sl = fabs(ul)+sqrt(GAMMA*pl/UL.d);
  sr = fabs(ur)+sqrt(GAMMA*pr/UR.d);
  s  = sl > sr ? sl : sr;

  sl = 0.5*(ul+s);
  sr = 0.5*(ur-s);

  F.d = sl*UL.d+sr*UR.d;
  F.M = UL.M*sl; F.M += UR.M*sr; F.M += normal*(0.5*(pl+pr));
  F.e = sl*UL.e+sr*UR.e+0.5*(ul*pl+ur*pr);

  return s;
}
