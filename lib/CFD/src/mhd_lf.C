#ifndef WEAK_ENCAPSULATION
#define WEAK_ENCAPSULATION // For performance
#endif

#include <mhd3d.H>

typedef PRECISION real; // PRECISION is defined outside

real MHD::LaxFriedrichsFlux(cState<real>& F,
			    const pState<real>& WL,
			    const pState<real>& WR,
			    const Vector3D<real>& normal) {
  static real s,sl,sr,ul,ur,bl,br,pl,pr;

  ul = WL.V*normal; ur = WR.V*normal;
  bl = WL.B*normal; br = WR.B*normal;
  pl = WL.p+0.5*WL.B.sqr(); pr = WR.p+0.5*WR.B.sqr();

  sl = pl+0.5*GAMMA_2*WL.p; sr = pr+0.5*GAMMA_2*WR.p;
  sl = fabs(ul)+sqrt((sl+sqrt(fabs(sl*sl-GAMMA*WL.p*bl*bl)))/WL.d);
  sr = fabs(ur)+sqrt((sr+sqrt(fabs(sr*sr-GAMMA*WR.p*br*br)))/WR.d);
  s  = sl > sr ? sl : sr;

  sl = 0.5*WL.d*(ul+s); sr = 0.5*WR.d*(ur-s);

  F.d  = sl+sr;
  F.M  = WL.V*sl;
  F.M += WR.V*sr;
  F.M += normal*(0.5*(pl+pr));
  F.M -= WL.B*(0.5*bl);
  F.M -= WR.B*(0.5*br);
  F.e  = sl*(0.5*WL.V.sqr()+(pl-GAMMA_2*INV_GAMMA_1*WL.p)/WL.d);
  F.e += sr*(0.5*WR.V.sqr()+(pr-GAMMA_2*INV_GAMMA_1*WR.p)/WR.d);
  F.e += 0.5*(ul*pl+ur*pr-bl*(WL.B*WL.V)-br*(WR.B*WR.V));
  F.B  = WL.B*(sl/WL.d);
  F.B += WR.B*(sr/WR.d);
  F.B -= WL.V*(0.5*bl);
  F.B -= WR.V*(0.5*br);

  return s;
}

real MHD::LaxFriedrichsFlux(cState<real>& F,
			    const cState<real>& UL,
			    const cState<real>& UR,
			    const Vector3D<real>& normal) {
  static real s,sl,sr,ul,ur,bl,br,pl,pr;

  bl = UL.B*normal; br = UR.B*normal;
  ul = UL.Pressure(); ur = UR.Pressure();
  pl = ul+0.5*UL.B.sqr(); pr = ur+0.5*UR.B.sqr();

  sl = pl+0.5*GAMMA_2*ul; sr = pr+0.5*GAMMA_2*ur;
  sl = sqrt((sl+sqrt(fabs(sl*sl-GAMMA*ul*bl*bl)))/UL.d);
  sr = sqrt((sr+sqrt(fabs(sr*sr-GAMMA*ur*br*br)))/UR.d);

  ul = (UL.M*normal)/UL.d; ur = (UR.M*normal)/UR.d;
  sl += fabs(ul); sr += fabs(ur);
  s  = sl > sr ? sl : sr;

  sl = 0.5*(ul+s); sr = 0.5*(ur-s);

  F.d  = sl*UL.d+sr*UR.d;
  F.M  = UL.M*sl;
  F.M += UR.M*sr;
  F.M += normal*(0.5*(pl+pr));
  F.M -= UL.B*(0.5*bl);
  F.M -= UR.B*(0.5*br);
  F.e  = sl*UL.e+sr*UR.e+0.5*(ul*pl+ur*pr-bl*(UL.B*UL.M)/UL.d-br*(UR.B*UR.M)/UR.d);
  F.B  = UL.B*sl;
  F.B += UR.B*sr;
  F.B -= UL.M*(0.5*bl/UL.d);
  F.B -= UR.M*(0.5*br/UR.d);

  return s;
}
