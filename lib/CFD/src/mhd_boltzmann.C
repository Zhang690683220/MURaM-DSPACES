#ifndef WEAK_ENCAPSULATION
#define WEAK_ENCAPSULATION // For performance
#endif

#include <mhd3d.H>

typedef PRECISION real; // PRECISION is defined outside

real MHD::BoltzmannFlux(cState<real>& F,
			const pState<real>& WL,
			const pState<real>& WR,
			const Vector3D<real>& normal) {
  static real a,s,bc,p,Vn,Bn,dBn;

  // Left flux
  p  = WL.p+0.5*WL.B.sqr();
  bc = sqrt(2.*p/WL.d);
  Vn = WL.V*normal;

  a  = -Vn/bc;
  Bn = 0.28209479177387814348*exp(-a*a);
  a  = 0.5*(1.-erf(a));

  bc *= Bn; Vn *= a; Vn += bc; s = Vn>0. ? Vn:0.;
  Bn = WL.B*normal; dBn = -Bn;

  F.d = Vn*WL.d;
  F.M = WL.V*F.d; F.M += normal*(a*p); a *= Bn; F.M -= WL.B*a;
  F.e = Vn*(0.5*WL.d*WL.V.sqr()+2.*p-GAMMA_2*INV_GAMMA_1*WL.p)-a*(WL.B*WL.V)-0.5*bc*(p+Bn*Bn);
  F.B = WL.B*Vn; F.B -= WL.V*a; F.B -= normal*(bc*Bn);

  // Right flux
  p  = WR.p+0.5*WR.B.sqr();
  bc = sqrt(2.*p/WR.d);
  Vn = -(WR.V*normal);

  a  = -Vn/bc;
  Bn = 0.28209479177387814348*exp(-a*a);
  a  = 0.5*(1.-erf(a));

  bc *= Bn; Vn *= a; Vn += bc; s += Vn>0. ? Vn:0.;
  Bn = -(WR.B*normal); dBn += -Bn;

  F.d -= Vn*WR.d;
  F.M -= WR.V*(Vn*WR.d); F.M += normal*(a*p); a *= Bn; F.M += WR.B*a;
  F.e -= Vn*(0.5*WR.d*WR.V.sqr()+2.*p-GAMMA_2*INV_GAMMA_1*WR.p)-a*(WR.B*WR.V)-0.5*bc*(p+Bn*Bn);
  F.B -= WR.B*Vn; F.B += WR.V*a; F.B -= normal*(bc*Bn);

  F.B -= normal*(0.5*s*dBn);

  return s;
}

real MHD::BoltzmannFlux(cState<real>& F,
			const cState<real>& UL,
			const cState<real>& UR,
			const Vector3D<real>& normal) {
  static real a,s,bc,p,Vn,Bn,dBn;

  // Left flux
  p  = GAMMA_1*(UL.e-0.5*UL.M.sqr()/UL.d)-0.5*GAMMA_2*UL.B.sqr();
  bc = sqrt(2.*p/UL.d);
  Vn = (UL.M*normal)/UL.d;

  a  = -Vn/bc;
  Bn = 0.28209479177387814348*exp(-a*a);
  a  = 0.5*(1.-erf(a));

  bc *= Bn; Vn *= a; Vn += bc; s = Vn>0. ? Vn:0.;
  Bn = UL.B*normal; dBn = -Bn;

  F.d = Vn*UL.d;
  F.M = UL.M*Vn; F.M += normal*(a*p); a *= Bn; F.M -= UL.B*a; a /= UL.d;
  F.e = Vn*(UL.e+p)-a*(UL.B*UL.M)-0.5*bc*(p+Bn*Bn);
  F.B = UL.B*Vn; F.B -= UL.M*a; F.B -= normal*(bc*Bn);

  // Right flux
  p  = GAMMA_1*(UR.e-0.5*UR.M.sqr()/UR.d)-0.5*GAMMA_2*UR.B.sqr();
  bc = sqrt(2.*p/UR.d);
  Vn = -(UR.M*normal)/UR.d;

  a  = -Vn/bc;
  Bn = 0.28209479177387814348*exp(-a*a);
  a  = 0.5*(1.-erf(a));

  bc *= Bn; Vn *= a; Vn += bc; s += Vn>0. ? Vn:0.;
  Bn = -(UR.B*normal); dBn += -Bn;

  F.d -= Vn*UR.d;
  F.M -= UR.M*Vn; F.M += normal*(a*p); a *= Bn; F.M += UR.B*a; a /= UR.d;
  F.e -= Vn*(UR.e+p)-a*(UR.B*UR.M)-0.5*bc*(p+Bn*Bn);
  F.B -= UR.B*Vn; F.B += UR.M*a; F.B -= normal*(bc*Bn);

  F.B -= normal*(0.5*s*dBn);

  return s;
}
