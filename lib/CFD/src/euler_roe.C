#ifndef WEAK_ENCAPSULATION
#define WEAK_ENCAPSULATION // For performance
#endif

#include <euler3d.H>

typedef PRECISION real; // PRECISION is defined outside

inline real EntropyFixPlus(real al, real ar, real a) {
  ar = 4.*(ar-al); ar = ar > 0. ? ar : 0.;
  return fabs(a)>0.5*ar ? (a > 0. ? a : 0.) : 0.5*(a*a/ar+a+0.25*ar);
}

inline real EntropyFixMinus(real al, real ar, real a) {
  ar = 4.*(ar-al); ar = ar > 0. ? ar : 0.;
  return fabs(a)>0.5*ar ? (a < 0. ? a : 0.) : -0.5*(a*a/ar-a+0.25*ar);
}

real EULER::RoeFlux(cState<real>& F,
		    const pState<real>& WL,
		    const pState<real>& WR,
		    const Vector3D<real>& normal) {
  static real d,c,cl,cr,u,ul,ur,h,hl,hr;
  static Vector3D<real> V;

  ul = WL.V*normal;           ur = WR.V*normal;
  cl = WL.SoundSpeed();       cr = WR.SoundSpeed();
  hl = WL.SpecificEnthalpy(); hr = WR.SpecificEnthalpy();

  c = sqrt(WR.d/WL.d); d = c*WL.d;
  c = 1./(1.+c); h = 1.-c;
  V = WL.V*c; V += WR.V*h; u = c*ul+h*ur;
  h = c*hl+h*hr; c = sqrt(GAMMA_1*(h-0.5*V.sqr()));

  if( u > 0. ) {
    d = EntropyFixMinus(ul-cl,ur-cr,u-c)*(WR.p-WL.p-d*c*(ur-ul))/(2.*c*c);
    ul *= WL.d;
    F.d = ul+d;
    F.M = WL.V*ul; F.M += V*d; F.M += normal*(WL.p-c*d);
    F.e = hl*ul+(h-u*c)*d;
  }
  else {
    d = EntropyFixPlus(ul+cl,ur+cr,u+c)*(WR.p-WL.p+d*c*(ur-ul))/(2.*c*c);
    ur *= WR.d;
    F.d = ur-d;
    F.M = WR.V*ur; F.M -= V*d; F.M += normal*(WR.p-c*d);
    F.e = hr*ur-(h+u*c)*d;
  }

  return fabs(u)+c;
}

real EULER::RoeFlux(cState<real>& F,
		    const cState<real>& UL,
		    const cState<real>& UR,
		    const Vector3D<real>& normal) {
  static real d,c,cl,cr,u,ul,ur,h,pl,pr;
  static Vector3D<real> V;

  ul = (UL.M*normal)/UL.d; ur = (UR.M*normal)/UR.d;
  pl = UL.Pressure(); pr = UR.Pressure();
  cl = sqrt(GAMMA*pl/UL.d); cr = sqrt(GAMMA*pr/UR.d);

  c = sqrt(UR.d/UL.d); d = c*UL.d;
  c = 1./(1.+c); h = (1.-c)/d; c /= d;
  V = UR.M*c; V += UL.M*h; u = d*(c*ul+h*ur);
  h = c*(UR.e+pr)+h*(UL.e+pl); c = sqrt(GAMMA_1*(h-0.5*V.sqr()));

  if( u > 0. ) {
    d = EntropyFixMinus(ul-cl,ur-cr,u-c)*(pr-pl-d*c*(ur-ul))/(2.*c*c);
    F.d = ul*UL.d+d;
    F.M = UL.M*ul; F.M += V*d; F.M += normal*(pl-c*d);
    F.e = ul*(UL.e+pl)+(h-u*c)*d;
  }
  else {
    d = EntropyFixPlus(ul+cl,ur+cr,u+c)*(pr-pl+d*c*(ur-ul))/(2.*c*c);
    F.d = ur*UR.d-d;
    F.M = UR.M*ur; F.M -= V*d; F.M += normal*(pr-c*d);
    F.e = ur*(UR.e+pr)-(h+u*c)*d;
  }

  return fabs(u)+c;
}
