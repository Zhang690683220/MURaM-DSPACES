#ifndef WEAK_ENCAPSULATION
#define WEAK_ENCAPSULATION // For performance
#endif

#include <euler3d.H>

typedef PRECISION real; // PRECISION is defined outside

inline real min(const real& a, const real& b) { return a<b ? a : b; }
inline real max(const real& a, const real& b) { return a>b ? a : b; }

inline real InnerProduct(const EULER::cState<real>& U,
			 const EULER::cState<real>& V,
			 const real& c) {
  return c*c*(c*c*U.d*V.d+U.M*V.M)+U.e*V.e;
}

/*
real EULER::TestFlux(cState<real>& F,
		     const cState<real>& UL,
		     const cState<real>& UR,
		     const Vector3D<real>& normal) {
  static real lambda,a,s,cl,cr,du,df;
  static cState<real> DU,DF;

  cl = UL.SoundSpeed();
  cr = UR.SoundSpeed();

  s = 0.5*(cl+cr);

  lambda = max(fabs(normal*UL.Velocity())+cl,
	       fabs(normal*UR.Velocity())+cr);

  DU = UR-UL;
  du = InnerProduct(DU,DU,s);

  DF = UR.Flux(normal)-UL.Flux(normal);
  df = InnerProduct(DF,DF,s);

  if( du==0. || df==0. )
    UL.SetFlux(F,normal);
  else {
    s = fabs(InnerProduct(DU,DF,s))/du;
    a = s*sqrt(du/df);

    F = 0.5*(UL.Flux(normal)+UR.Flux(normal)-((1.-a)*lambda+a*s)*DU);
  }

  return lambda;
}
*/

real EULER::TestFlux(cState<real>& F,
		     const cState<real>& UL,
		     const cState<real>& UR,
		     const Vector3D<real>& normal) {
  static real a,s,sp,sm,cl,cr,ul,ur,du2,df2;
  static cState<real> DU,DF;

  ul = normal*UL.Velocity();
  ur = normal*UR.Velocity();

  cl = UL.SoundSpeed();
  cr = UR.SoundSpeed();

  sm = min(ul-cl,ur-cr);
  sp = max(ul+cl,ur+cr);

  s = max(sp,-sm);

  DU  = UR-UL;
  du2 = InnerProduct(DU,DU,s);

  DF  = UR.Flux(normal)-UL.Flux(normal);
  df2 = InnerProduct(DF,DF,s);

  if( du2==0. || df2==0. )
    UL.SetFlux(F,normal);
  else {
    s = InnerProduct(DU,DF,s)/du2;
    a = fabs(s)*sqrt(du2/df2);

    sm = min(sm,s); sm = min(sm,0.);
    sp = max(sp,s); sp = max(sp,0.);

    if( s>0. )
      F = UL.Flux(normal)+(sm/(sp-sm))*(((1.-a)*sp+a*s)*DU-DF);
    else
      F = UR.Flux(normal)+(sp/(sp-sm))*(((1.-a)*sm+a*s)*DU-DF);
  }
  return max(sp,-sm);
}


real EULER::TestFlux(cState<real>& F,
		     const pState<real>& UL,
		     const pState<real>& UR,
		     const Vector3D<real>& normal) {
  return TestFlux(F,UL.Cons(),UR.Cons(),normal);
}
