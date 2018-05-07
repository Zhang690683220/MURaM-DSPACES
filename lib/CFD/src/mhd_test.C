#ifndef WEAK_ENCAPSULATION
#define WEAK_ENCAPSULATION // For performance
#endif

#include <mhd3d.H>

typedef PRECISION real; // PRECISION is defined outside

inline real min(const real& a, const real& b) { return a<b ? a : b; }
inline real max(const real& a, const real& b) { return a>b ? a : b; }

inline real H(const real& a) { return a < 0. ? 0. : 1.; }

inline real InnerProduct(const MHD::cState<real>& U,
			 const MHD::cState<real>& V) {
  return U.d*V.d+U.M*V.M+U.e*V.e+U.B*V.B;
}

real MHD::TestFlux(cState<real>& F,
		   const cState<real>& UL,
		   const cState<real>& UR,
		   const Vector3D<real>& normal) {
  static real a,s,sp,sm,cl,cr,ul,ur,du2,df2;
  static cState<real> DU,DF;

  ul = normal*UL.Velocity();
  ur = normal*UR.Velocity();

  cl = UL.FastSpeed(normal);
  cr = UR.FastSpeed(normal);

  sm = min(ul-cl,ur-cr);
  sp = max(ul+cl,ur+cr);

  DU  = UR-UL;
  du2 = InnerProduct(DU,DU);

  if( du2==0. )
    UL.SetFlux(F,normal);
  else {

    DF  = UR.Flux(normal)-UL.Flux(normal);
    df2 = InnerProduct(DF,DF);

    if( df2==0. )
      if( ur*UR.Entropy() < ul*UL.Entropy() )
	UL.SetFlux(F,normal);
      else
	F = (UL.Flux(normal)*sp-UR.Flux(normal)*sm+DU*(sm*sp))*(1./(sp-sm));
    else {
      s = InnerProduct(DU,DF)/du2;
      a = (s*s)*(du2/df2)*H((s-ur)*UR.Entropy()-(s-ul)*UL.Entropy());

      sm = min(sm,s); sm = min(sm,0.);
      sp = max(sp,s); sp = max(sp,0.);

      if( s>0. )
	F = UL.Flux(normal)+(sm/(sp-sm))*(((1.-a)*sp+a*s)*DU-DF);
      else
	F = UR.Flux(normal)+(sp/(sp-sm))*(((1.-a)*sm+a*s)*DU-DF);
    }
  }
  return max(sp,-sm);
}

real MHD::TestFlux(cState<real>& F,
		   const pState<real>& UL,
		   const pState<real>& UR,
		   const Vector3D<real>& normal) {
  return TestFlux(F,UL.Cons(),UR.Cons(),normal);
}
