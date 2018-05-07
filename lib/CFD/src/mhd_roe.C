#ifndef WEAK_ENCAPSULATION
#define WEAK_ENCAPSULATION // For performance
#endif

#include <mhd3d.H>

typedef PRECISION real; // PRECISION is defined outside

inline real EntropyFixPlus(real al, real ar, real a) {
  ar = 4.*(ar-al); ar = ar > 0. ? ar : 0.;
  return fabs(a)>0.5*ar ? (a > 0. ? a : 0.) : 0.5*(a*a/ar+a+0.25*ar);
}

inline real EntropyFixMinus(real al, real ar, real a) {
  ar = 4.*(ar-al); ar = ar > 0. ? ar : 0.;
  return fabs(a)>0.5*ar ? (a < 0. ? a : 0.) : -0.5*(a*a/ar-a+0.25*ar);
}

real MHD::RoeFlux(cState<real>& F,
		  const pState<real>& WL,
		  const pState<real>& WR,
		  const Vector3D<real>& normal) {
  static const real sqhalf = 0.7071067811865475244;

  static real s,d,c,c2,un,cs,cf,as,af,bn,sbn;
  static real unl,unr,csl,csr,cfl,cfr,bnl,bnr;
  static real pl,pr,el,er;
  static real sqrtd,B2,V2;

  static Vector3D<real> V,B;
  static Vector3D<real> b,t;

  unl = WL.V*normal; unr = WR.V*normal;
  bnl = WL.B*normal; bnr = WR.B*normal;

  pl = WL.p+0.5*WL.B.sqr(); pr = WR.p+0.5*WR.B.sqr();
  el = 0.5*WL.V.sqr()+(pl-GAMMA_2*INV_GAMMA_1*WL.p)/WL.d;
  er = 0.5*WR.V.sqr()+(pr-GAMMA_2*INV_GAMMA_1*WR.p)/WR.d;

  s  = sqrt(WR.d/WL.d);
  d = s*WL.d; sqrtd = sqrt(d);
  s  = 1./(1.+s); c = 1.-s;
  V  = WL.V*s; V += WR.V*c; un = unl*s+unr*c; V2 = V.sqr();
  B  = WR.B*s; B += WL.B*c; bn = bnr*s+bnl*c; B2 = B.sqr();
  c2 = GAMMA_1*((el+pl/WL.d)*s+(er+pr/WR.d)*c-0.5*V2-B2/d);

  sbn = bn < 0. ? -1. : 1.;

  cf  = 0.5*(c2+B2/d);
  cfl = pl+0.5*GAMMA_2*WL.p;
  cfr = pr+0.5*GAMMA_2*WR.p;
#ifndef NDEBUG
  cs  = fabs(cf -sqrt(fabs(cf*cf-c2*bn*bn/d)));
  csl = fabs(cfl-sqrt(fabs(cfl*cfl-GAMMA*WL.p*bnl*bnl)))/WL.d;
  csr = fabs(cfr-sqrt(fabs(cfr*cfr-GAMMA*WR.p*bnr*bnr)))/WR.d;
#else
  cs  = cf-sqrt(cf*cf-c2*bn*bn/d);
  csl = (cfl-sqrt(cfl*cfl-GAMMA*WL.p*bnl*bnl))/WL.d;
  csr = (cfr-sqrt(cfr*cfr-GAMMA*WR.p*bnr*bnr))/WR.d;
#endif
  cf  = 2.*cf-cs;
  cfl = 2.*cfl/WL.d-csl;
  cfr = 2.*cfr/WR.d-csr;

  if( cf > cs ) {
#ifndef NDEBUG
    af = fabs((c2-cs)/(cf-cs));
    as = sqrt(fabs(1.-af)); af = sqrt(af);
#else
    af = (c2-cs)/(cf-cs);
    as = sqrt(1.-af); af = sqrt(af);
#endif
  }
  else
    af = as = sqhalf;

  c   = sqrt(c2);
  cs  = sqrt(cs);  cf  = sqrt(cf);
  csl = sqrt(csl); cfl = sqrt(cfl);
  csr = sqrt(csr); cfr = sqrt(cfr);

#ifndef NDEBUG
  s = fabs(B2-bn*bn);
#else
  s = B2-bn*bn;
#endif
  if( s > 0 ) {
    b = B; b -= normal*bn; b *= 1./sqrt(s);
  }
  else {
    s = sqrt(normal.x*normal.x+normal.y*normal.y);
    b = s==0. ? Vector3D<real>(sqhalf,sqhalf,0.) :
      Vector3D<real>(normal.y/s,-normal.x/s,0.);
  }

  if( un > 0. ) { // Left flux
    F.d = WL.d*unl;
    F.M = WL.V*F.d; F.M += normal*pl; F.M -= WL.B*bnl;
    F.e = unl*(WL.d*el+pl)-bnl*(WL.V*WL.B);
    F.B = WL.B*unl; F.B -= WL.V*bnl;

    s = EntropyFixMinus(unl-cfl,unr-cfr,un-cf);
    if( s < 0. ) { // Left fast wave
      s *= (0.5/c2)*(af*(WR.p-WL.p-d*cf*(unr-unl))+
		     as*(b*((WR.B-WL.B)*(c*sqrtd)+(WR.V-WL.V)*(d*cs*sbn))));
      F.d += s*af;
      F.M += V*(s*af); F.M -= normal*(s*af*cf); F.M += b*(s*as*cs*sbn);
      F.e += s*(af*(0.5*V2-cf*un+INV_GAMMA_1*c2)+as*(b*(B*(c/sqrtd)+V*(cs*sbn))));
      F.B += b*(s*as*c/sqrtd);

      s = un-sbn*bn/sqrtd;
      if( s < 0. ) { // Left Alfven wave
	t = normal^b;
	s *= 0.5*(t*((WR.V-WL.V)*d+(WR.B-WL.B)*(sbn*sqrtd)));
	F.M += t*s;
	F.e += s*(t*(V+B*(sbn/sqrtd)));
	F.B += t*(s*sbn/sqrtd);

	s = EntropyFixMinus(unl-csl,unr-csr,un-cs);
	if( s < 0. ) { // Left slow wave
	  s *= (0.5/c2)*(as*(WR.p-WL.p-d*cs*(unr-unl))-
			 af*(b*((WR.B-WL.B)*(c*sqrtd)+(WR.V-WL.V)*(d*cf*sbn))));
	  F.d += s*as;
	  F.M += V*(s*as); F.M -= normal*(s*as*cs); F.M -= b*(s*af*cf*sbn);
	  F.e += s*(as*(0.5*V2-cs*un+INV_GAMMA_1*c2)-af*(b*(B*(c/sqrtd)+V*(cf*sbn))));
	  F.B -= b*(s*af*c/sqrtd);

	}
      }
    }
  }
  else { // Right flux
    F.d = WR.d*unr;
    F.M = WR.V*F.d; F.M += normal*pr; F.M -= WR.B*bnr;
    F.e = unr*(WR.d*er+pr)-bnr*(WR.V*WR.B);
    F.B = WR.B*unr; F.B -= WR.V*bnr;

    s = EntropyFixPlus(unl+cfl,unr+cfr,un+cf);
    if( s > 0. ) { // Right fast wave
      s *= (0.5/c2)*(af*(WR.p-WL.p+d*cf*(unr-unl))+
		     as*(b*((WR.B-WL.B)*(c*sqrtd)-(WR.V-WL.V)*(d*cs*sbn))));
      F.d -= s*af;
      F.M -= V*(s*af); F.M -= normal*(s*af*cf); F.M += b*(s*as*cs*sbn);
      F.e -= s*(af*(0.5*V2+cf*un+INV_GAMMA_1*c2)+as*(b*(B*(c/sqrtd)-V*(cs*sbn))));
      F.B -= b*(s*as*c/sqrtd);

      s = un+sbn*bn/sqrtd;
      if( s > 0. ) { // Right Alfven wave
	t = normal^b;
	s *= 0.5*(t*((WR.V-WL.V)*d-(WR.B-WL.B)*(sbn*sqrtd)));
	F.M -= t*s;
	F.e -= s*(t*(V-B*(sbn/sqrtd)));
	F.B += t*(s*sbn/sqrtd);

	s = EntropyFixPlus(unl+csl,unr+csr,un+cs);
	if( s > 0. ) { // Right slow wave
	  s *= (0.5/c2)*(as*(WR.p-WL.p+d*cs*(unr-unl))-
			 af*(b*((WR.B-WL.B)*(c*sqrtd)-(WR.V-WL.V)*(d*cf*sbn))));
	  F.d -= s*as;
	  F.M -= V*(s*as); F.M -= normal*(s*as*cs); F.M -= b*(s*af*cf*sbn);
	  F.e -= s*(as*(0.5*V2+cs*un+INV_GAMMA_1*c2)-af*(b*(B*(c/sqrtd)-V*(cf*sbn))));
	  F.B += b*(s*af*c/sqrtd);

	}
      }
    }
  }

  // Div-B wave
  un = fabs(un);
  s = 0.5*un*(bnr-bnl);
  F.e -= s*bn; F.B -= normal*s;

  return un+cf;
}


real MHD::RoeFlux(cState<real>& F,
		  const cState<real>& UL,
		  const cState<real>& UR,
		  const Vector3D<real>& normal) {
  return RoeFlux(F,UL.Prim(),UR.Prim(),normal);
}
