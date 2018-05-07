#ifndef __MHD3D_INCLUDED__
#define __MHD3D_INCLUDED__

using namespace std;
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <vector3d.H>

namespace MHD {

  // Declare adiabatic constant
  extern double GAMMA;
  extern double GAMMA_1;
  extern double GAMMA_2;
  extern double INV_GAMMA;
  extern double INV_GAMMA_1;
  extern double GAMMA_INV_GAMMA_1;

  inline void SetGamma(const double& gamma) {
    GAMMA = gamma;
    GAMMA_1 = gamma-1.;
    GAMMA_2 = gamma-2.;
    INV_GAMMA = 1./gamma;
    INV_GAMMA_1 = 1./(gamma-1.);
    GAMMA_INV_GAMMA_1 = gamma/(gamma-1.);
  }

  // Declare variable classes
  template <class real> class pState; // primitive
  template <class real> class cState; // conserved
  template <class real> class rState; // rt

  // Define primitive variable class
  template <class real> class pState{

#ifndef WEAK_ENCAPSULATION
  private:
#else
  public:
#endif

    real           d; // Density
    Vector3D<real> V; // Velocity
    real           p; // Pressure
    Vector3D<real> B; // Magnetic field

#ifndef WEAK_ENCAPSULATION
  public:
#endif

    // Possible constructors
    pState() { d = 0.; V.zero(); p = 0.; B.zero(); }
    pState(const real& dd, const Vector3D<real>& VV,
	   const real& pp, const Vector3D<real>& BB):
      d(dd),V(VV),p(pp),B(BB) {}
    pState(const real& dd,
	   const real& vx, const real& vy, const real& vz,
	   const real& pp,
	   const real& bx, const real& by, const real& bz):
      d(dd),V(vx,vy,vz),p(pp),B(bx,by,bz) {}
    pState(const cState<real>&);

    // Casts
    const pState<real>& Prim() const { return *this; }
    const cState<real> Cons() const;

    // Utility functions
    void zero() { d = 0.; V.zero(); p = 0.; B.zero(); }
    int positive() const { return d>0. && p>0.; }

    // Public interfaces for private members
    real& Density() { return d; }
    const real& Density() const { return d; }

    Vector3D<real>& Velocity() { return V; }
    const Vector3D<real>& Velocity() const { return V; }

    real& Pressure() { return p; }
    const real& Pressure() const { return p; }

    Vector3D<real>& MagneticField() { return B; }
    const Vector3D<real>& MagneticField() const { return B; }

    // Various physical quantities
    const Vector3D<real> Momentum() { return V*d; }

    real Energy() const { return 0.5*(d*V.sqr()+B.sqr())+INV_GAMMA_1*p; }

    real Enthalpy() const { return 0.5*d*V.sqr()+B.sqr()+GAMMA_INV_GAMMA_1*p; }

    //          Sum of relative mass fractions
    // x = -----------------------------------------
    //      Sum of relative concentration fractions
    // Default value is for fully ionized hydrogen plasma

    real Temperature(const real& x=1.) const { return 0.5*x*p/d; }

    real SoundSpeed() const {return sqrt(GAMMA*p/d); }

    real SoundSpeedSqr() const {return GAMMA*p/d; }

    real SpecificEnergy() const { return 0.5*V.sqr()+(0.5*B.sqr()+INV_GAMMA_1*p)/d; }

    real SpecificEnthalpy() const { return 0.5*V.sqr()+(B.sqr()+GAMMA_INV_GAMMA_1*p)/d; }

    real KineticEnergy() const { return 0.5*d*V.sqr(); }

    real MagneticEnergy() const { return 0.5*B.sqr(); }

    real Entropy() const { return p/pow(d,GAMMA); }

    const Vector3D<real> AlfvenVelocity() const { B/sqrt(d); }

    real AlfvenSpeed(const Vector3D<real>& n) const { return fabs(B*n)/sqrt(d); }

    real AlfvenSpeedSqr(const Vector3D<real>& n) const { real Bn = B*n; return Bn*Bn/d; }

    real SlowSpeed(const Vector3D<real>& n) const {
      real aux = 0.5*(GAMMA*p+B.sqr());
      real Bn = B*n;
#ifndef NDEBUG
      return sqrt(fabs(aux-sqrt(fabs(aux*aux-GAMMA*p*Bn*Bn)))/d);
#else
      return sqrt((aux-sqrt(aux*aux-GAMMA*p*Bn*Bn))/d);
#endif
    }

    real SlowSpeedSqr(const Vector3D<real>& n) const {
      real aux = 0.5*(GAMMA*p+B.sqr());
      real Bn = B*n;
#ifndef NDEBUG
      return fabs(aux-sqrt(fabs(aux*aux-GAMMA*p*Bn*Bn)))/d;
#else
      return (aux-sqrt(aux*aux-GAMMA*p*Bn*Bn))/d;
#endif
    }

    real FastSpeed(const Vector3D<real>& n) const {
      real aux = 0.5*(GAMMA*p+B.sqr());
      real Bn = B*n;
#ifndef NDEBUG
      return sqrt((aux+sqrt(fabs(aux*aux-GAMMA*p*Bn*Bn)))/d);
#else
      return sqrt((aux+sqrt(aux*aux-GAMMA*p*Bn*Bn))/d);
#endif
    }

    real FastSpeedSqr(const Vector3D<real>& n) const {
      real aux = 0.5*(GAMMA*p+B.sqr());
      real Bn = B*n;
#ifndef NDEBUG
      return (aux+sqrt(fabs(aux*aux-GAMMA*p*Bn*Bn)))/d;
#else
      return (aux+sqrt(aux*aux-GAMMA*p*Bn*Bn))/d;
#endif
    }

    void SlowFastSpeeds(real& cs, real& cf, const Vector3D<real>& n) const {
      real aux = 0.5*(GAMMA*p+B.sqr()); cf = B*n;
#ifndef NDEBUG
      cs = fabs(aux-sqrt(fabs(aux*aux-GAMMA*p*cf*cf)))/d;
#else
      cs = (aux-sqrt(aux*aux-GAMMA*p*cf*cf))/d;
#endif
      cf = 2.*aux/d-cs;
      cs = sqrt(cs);
      cf = sqrt(cf);
    }

    void SlowFastSpeedsSqr(real& cs, real& cf, const Vector3D<real>& n) const {
      real aux = 0.5*(GAMMA*p+B.sqr()); cf = B*n;
#ifndef NDEBUG
      cs = fabs(aux-sqrt(fabs(aux*aux-GAMMA*p*cf*cf)))/d;
#else
      cs = (aux-sqrt(aux*aux-GAMMA*p*cf*cf))/d;
#endif
      cf = 2.*aux/d-cs;
    }

    // Flux functions
    const cState<real> Flux(const Vector3D<real>&) const;
    void SetFlux(cState<real>&, const Vector3D<real>&) const;
    const cState<real> Source() const;
    void SetSource(cState<real>&) const;

    // Overloaded operators
    const pState& operator+() const {
      return *this;
    }
    const pState operator+(const pState& W) const {
      return pState(d+W.d,V+W.V,p+W.p,B+W.B);
    }
    const pState operator-() const {
      return pState(-d,-V,-p,-B);
    }
    const pState operator-(const pState& W) const {
      return pState(d-W.d,V-W.V,p-W.p,B-W.B);
    }
    const pState operator*(const real& a) const {
      return pState(d*a,V*a,p*a,B*a);
    }
    const pState operator/(const real& a) const {
      return pState(d/a,V/a,p/a,B/a);
    }
    const pState& operator+=(const pState& W) {
      d += W.d; V += W.V; p += W.p; B += W.B; return *this;
    }
    const pState& operator-=(const pState& W) {
      d -= W.d; V -= W.V; p -= W.p; B -= W.B; return *this;
    }
    const pState& operator*=(const real& a) {
      d *= a; V *= a; p *= a; B *= a; return *this;
    }
    const pState& operator/=(const real& a) {
      d /= a; V /= a; p /= a; B /= a; return *this;
    }
    int operator==(const pState& W) const {
      return d==W.d && V==W.V && p==W.p && B==W.B;
    }
    int operator!=(const pState& W) const {
      return d!=W.d || V!=W.V || p!=W.p || B!=W.B;
    }
  };


  // Define conserved variable class
  template <class real> class cState{

#ifndef WEAK_ENCAPSULATION
  private:
#else
  public:
#endif

    real           d; // Density
    Vector3D<real> M; // Momentum
    real           e; // Energy
    Vector3D<real> B; // Magnetic field
                      // population densities & time derivatives of population densities ...
#ifndef WEAK_ENCAPSULATION
  public:
#endif

    cState() { d = 0.; M.zero(); e = 0.; B.zero(); }
    cState(const real& dd, const Vector3D<real>& mm,
	   const real& ee, const Vector3D<real>& bb):
      d(dd),M(mm),e(ee),B(bb) {}
    cState(const real& dd, const real& mx, const real& my, const real& mz,
	   const real& ee, const real& bx, const real& by, const real& bz):
      d(dd),M(mx,my,mz),e(ee),B(bx,by,bz) {}
    cState(const pState<real>& W) {
#ifdef WEAK_ENCAPSULATION
      d = W.d; M = W.V*d; B = W.B; e = 0.5*(M.sqr()/d+B.sqr())+INV_GAMMA_1*W.p;
#else
      d = W.Density(); M = W.Momentum(); e = W.Energy(); B = W.MagneticField();
#endif
    }

    // Cast to primitive variables
    const cState<real>& Cons() const { return *this; }
    const pState<real> Prim() const {
      return pState<real>(d,M/d,GAMMA_1*(e-0.5*(M.sqr()/d+B.sqr())),B);
    }

    // Utility functions
    void zero() { d = 0.; M.zero(); e = 0.; B.zero(); }
    int positive() const { return d>0. && (e-0.5*(M.sqr()/d+B.sqr()))>0.; }

    // Public interfaces for private members
    real& Density() { return d; }
    const real& Density() const { return d; }

    Vector3D<real>& Momentum() { return M; }
    const Vector3D<real>& Momentum() const { return M; }

    real& Energy() { return e; }
    const real& Energy() const { return e; }

    Vector3D<real>& MagneticField() { return B; }
    const Vector3D<real>& MagneticField() const { return B; }

    // Various physical quantities
    const Vector3D<real> Velocity() const { return M/d; }

    real Pressure() const { return GAMMA_1*(e-0.5*(M.sqr()/d+B.sqr())); }

    real Enthalpy() const { return GAMMA*e-0.5*(GAMMA_1*M.sqr()/d+GAMMA_2*B.sqr()); }

    real Temperature(const real& x = 1.) const { return 0.5*x*GAMMA_1*(e-0.5*(M.sqr()/d+B.sqr()))/d; }

    real SoundSpeed() const { return sqrt(GAMMA*GAMMA_1*(e-0.5*(M.sqr()/d+B.sqr()))/d); }

    real SoundSpeedSqr() const { return GAMMA*GAMMA_1*(e-0.5*(M.sqr()/d+B.sqr()))/d; }

    real SpecificEnthalpy() const { return (GAMMA*e-0.5*(GAMMA_1*M.sqr()/d+GAMMA_2*B.sqr()))/d; }

    real SpecificEnergy() const { return e/d; }

    real KineticEnergy() const { return 0.5*M.sqr()/d; }

    real MagneticEnergy() const { return 0.5*B.sqr(); }

    real Entropy() const { return GAMMA_1*(e-0.5*(M.sqr()/d+B.sqr()))/pow(d,GAMMA); }

    const Vector3D<real> AlfvenVelocity() const { B/sqrt(d); }

    real AlfvenSpeed(const Vector3D<real>& n) const { return fabs(B*n)/sqrt(d); }

    real AlfvenSpeedSqr(const Vector3D<real>& n) const { real Bn = B*n; return Bn*Bn/d; }

    real SlowSpeed(const Vector3D<real>& n) const {
      real Bn = B.sqr();
      real p = GAMMA_1*(e-0.5*(M.sqr()/d+Bn));
      real aux = 0.5*(GAMMA*p+Bn); Bn = B*n;
#ifndef NDEBUG
      return sqrt(fabs(aux-sqrt(fabs(aux*aux-GAMMA*p*Bn*Bn)))/d);
#else
      return sqrt((aux-sqrt(aux*aux-GAMMA*p*Bn*Bn))/d);
#endif
    }

    real SlowSpeedSqr(const Vector3D<real>& n) const {
      real Bn = B.sqr();
      real p = GAMMA_1*(e-0.5*(M.sqr()/d+Bn));
      real aux = 0.5*(GAMMA*p+Bn); Bn = B*n;
#ifndef NDEBUG
      return fabs(aux-sqrt(fabs(aux*aux-GAMMA*p*Bn*Bn)))/d;
#else
      return (aux-sqrt(aux*aux-GAMMA*p*Bn*Bn))/d;
#endif
    }

    real FastSpeed(const Vector3D<real>& n) const {
      real Bn = B.sqr();
      real p = GAMMA_1*(e-0.5*(M.sqr()/d+Bn));
      real aux = 0.5*(GAMMA*p+Bn); Bn = B*n;
#ifndef NDEBUG
      return sqrt((aux+sqrt(fabs(aux*aux-GAMMA*p*Bn*Bn)))/d);
#else
      return sqrt((aux+sqrt(aux*aux-GAMMA*p*Bn*Bn))/d);
#endif
    }

    real FastSpeedSqr(const Vector3D<real>& n) const {
      real Bn = B.sqr();
      real p = GAMMA_1*(e-0.5*(M.sqr()/d+Bn));
      real aux = 0.5*(GAMMA*p+Bn); Bn = B*n;
#ifndef NDEBUG
      return (aux+sqrt(fabs(aux*aux-GAMMA*p*Bn*Bn)))/d;
#else
      return (aux+sqrt(aux*aux-GAMMA*p*Bn*Bn))/d;
#endif
    }

    void SlowFastSpeeds(real& cs, real& cf, const Vector3D<real>& n) const {
      real Bn = B.sqr();
      real p = GAMMA_1*(e-0.5*(M.sqr()/d+Bn));
      real aux = 0.5*(GAMMA*p+Bn); Bn = B*n;
#ifndef NDEBUG
      cs = fabs(aux-sqrt(fabs(aux*aux-GAMMA*p*Bn*Bn)))/d;
#else
      cs = (aux-sqrt(aux*aux-GAMMA*p*Bn*Bn))/d;
#endif
      cf = 2.*aux/d-cs;
      cs = sqrt(cs);
      cf = sqrt(cf);
    }

    void SlowFastSpeedsSqr(real& cs, real& cf, const Vector3D<real>& n) const {
      real Bn = B.sqr();
      real p = GAMMA_1*(e-0.5*(M.sqr()/d+Bn));
      real aux = 0.5*(GAMMA*p+B); Bn = B*n;
#ifndef NDEBUG
      cs = fabs(aux-sqrt(fabs(aux*aux-GAMMA*p*Bn*Bn)))/d;
#else
      cs = (aux-sqrt(aux*aux-GAMMA*p*Bn*Bn))/d;
#endif
      cf = 2.*aux/d-cs;
    }

    // Flux functions
    const cState Flux(const Vector3D<real>& n) const {
      real Bn = B*n;
      real un = (M*n)/d;
      real pstar = GAMMA_1*(e-0.5*M.sqr()/d)-0.5*GAMMA_2*B.sqr();
      return cState(d*un,M*un+n*pstar-B*Bn,un*(e+pstar)-(Bn/d)*(M*B),B*un-M*(Bn/d));
    }
    void SetFlux(cState& F, const Vector3D<real>& n) const {
      real Bn = B*n;
      F.d  = (M*n)/d;
      F.e  = GAMMA_1*(e-0.5*M.sqr()/d)-0.5*GAMMA_2*B.sqr();
      F.M  = M*F.d; F.M += n*F.e; F.M -= B*Bn;
      F.e += e; F.e *= F.d; Bn /= d; F.e -= Bn*(M*B);
      F.B  = B*F.d-M*Bn;
      F.d *= d;
    }

    const cState Source() const {
      return cState(0.,B,(M*B)/d,M/d);
    }
    void SetSource(cState& S) const {
      S.d = 0.; S.M = B; S.B = M/d; S.e = B*S.B;
    }

    // Overloaded operators
    const cState& operator+() const {
      return *this;
    }
    const cState operator+(const cState& U) const {
      return cState(d+U.d,M+U.M,e+U.e,B+U.B);
    }
    const cState operator-(const cState& U) const {
      return cState(d-U.d,M-U.M,e-U.e,B-U.B);
    }
    const cState operator-() const {
      return cState(-d,-M,-e,-B);
    }
    const cState operator*(const real& a) const {
      return cState(d*a,M*a,e*a,B*a);
    }
    const cState operator/(const real& a) const {
      return cState(d/a,M/a,e/a,B/a);
    }
    const cState& operator+=(const cState& U) {
      d += U.d; M += U.M; e += U.e; B += U.B; return *this;
    }
    const cState& operator-=(const cState& U) {
      d -= U.d; M -= U.M; e -= U.e; B -= U.B; return *this;
    }
    const cState& operator*=(const real& a) {
      d *= a; M *= a; e *= a; B *= a; return *this;
    }
    const cState& operator/=(const real& a) {
      d /= a; M /= a; e /= a; B /= a; return *this;
    }
    int operator==(const cState& U) const {
      return d==U.d && M==U.M && e==U.e && B==U.B;
    }
    int operator!=(const cState& U) const {
      return d!=U.d || M!=U.M || e!=U.e || B!=U.B;
    }
  };

//
// end RT mod
//

  // Define RT population class
  template <class real> class rState{

#ifndef WEAK_ENCAPSULATION
  private:
#else
  public:
#endif

    int            n; // number of populations
    real           *N; // populations

#ifndef WEAK_ENCAPSULATION
  public:
#endif

    rState(void){ n=0; N=0; }
    void init(real* NN,const int nn){ n=nn; N=NN; };

// Utility functions
    size_t size(void){ return n*sizeof(real); }
    void zero(void){ memset(N,0,n*sizeof(real)); }

// Flux functions (MvN: do we need these?)
/*
    const cState Flux(const Vector3D<real>& n) const {
      real Bn = B*n;
      real un = (M*n)/d;
      real pstar = GAMMA_1*(e-0.5*M.sqr()/d)-0.5*GAMMA_2*B.sqr();
      return cState(d*un,M*un+n*pstar-B*Bn,un*(e+pstar)-(Bn/d)*(M*B),B*un-M*(Bn/d));
    }
    void SetFlux(cState& F, const Vector3D<real>& n) const {
      real Bn = B*n;
      F.d  = (M*n)/d;
      F.e  = GAMMA_1*(e-0.5*M.sqr()/d)-0.5*GAMMA_2*B.sqr();
      F.M  = M*F.d; F.M += n*F.e; F.M -= B*Bn;
      F.e += e; F.e *= F.d; Bn /= d; F.e -= Bn*(M*B);
      F.B  = B*F.d-M*Bn;
      F.d *= d;
    }
    const cState Source() const {
      return cState(0.,B,(M*B)/d,M/d);
    }
    void SetSource(cState& S) const {
      S.d = 0.; S.M = B; S.B = M/d; S.e = B*S.B;
    }

    // Overloaded operators
    const cState& operator+() const {
      return *this;
    }
    const cState operator+(const cState& U) const {
      return cState(d+U.d,M+U.M,e+U.e,B+U.B);
    }
    const cState operator-(const cState& U) const {
      return cState(d-U.d,M-U.M,e-U.e,B-U.B);
    }
    const cState operator-() const {
      return cState(-d,-M,-e,-B);
    }
    const cState operator*(const real& a) const {
      return cState(d*a,M*a,e*a,B*a);
    }
    const cState operator/(const real& a) const {
      return cState(d/a,M/a,e/a,B/a);
    }
    const cState& operator+=(const cState& U) {
      d += U.d; M += U.M; e += U.e; B += U.B; return *this;
    }
    const cState& operator-=(const cState& U) {
      d -= U.d; M -= U.M; e -= U.e; B -= U.B; return *this;
    }
    const cState& operator*=(const real& a) {
      d *= a; M *= a; e *= a; B *= a; return *this;
    }
    const cState& operator/=(const real& a) {
      d /= a; M /= a; e /= a; B /= a; return *this;
    }
    int operator==(const cState& U) const {
      return d==U.d && M==U.M && e==U.e && B==U.B;
    }
    int operator!=(const cState& U) const {
      return d!=U.d || M!=U.M || e!=U.e || B!=U.B;
    }
*/
  };

//
// end RT mod
//

  template <class real>
  inline pState<real>::pState(const cState<real>& U) {
#ifdef WEAK_ENCAPSULATION
    d = U.d; V = U.M/d; B = U.B; p = GAMMA_1*(U.e-0.5*(U.M*V+B.sqr()));
#else
    d = U.Density(); V = U.Velocity(); p = U.Pressure(); B = U.MagneticField();
#endif
  }

  template <class real>
  inline const cState<real> pState<real>::Cons() const {
    return cState<real>(d,V*d,0.5*(d*V.sqr()+B.sqr())+INV_GAMMA_1*p,B);
  }

  template <class real>
  inline const cState<real> pState<real>::Flux(const Vector3D<real>& n) const {
    real dun = d*(V*n);
    real Bn  = B*n;
    real pstar = p+0.5*B.sqr();
    return cState<real>(dun,
			V*dun+n*pstar-B*Bn,
			dun*(0.5*V.sqr()+(2.*pstar-GAMMA_2*INV_GAMMA_1*p)/d)-Bn*(B*V),
			B*(dun/d)-V*Bn);
  }

  template <class real>
  inline void pState<real>::SetFlux(cState<real>& F, const Vector3D<real>& n) const {
    real dun = d*(V*n);
    real Bn  = B*n;
    real pstar = p+0.5*B.sqr();
#ifdef WEAK_ENCAPSULATION
    F.d = dun;
    F.M = V*dun; F.M += n*pstar; F.M -= B*Bn;
    F.e = dun*(0.5*V.sqr()+(2.*pstar-GAMMA_2*INV_GAMMA_1*p)/d)-Bn*(B*V);
    F.B = B*(dun/d); F.B -= V*Bn;
#else
    F.Density() = dun;
    F.Momentum() = V*dun+n*pstar-B*Bn;
    F.Energy() = dun*(0.5*V.sqr()+(2.*pstar-GAMMA_2*INV_GAMMA_1*p)/d)-Bn*(B*V);
    F.MagneticField() = B*(dun/d)-V*Bn;
#endif
  }

  template <class real>
  inline const cState<real> pState<real>::Source() const {
    return cState<real>(0.,B,V*B,V);
  }

  template <class real>
  inline void pState<real>::SetSource(cState<real>& S) const {
#ifdef WEAK_ENCAPSULATION
    S.d = 0.; S.M = B; S.e = V*B; S.B = V;
#else
    S.Density() = 0.; S.Momentum() = B; S.Energy() = V*B; S.MagneticField() = V;
#endif
  }

  template <class real>
  inline real SignalSpeed(const pState<real>& WL, const pState<real>& WR,
			  const Vector3D<real>& n) {
#ifdef WEAK_ENCAPSULATION
    real sl = fabs(WL.V*n)+WL.FastSpeed(n);
    real sr = fabs(WR.V*n)+WR.FastSpeed(n);
#else
    real sl = fabs(WL.Velocity()*n)+WL.FastSpeed(n);
    real sr = fabs(WR.Velocity()*n)+WR.FastSpeed(n);
#endif
    return sl > sr ? sl : sr;
  }

  template <class real>
  inline real SignalSpeed(const cState<real>& UL, const cState<real>& UR,
			  const Vector3D<real>& n) {
#ifdef WEAK_ENCAPSULATION
    real sl = fabs(UL.M*n)/UL.d+UL.FastSpeed(n);
    real sr = fabs(UR.M*n)/UR.d+UR.FastSpeed(n);
#else
    real sl = fabs(UL.Velocity()*n)+UL.FastSpeed(n);
    real sr = fabs(UR.Velocity()*n)+UR.FastSpeed(n);
#endif
    return sl > sr ? sl : sr;
  }

  // Finite-volume flux functions
  extern float RoeFlux(cState<float>&,const pState<float>&,
		       const pState<float>&,const Vector3D<float>&);
  extern double RoeFlux(cState<double>&,const pState<double>&,
			const pState<double>&,const Vector3D<double>&);
  extern float RoeFlux(cState<float>&,const cState<float>&,
		       const cState<float>&,const Vector3D<float>&);
  extern double RoeFlux(cState<double>&,const cState<double>&,
			const cState<double>&,const Vector3D<double>&);

  extern float HLLEFlux(cState<float>&,const pState<float>&,
			const pState<float>&,const Vector3D<float>&);
  extern double HLLEFlux(cState<double>&,const pState<double>&,
			 const pState<double>&,const Vector3D<double>&);
  extern float HLLEFlux(cState<float>&,const cState<float>&,
			const cState<float>&,const Vector3D<float>&);
  extern double HLLEFlux(cState<double>&,const cState<double>&,
			 const cState<double>&,const Vector3D<double>&);

  extern float HLLELFlux(cState<float>&,const pState<float>&,
			 const pState<float>&,const Vector3D<float>&);
  extern double HLLELFlux(cState<double>&,const pState<double>&,
			  const pState<double>&,const Vector3D<double>&);
  extern float HLLELFlux(cState<float>&,const cState<float>&,
			 const cState<float>&,const Vector3D<float>&);
  extern double HLLELFlux(cState<double>&,const cState<double>&,
			  const cState<double>&,const Vector3D<double>&);

  extern float LaxFriedrichsFlux(cState<float>&,const pState<float>&,
				 const pState<float>&,const Vector3D<float>&);
  extern double LaxFriedrichsFlux(cState<double>&,const pState<double>&,
				  const pState<double>&,const Vector3D<double>&);
  extern float LaxFriedrichsFlux(cState<float>&,const cState<float>&,
				 const cState<float>&,const Vector3D<float>&);
  extern double LaxFriedrichsFlux(cState<double>&,const cState<double>&,
				  const cState<double>&,const Vector3D<double>&);

  extern float BoltzmannFlux(cState<float>&,const pState<float>&,
			   const pState<float>&,const Vector3D<float>&);
  extern double BoltzmannFlux(cState<double>&,const pState<double>&,
			    const pState<double>&,const Vector3D<double>&);
  extern float BoltzmannFlux(cState<float>&,const cState<float>&,
			   const cState<float>&,const Vector3D<float>&);
  extern double BoltzmannFlux(cState<double>&,const cState<double>&,
			    const cState<double>&,const Vector3D<double>&);

  // Experimental part
  extern float TestFlux(cState<float>&,const pState<float>&,
			const pState<float>&,const Vector3D<float>&);
  extern double TestFlux(cState<double>&,const pState<double>&,
			 const pState<double>&,const Vector3D<double>&);
  extern float TestFlux(cState<float>&,const cState<float>&,
			const cState<float>&,const Vector3D<float>&);
  extern double TestFlux(cState<double>&,const cState<double>&,
			 const cState<double>&,const Vector3D<double>&);

}

template <class real>
inline const MHD::pState<real> operator*(const int a,
					 const MHD::pState<real>& W) {
#ifdef WEAK_ENCAPSULATION
  return MHD::pState<real>(W.d*a,W.V*a,W.p*a,W.B*a);
#else
  return MHD::pState<real>(a*W.Density(),W.Velocity()*a,
			   a*W.Pressure(),W.MagneticField()*a);
#endif
}

template <class real>
inline const MHD::pState<real> operator*(const float a,
					 const MHD::pState<real>& W) {
#ifdef WEAK_ENCAPSULATION
  return MHD::pState<real>(W.d*a,W.V*a,W.p*a,W.B*a);
#else
  return MHD::pState<real>(a*W.Density(),W.Velocity()*a,
			   a*W.Pressure(),W.MagneticField()*a);
#endif
}

template <class real>
inline const MHD::pState<real> operator*(const double& a,
					 const MHD::pState<real>& W) {
#ifdef WEAK_ENCAPSULATION
  return MHD::pState<real>(W.d*a,W.V*a,W.p*a,W.B*a);
#else
  return MHD::pState<real>(a*W.Density(),W.Velocity()*a,
			   a*W.Pressure(),W.MagneticField()*a);
#endif
}

template <class real>
inline const MHD::cState<real> operator*(const int a,
					 const MHD::cState<real>& U) {
#ifdef WEAK_ENCAPSULATION
  return MHD::cState<real>(U.d*a,U.M*a,U.e*a,U.B*a);
#else
  return MHD::cState<real>(a*U.Density(),U.Momentum()*a,
			   a*U.Energy(),U.MagneticField()*a);
#endif
}

template <class real>
inline const MHD::cState<real> operator*(const float a,
					 const MHD::cState<real>& U) {
#ifdef WEAK_ENCAPSULATION
  return MHD::cState<real>(U.d*a,U.M*a,U.e*a,U.B*a);
#else
  return MHD::cState<real>(a*U.Density(),U.Momentum()*a,
			   a*U.Energy(),U.MagneticField()*a);
#endif
}

template <class real>
inline const MHD::cState<real> operator*(const double& a,
					 const MHD::cState<real>& U) {
#ifdef WEAK_ENCAPSULATION
  return MHD::cState<real>(U.d*a,U.M*a,U.e*a,U.B*a);
#else
  return MHD::cState<real>(a*U.Density(),U.Momentum()*a,
			   a*U.Energy(),U.MagneticField()*a);
#endif
}

template <class real>
inline ostream& operator<<(ostream& s, const MHD::pState<real>& W) {
#ifdef WEAK_ENCAPSULATION
  return s << W.d << ' ' << W.V << ' ' << W.p << ' ' << W.B;
#else
  return s << W.Density()  << ' ' << W.Velocity() << ' '
	   << W.Pressure() << ' ' << W.MagneticField();
#endif
}

template <class real>
inline istream& operator>>(istream& s, MHD::pState<real>& W) {
#ifdef WEAK_ENCAPSULATION
  return s >> W.d >> W.V >> W.p >> W.B;
#else
  return s >> W.Density() >> W.Velocity() >> W.Pressure() >> W.MagneticField();
#endif
}

template <class real>
inline ostream& operator<<(ostream& s, const MHD::cState<real>& U) {
#ifdef WEAK_ENCAPSULATION
  return s << U.d << ' ' << U.M << ' ' << U.e << ' ' << U.B;
#else
  return s << U.Density() << ' ' << U.Momentum() << ' '
	   << U.Energy()  << ' ' << U.MagneticField();
#endif
}

template <class real>
inline istream& operator>>(istream& s, MHD::cState<real>& U) {
#ifdef WEAK_ENCAPSULATION
  return s >> U.d >> U.M >> U.e >> U.B;
#else
  return s >> U.Density() >> U.Momentum() >> U.Energy() >> U.MagneticField();
#endif
}

#endif // __MHD3D_INCLUDED__
