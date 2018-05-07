#ifndef __EULER3D_INCLUDED__
#define __EULER3D_INCLUDED__

using namespace std;
#include <iostream>
#include <math.h>
#include <vector3d.H>

namespace EULER {

  // Declare adiabatic constant
  extern double GAMMA;
  extern double GAMMA_1;
  extern double INV_GAMMA;
  extern double INV_GAMMA_1;
  extern double GAMMA_INV_GAMMA_1;

  inline void SetGamma(const double& gamma) {
    GAMMA = gamma;
    GAMMA_1 = gamma-1.;
    INV_GAMMA = 1./gamma;
    INV_GAMMA_1 = 1./(gamma-1.);
    GAMMA_INV_GAMMA_1 = gamma/(gamma-1.);
  }

  // Declare conserved variable class
  template <class real> class cState;

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

#ifndef WEAK_ENCAPSULATION
  public:
#endif

    // Possible constructors
    pState() { d = 0.; V.zero(); p = 0.; }
    pState(const real& dd, const Vector3D<real>& VV,
	   const real& pp) : d(dd),V(VV),p(pp) {}
    pState(const real& dd, const real& vx, const real& vy,
	   const real& vz,const real& pp): d(dd),V(vx,vy,vz),p(pp) {}
    pState(const cState<real>&);

    // Casts
    const pState<real>& Prim() const { return *this; }
    const cState<real> Cons() const;

    // Utility functions
    void zero() { d = 0.; V.zero(); p = 0.; }
    int positive() const { return d > 0. && p > 0.; }

    // Public interfaces for private members
    real& Density() { return d; }
    const real& Density() const { return d; }
    Vector3D<real>& Velocity() { return V; }
    const Vector3D<real>& Velocity() const { return V; }
    real& Pressure() { return p; }
    const real& Pressure() const { return p; }

    // Various physical quantities
    const Vector3D<real> Momentum() const { return V*d; }
    real Energy() const { return 0.5*d*V.sqr()+INV_GAMMA_1*p; }
    real Enthalpy() const { return 0.5*d*V.sqr()+GAMMA_INV_GAMMA_1*p; }
    real Temperature() const { return p/d; }
    real SoundSpeed() const { return sqrt(GAMMA*p/d); }
    real SoundSpeedSqr() const { return GAMMA*p/d; }
    real SpecificEnergy() const { return 0.5*V.sqr()+INV_GAMMA_1*p/d; }
    real SpecificEnthalpy() const { return 0.5*V.sqr()+GAMMA_INV_GAMMA_1*p/d; }
    real Entropy() const { return p/pow(d,GAMMA); }

    // Flux functions
    const cState<real> Flux(const Vector3D<real>&) const;
    void SetFlux(cState<real>&,const Vector3D<real>&) const;

    // Overloaded operators
    const pState& operator+() const {
      return *this;
    }
    const pState operator+(const pState& W) const {
      return pState(d+W.d,V+W.V,p+W.p);
    }
    const pState operator-() const {
      return pState(-d,-V,-p);
    }
    const pState operator-(const pState& W) const {
      return pState(d-W.d,V-W.V,p-W.p);
    }
    const pState operator*(const real& a) const {
      return pState(d*a,V*a,p*a);
    }
    const pState operator/(const real& a) const {
      return pState(d/a,V/a,p/a);
    }
    const pState& operator+=(const pState& W) {
      d += W.d; V += W.V; p += W.p; return *this;
    }
    const pState& operator-=(const pState& W) {
      d -= W.d; V -= W.V; p -= W.p; return *this;
    }
    const pState& operator*=(const real& a) {
      d *= a; V *= a; p *= a; return *this;
    }
    const pState& operator/=(const real& a) {
      d /= a; V /= a; p /= a; return *this;
    }
    int operator==(const pState& W) const {
      return d==W.d && V==W.V && p==W.p;
    }
    int operator!=(const pState& W) const {
      return d!=W.d || V!=W.V || p!=W.p;
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

#ifndef WEAK_ENCAPSULATION
  public:
#endif

    // Possible constructors
    cState() { d = 0.; M.zero(); e = 0.; }
    cState(const real& dd, const Vector3D<real>& mm,
	   const real& ee): d(dd),M(mm),e(ee) {}
    cState(const real& dd, const real& mx, const real& my,
	   const real& mz, const real& ee): d(dd),M(mx,my,mz),e(ee) {}
    cState(const pState<real>& W) {
#ifdef WEAK_ENCAPSULATION
      d = W.d; M = W.V*d; e = M.sqr()/(2.*d)+INV_GAMMA_1*W.p;
#else
      d = W.Density(); M = W.Momentum(); e = W.Energy();
#endif
    }

    // Cast to primitive variables
    const cState<real>& Cons() const { return *this; }
    const pState<real> Prim() const {
      return pState<real>(d,M/d,GAMMA_1*(e-M.sqr()/(2.*d)));
    }

    // Utility functions
    void zero() { d = 0.; M.zero(); e = 0.; }
    int positive() const { return d>0. && GAMMA_1*(e-M.sqr()/(2.*d)); }

    // Public interfaces for private members
    real& Density() { return d; }
    const real& Density() const { return d; }
    Vector3D<real>& Momentum() { return M; }
    const Vector3D<real>& Momentum() const { return M; }
    real& Energy() { return e; }
    const real& Energy() const { return e; }

    // Various physical quantities
    const Vector3D<real> Velocity() const { return M/d; }
    real Pressure() const { return GAMMA_1*(e-M.sqr()/(2.*d)); }
    real Enthalpy() const { return GAMMA*e-GAMMA_1*M.sqr()/(2.*d); }
    real Temperature() const { return GAMMA_1*(e-M.sqr()/(2.*d))/d; }
    real SoundSpeed() const { return sqrt(GAMMA*GAMMA_1*(e-M.sqr()/(2.*d))/d); }
    real SoundSpeedSqr() const { return GAMMA*GAMMA_1*(e-M.sqr()/(2.*d))/d; }
    real SpecificEnthalpy() const { return (GAMMA*e-GAMMA_1*M.sqr()/(2.*d))/d; }
    real SpecificEnergy() const { return e/d; }
    real Entropy() const { return GAMMA_1*(e-M.sqr()/(2.*d))/pow(d,GAMMA); }

    // Flux function
    const cState Flux(const Vector3D<real>& n) const {
      real un = (M*n)/d;
      real p  = GAMMA_1*(e-M.sqr()/(2.*d));
      return cState(d*un,M*un+n*p,un*(e+p));
    }
    void SetFlux(cState& F, const Vector3D<real>& n) const {
      F.d  = (M*n)/d;
      F.e  = GAMMA_1*(e-M.sqr()/(2.*d));
      F.M  = M*F.d+n*F.e;
      F.e += e; F.e *= F.d;
      F.d *= d;
    }

    // Overloaded operators
    const cState& operator+() const {
      return *this;
    }
    const cState operator+(const cState& U) const {
      return cState(d+U.d,M+U.M,e+U.e);
    }
    const cState operator-(const cState& U) const {
      return cState(d-U.d,M-U.M,e-U.e);
    }
    const cState operator-() const {
      return cState(-d,-M,-e);
    }
    const cState operator*(const real& a) const {
      return cState(d*a,M*a,e*a);
    }
    const cState operator/(const real& a) const {
      return cState(d/a,M/a,e/a);
    }
    const cState& operator+=(const cState& U) {
      d += U.d; M += U.M; e += U.e; return *this;
    }
    const cState& operator-=(const cState& U) {
      d -= U.d; M -= U.M; e -= U.e; return *this;
    }
    const cState& operator*=(const real& a) {
      d *= a; M *= a; e *= a; return *this;
    }
    const cState& operator/=(const real& a) {
      d /= a; M /= a; e /= a; return *this;
    }
    int operator==(const cState& U) const {
      return d==U.d && M==U.M && e==U.e;
    }
    int operator!=(const cState& U) const {
      return d!=U.d || M!=U.M || e!=U.e;
    }
  };

  template <class real>
  inline pState<real>::pState(const cState<real>& U) {
#ifdef WEAK_ENCAPSULATION
    d = U.d; V = U.M/d; p = GAMMA_1*(U.e-0.5*d*V.sqr());
#else
    d = U.Density(); V = U.Velocity(); p = U.Pressure();
#endif
  }

  template <class real>
  inline const cState<real> pState<real>::Cons() const {
    return cState<real>(d,V*d,0.5*d*V.sqr()+INV_GAMMA_1*p);
  }

  template <class real>
  inline const cState<real> pState<real>::Flux(const Vector3D<real>& n) const {
    real dun = d*(V*n);
    return cState<real>(dun,V*dun+n*p,
			dun*(0.5*V.sqr()+GAMMA_INV_GAMMA_1*p/d));
  }

  template <class real>
  inline void pState<real>::SetFlux(cState<real>& F, const Vector3D<real>& n) const {
    real dun = d*(V*n);
#ifdef WEAK_ENCAPSULATION
    F.d = dun;
    F.M = V*dun+n*p;
    F.e = dun*(0.5*V.sqr()+GAMMA_INV_GAMMA_1*p/d);
#else
    F.Density() = dun;
    F.Momentum() = V*dun+n*p;
    F.Energy() = dun*(0.5*V.sqr()+GAMMA_INV_GAMMA_1*p/d);
#endif
  }

  // Speed of the fastest wave
  template <class real>
  inline real SignalSpeed(const pState<real>& WL, const pState<real>& WR,
			  const Vector3D<real>& n) {
#ifdef WEAK_ENCAPSULATION
    real sl = fabs(WL.V*n)+WL.SoundSpeed();
    real sr = fabs(WR.V*n)+WR.SoundSpeed();
#else
    real sl = fabs(WL.Velocity()*n)+WL.SoundSpeed();
    real sr = fabs(WR.Velocity()*n)+WR.SoundSpeed();
#endif
    return sl > sr ? sl : sr;
  }

  template <class real>
  inline real SignalSpeed(const cState<real>& UL, const cState<real>& UR,
			  const Vector3D<real>& n) {
#ifdef WEAK_ENCAPSULATION
    real sl = fabs(UL.M*n)/UL.d+UL.SoundSpeed();
    real sr = fabs(UR.M*n)/UR.d+UR.SoundSpeed();
#else
    real sl = fabs(UL.Velocity()*n)+UL.SoundSpeed();
    real sr = fabs(UR.Velocity()*n)+UR.SoundSpeed();
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

  extern float GodunovFlux(cState<float>&,const pState<float>&,
			   const pState<float>&,const Vector3D<float>&);
  extern double GodunovFlux(cState<double>&,const pState<double>&,
			    const pState<double>&,const Vector3D<double>&);
  extern float GodunovFlux(cState<float>&,const cState<float>&,
			   const cState<float>&,const Vector3D<float>&);
  extern double GodunovFlux(cState<double>&,const cState<double>&,
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
inline const EULER::pState<real> operator*(const int a,
					   const EULER::pState<real>& W) {
#ifdef WEAK_ENCAPSULATION
  return EULER::pState<real>(W.d*a,W.V*a,W.p*a);
#else
  return EULER::pState<real>(a*W.Density(),W.Velocity()*a,a*W.Pressure());
#endif
}

template <class real>
inline const EULER::pState<real> operator*(const float a,
					   const EULER::pState<real>& W) {
#ifdef WEAK_ENCAPSULATION
  return EULER::pState<real>(W.d*a,W.V*a,W.p*a);
#else
  return EULER::pState<real>(a*W.Density(),W.Velocity()*a,a*W.Pressure());
#endif
}

template <class real>
inline const EULER::pState<real> operator*(const double& a,
					   const EULER::pState<real>& W) {
#ifdef WEAK_ENCAPSULATION
  return EULER::pState<real>(W.d*a,W.V*a,W.p*a);
#else
  return EULER::pState<real>(a*W.Density(),W.Velocity()*a,a*W.Pressure());
#endif
}

template <class real>
inline const EULER::cState<real> operator*(const int a,
					   const EULER::cState<real>& U) {
#ifdef WEAK_ENCAPSULATION
  return EULER::cState<real>(U.d*a,U.M*a,U.e*a);
#else
  return EULER::cState<real>(a*U.Density(),U.Momentum()*a,a*U.Energy());
#endif
}

template <class real>
inline const EULER::cState<real> operator*(const float a,
					   const EULER::cState<real>& U) {
#ifdef WEAK_ENCAPSULATION
  return EULER::cState<real>(U.d*a,U.M*a,U.e*a);
#else
  return EULER::cState<real>(a*U.Density(),U.Momentum()*a,a*U.Energy());
#endif
}

template <class real>
inline const EULER::cState<real> operator*(const double& a,
					   const EULER::cState<real>& U) {
#ifdef WEAK_ENCAPSULATION
  return EULER::cState<real>(U.d*a,U.M*a,U.e*a);
#else
  return EULER::cState<real>(a*U.Density(),U.Momentum()*a,a*U.Energy());
#endif
}

template <class real>
inline ostream& operator<<(ostream& s, const EULER::pState<real>& W) {
#ifdef WEAK_ENCAPSULATION
  return s << W.d << ' ' << W.V << ' ' << W.p;
#else
  return s << W.Density() << ' ' << W.Velocity() << ' ' << W.Pressure();
#endif
}

template <class real>
inline istream& operator>>(istream& s, EULER::pState<real>& W) {
#ifdef WEAK_ENCAPSULATION
  return s >> W.d >> W.V >> W.p;
#else
  return s >> W.Density() >> W.Velocity() >> W.Pressure();
#endif
}

template <class real>
inline ostream& operator<<(ostream& s, const EULER::cState<real>& U) {
#ifdef WEAK_ENCAPSULATION
  return s << U.d << ' ' << U.M << ' ' << U.e;
#else
  return s << U.Density() << ' ' << U.Momentum() << ' ' << U.Energy();
#endif
}

template <class real>
inline istream& operator>>(istream& s, EULER::cState<real>& U) {
#ifdef WEAK_ENCAPSULATION
  return s >> U.d >> U.M >> U.e;
#else
  return s >> U.Density() >> U.Momentum() >> U.Energy();
#endif
}

#endif // __EULER3D_INCLUDED__
