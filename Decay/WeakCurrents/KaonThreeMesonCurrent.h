// -*- C++ -*-
//
// KaonThreeMesonCurrent.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_KaonThreeMesonCurrent_H
#define HERWIG_KaonThreeMesonCurrent_H
//
// This is the declaration of the KaonThreeMesonCurrent class.
//

#include "ThreeMesonCurrentBase.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the KaonThreeMesonCurrent class.
 *
 * @see \ref KaonThreeMesonCurrentInterfaces "The interfaces"
 * defined for KaonThreeMesonCurrent.
 */
class KaonThreeMesonCurrent: public ThreeMesonCurrentBase {

public:

  /**
   * The default constructor.
   */
  KaonThreeMesonCurrent();

  /** @name Methods for the construction of the phase space integrator. */
  //@{  
  /**
   * Complete the construction of the decay mode for integration.
   * This version addes the mesons for the current
   * @param icharge The total charge of the outgoing particles in the current.
   * @param imode   The mode in the current being asked for.
   * @param mode    The phase space mode for the integration
   * @param iloc    The location of the of the first particle from the current in
   *                the list of outgoing particles.
   * @param ires    The location of the first intermediate for the current.
   * @param phase   The prototype phase space channel for the integration.
   * @param upp     The maximum possible mass the particles in the current are
   *                allowed to have.
   * @return Whether the current was sucessfully constructed.
   */
  virtual bool createMode(int icharge,unsigned int imode,DecayPhaseSpaceModePtr mode,
			  unsigned int iloc,unsigned int ires,
			  DecayPhaseSpaceChannelPtr phase,Energy upp);
  //@}

  /**
   * Output the setup information for the particle database
   * @param os The stream to output the information to
   * @param header Whether or not to output the information for MySQL
   * @param create Whether or not to add a statement creating the object
   */
  virtual void dataBaseOutput(ofstream & os,bool header,bool create) const;
  
  /**
   * the matrix element for the \f$a_1\f$ decay to calculate the running width
   * @param q2 The mass of the decaying off-shell \f$a_1\f$, \f$q^2\f$.
   * @param s3 The invariant mass squared of particles 1 and 2, \f$s_3=m^2_{12}\f$.
   * @param s2 The invariant mass squared of particles 1 and 3, \f$s_2=m^2_{13}\f$.
   * @param s1 The invariant mass squared of particles 2 and 3, \f$s_1=m^2_{23}\f$.
   * @return The matrix element squared summed over spins.
   */
  double threeBodyMatrixElement(const int ,  const Energy2 q2,
				const Energy2 s3, const Energy2 s2, 
				const Energy2 s1, const Energy  , 
				const Energy  , const Energy  ) const {
    Energy2 mpi2(sqr(_mpi));
    Complex propb(Trho1(s1,-1)),propa(Trho1(s2,-1)); 
    // the matrix element
    Energy2 output(0.*MeV2); 
    // first resonance
    output+= ((s1-4.*mpi2)+0.25*(s3-s2)*(s3-s2)/q2)*real(propb*conj(propb)); 
    // second resonance
    output+= ((s2-4.*mpi2)+0.25*(s3-s1)*(s3-s1)/q2)*real(propa*conj(propa)); 
    // the interference term 
    output+= (0.5*q2-s3-0.5*mpi2+0.25*(s3-s2)*(s3-s1)/q2)*real(propa*conj(propb)+
							       propb*conj(propa)); 
    return output / sqr(_rho1mass[0]);
  }

protected:

  /**
   * Can the current handle a particular set of mesons. 
   * As this current includes all the allowed modes this is always true.
   */
  virtual bool acceptMode(int) const;

  /**
   * Calculate the form factor for the current. Implements the form factors
   * described above.
   * @param ichan The phase space channel
   * @param imode The mode
   * @param q2 The scale \f$q^2\f$ for the current.
   * @param s1 The invariant mass squared of particles 2 and 3, \f$s_1=m^2_{23}\f$.
   * @param s2 The invariant mass squared of particles 1 and 3, \f$s_2=m^2_{13}\f$.
   * @param s3 The invariant mass squared of particles 1 and 2, \f$s_3=m^2_{12}\f$.
   */
  virtual FormFactors calculateFormFactors(const int ichan,const int imode,
					   Energy2 q2,Energy2 s1,Energy2 s2,Energy2 s3) const;

public:

  /** @name Functions used by the persistent I/O system. */
  //@{
  /**
   * Function used to write out object persistently.
   * @param os the persistent output stream written to.
   */
  void persistentOutput(PersistentOStream & os) const;

  /**
   * Function used to read in object persistently.
   * @param is the persistent input stream read from.
   * @param version the version number of the object when written.
   */
  void persistentInput(PersistentIStream & is, int version);
  //@}

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();

protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  virtual IBPtr clone() const {return new_ptr(*this);}

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const {return new_ptr(*this);}
  //@}

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving and
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit() throw(InitException);

  /**
   * Initialize this object to the begining of the run phase.
   */
  virtual void doinitrun();

  /**
   * Check sanity of the object during the setup phase.
   */
  virtual void doupdate() throw(UpdateException);
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<KaonThreeMesonCurrent> initKaonThreeMesonCurrent;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  KaonThreeMesonCurrent & operator=(const KaonThreeMesonCurrent &);

private:

  /**
   *  The \f$\rho\f$ lineshape for the axial-vector terms
   * @param q2 The scale \f$q^2\f$ for the lineshape
   * @param ires Which \f$\rho\f$ multiplet
   */
  Complex Trho1(Energy2 q2,int ires) const {
    Complex output(0.);
    double norm(0.);
    for(unsigned int ix=0,N=_rho1wgts.size();ix<N;++ix) norm+=_rho1wgts[ix];
    if(ires<0) {
      for(unsigned int ix=0,N=_rho1wgts.size();ix<N;++ix) {
	output+=_rho1wgts[ix]*BWrho1(q2,ix);
      }
    }
    else {
      unsigned int temp(ires);
      if(temp<_rho1wgts.size()) output=_rho1wgts[temp]*BWrho1(q2,temp);
    }
    return output/norm;
  }
  
  /**
   *  The \f$\rho\f$ lineshape for the vector terms
   * @param q2 The scale \f$q^2\f$ for the lineshape
   * @param ires Which \f$\rho\f$ multiplet
   */
  Complex Trho2(Energy2 q2,int ires) const {
    Complex output(0.);
    double norm(0.);
    for(unsigned int ix=0,N=_rho2wgts.size();ix<N;++ix) norm+=_rho2wgts[ix];
    if(ires<0) {
      for(unsigned int ix=0,N=_rho2wgts.size();ix<N;++ix) {
	output+=_rho2wgts[ix]*BWrho2(q2,ix);
      }
    }
    else {
      unsigned int temp(ires);
      if(temp<_rho2wgts.size()) output=_rho2wgts[temp]*BWrho2(q2,temp);
    }
    return output/norm;
  }
  
  /**
   *  The \f$K^*\f$ lineshape for the axial-vector terms
   * @param q2 The scale \f$q^2\f$ for the lineshape
   * @param ires Which \f$K^*\f$ multiplet
   */
  Complex TKstar1(Energy2 q2,int ires) const  {
    Complex output(0.);
    double norm(0.);
    for(unsigned int ix=0,N=_kstar1wgts.size();ix<N;++ix) norm+=_kstar1wgts[ix];
    if(ires<0) {
      for(unsigned int ix=0,N=_kstar1wgts.size();ix<N;++ix) {
	output+=_kstar1wgts[ix]*BWKstar1(q2,ix);
      }
    }
    else {
      unsigned int temp(ires);
      if(temp<_kstar1wgts.size()) output=_kstar1wgts[temp]*BWKstar1(q2,temp);
    }
    return output/norm;
  }
  
  /**
   *  The \f$\rho\f$ lineshape for the vector terms
   * @param q2 The scale \f$q^2\f$ for the lineshape
   * @param ires Which \f$K^*\f$ multiplet
   */
  Complex TKstar2(Energy2 q2,int ires) const {
    Complex output(0.);
    double norm(0.);
    for(unsigned int ix=0,N=_kstar2wgts.size();ix<N;++ix) norm+=_kstar2wgts[ix];
    if(ires<0) {
      for(unsigned int ix=0,N=_kstar2wgts.size();ix<N;++ix) {
	output+=_kstar2wgts[ix]*BWKstar2(q2,ix);
      }
    }
    else {
      unsigned int temp(ires);
      if(temp<_kstar2wgts.size()) output=_kstar2wgts[temp]*BWKstar2(q2,temp);
    }
    return output/norm;
  }

  /**
   *  The \f$\rho\f$ Breit-Wigner for the axial-vector terms
   * @param q2 The scale \f$q^2\f$ for the Breit-Wigner
   * @param ires Which \f$\rho\f$ multiplet
   */
  Complex BWrho1(Energy2 q2, unsigned int ires) const {
    if(ires>=_rho1mass.size()) return 0.;
    Energy mass  = _rho1mass [ires];
    Energy width = _rho1width[ires];
    Energy q=sqrt(q2);
    Energy pcm0 = Kinematics::pstarTwoBodyDecay(mass,_mpi,_mpi);
    Energy pcm  = q<=2.*_mpi ? 0.*GeV : Kinematics::pstarTwoBodyDecay(q,_mpi,_mpi);
    double ratio = Math::powi(pcm/pcm0, 3);
    Energy gam(width*mass*ratio/q);
    return sqr(mass)/(sqr(mass)-q2-Complex(0.,1.)*mass*gam);
  } 

  /**
   *  The \f$\rho\f$ Breit-Wigner for the vector terms
   * @param q2 The scale \f$q^2\f$ for the Breit-Wigner
   * @param ires Which \f$\rho\f$ multiplet
   */
  Complex BWrho2(Energy2 q2, unsigned int ires) const {
    if(ires>=_rho2mass.size()) return 0.;
    Energy mass  = _rho2mass [ires];
    Energy width = _rho2width[ires];
    Energy q=sqrt(q2);
    Energy pcm0 = Kinematics::pstarTwoBodyDecay(mass,_mpi,_mpi);
    Energy pcm  = q<=2.*_mpi ? 0.*GeV : Kinematics::pstarTwoBodyDecay(q,_mpi,_mpi);
    double ratio(pcm/pcm0);ratio*=ratio*ratio;
    Energy gam(width*mass*ratio/q);
    return sqr(mass)/(sqr(mass)-q2-Complex(0.,1.)*mass*gam);
  }
  
  /**
   * The \f$K^*\f$ Breit-Wigner for the axial-vector terms
   * @param q2 The scale \f$q^2\f$ for the Breit-Wigner
   * @param ires Which \f$K^*\f$ multiplet
   */
  Complex BWKstar1(Energy2 q2, unsigned int ires) const {
    if(ires>=_kstar1mass.size()) return 0.;
    Energy mass  = _kstar1mass [ires];
    Energy width = _kstar1width[ires];
    Energy q=sqrt(q2);
    Energy pcm0 = Kinematics::pstarTwoBodyDecay(mass,_mK,_mpi);
    Energy pcm  = q<=_mpi+_mK ? 0.*GeV : Kinematics::pstarTwoBodyDecay(q,_mK,_mpi);
    double ratio(pcm/pcm0);ratio*=ratio*ratio;
    Energy gam(width*mass*ratio/q);
    return sqr(mass)/(sqr(mass)-q2-Complex(0.,1.)*mass*gam);
  }

  /**
   * The \f$K^*\f$ Breit-Wigner for the vector terms
   * @param q2 The scale \f$q^2\f$ for the Breit-Wigner
   * @param ires Which \f$K^*\f$ multiplet
   */
  Complex BWKstar2(Energy2 q2, unsigned int ires) const  {
    if(ires>=_kstar2mass.size()) return 0.;
    Energy mass  = _kstar2mass [ires];
    Energy width = _kstar2width[ires];
    Energy q=sqrt(q2);
    Energy pcm0 = Kinematics::pstarTwoBodyDecay(mass,_mK,_mpi);
    Energy pcm  = q<=_mpi+_mK ? 0.*GeV : Kinematics::pstarTwoBodyDecay(q,_mK,_mpi);
    double ratio(pcm/pcm0);ratio*=ratio*ratio;
    Energy gam(width*mass*ratio/q);
    return sqr(mass)/(sqr(mass)-q2-Complex(0.,1.)*mass*gam);
  }

  /**
   * \f$a_1\f$ Breit-Wigner
   * @param q2 The scale \f$q^2\f$ for the Breit-Wigner
   * @return The Breit-Wigner
   */
  Complex a1BreitWigner(Energy2 q2) const {
    Complex ii(0.,1.);
    Energy2 m2(_a1mass*_a1mass);
    Energy  q(sqrt(q2));
    return m2/(m2-q2-ii*q*a1Width(q2));
  }
  
  /**
   *  The \f$K_1\f$ line shape
   * @param q2 The scale \f$q^2\f$ for the Breit-Wigner
   * @param iopt Whether this is \f$K^*\pi\f$ or \f$\rho K\f$.
   * @param ires the resonance
   */
  Complex TK1(Energy2 q2,unsigned int iopt,int ires) const {
    Complex denom(0),num(0.);
    if(iopt==0) {
      for(unsigned int ix=0;ix<_k1wgta.size();++ix) denom+=_k1wgta[ix];
      if(ires==-1) {
	for(unsigned int ix=0;ix<_k1wgta.size();++ix) 
	  num+=_k1wgta[ix]*K1BreitWigner(q2,ix);
      }
      else {	
	num+=_k1wgta[ires]*K1BreitWigner(q2,ires);
      }
    }
    else if(iopt==1) {
      for(unsigned int ix=0;ix<_k1wgtb.size();++ix) denom+=_k1wgtb[ix];
      if(ires==-1) {
	for(unsigned int ix=0;ix<_k1wgtb.size();++ix) 
	  num+=_k1wgtb[ix]*K1BreitWigner(q2,ix);
      }
      else {	
	num+=_k1wgtb[ires]*K1BreitWigner(q2,ires);
      }
    }
    else {
      return 0.;
    }
    return num/denom;
  }

  /**
   * The \f$K_1\f$ Breit-Wigner
   * @param q2 The scale \f$q^2\f$ for the Breit-Wigner
   * @param ires the resonance
   * @return The Breit-Wigner
   */
  Complex K1BreitWigner(Energy2 q2,unsigned int ires) const {
    if(ires>=_k1mass.size()) return 0.;
    Energy2 m2=sqr(_k1mass[ires]),mg=_k1mass[ires]*_k1width[ires];
    return (-m2+Complex(0.,1.)*mg)/(q2-m2+Complex(0.,1.)*mg);
  }

  /**
   * The \f$a_1\f$ running width
   * @param q2 The scale \f$q^2\f$ for the Breit-Wigner
   * @return The \f$a_1\f$ running width.
   */
  Energy a1Width(Energy2 q2) const {
    if(!_a1opt) return _a1mass*_a1width*g(q2)/g(_a1mass*_a1mass)/sqrt(q2);
    else        return (*_a1runinter)(q2);
  }
  
  /**
   *  The \f$g(Q^2)\f$ function of Kuhn and Santamaria
   */
  double g(Energy2 q2) const {
    double output;
    if(q2<9.*_mpi*_mpi) {
      output=0.;
    }
    else if(q2<sqr(_rho1mass[0]+_mpi)) {
      double diff=(q2-9.*_mpi*_mpi)/GeV2;
      
      output=4.1*sqr(diff)*diff*(1.-3.3*diff+5.8*sqr(diff));
    }
    else {
      double ratio = q2/GeV2;
      output = ratio*(1.623+10.38/ratio-9.32/sqr(ratio)+0.65/(ratio*sqr(ratio)));
    }
    return output;
  }
  
  /**
   * Initialize the \f$a_1\f$ running width
   * @param iopt Initialization option (-1 full calculation, 0 set up the interpolation)
   */
  void inita1Width(int iopt);

  /**
   *  The \f$T_\omega\f$ function
   * @param q2 The scale
   * @param ires the resonance
   */
  Complex Tomega(Energy2 q2, int ires) const {
    double denom=(1.+_epsomega);
    Complex num(0.);
    if(ires<0) num=OmegaPhiBreitWigner(q2,0)+_epsomega*OmegaPhiBreitWigner(q2,1);
    else if(ires==0) num=OmegaPhiBreitWigner(q2,0);
    else             num=OmegaPhiBreitWigner(q2,1);
    return num/denom;
  }

  /**
   *  The \f$\omega\f$ and \f$\phi\f$ Breit-Wigner
   * @param q2 The scale
   * @param ires the resonance
   */
  Complex OmegaPhiBreitWigner(Energy2 q2, unsigned int ires) const {
    Energy2 m2,mg;
    if(ires==0) {
      m2=sqr(_omegamass);
      mg=_omegamass*_omegawidth;
    }
    else {
      m2=sqr(_phimass);
      mg=_phimass*_phiwidth;
    }
    return (-m2+Complex(0.,1.)*mg)/(q2-m2+Complex(0.,1.)*mg);
  }

  /**
   * The \f$\omega-\phi\f$ \f$K^*\f$ form-factor for the \f$F_5\f$ form-factor
   * @param s1 The scale \f$s_1\f$.
   * @param s2 The scale \f$s_2\f$.
   * @param ires Which resonances to use
   * @return The mixed Breit-Wigner
   */
  Complex TOmegaKStar(Energy2 s1,Energy2 s2,int ires) const {
    Complex output;
    if(ires<0)         output = _omegaKstarwgt*TKstar1(s1,-1)+Tomega(s2,-1);
    else if(ires%2==0) output = _omegaKstarwgt*TKstar1(s1,ires/2);
    else if(ires%2==1) output = Tomega(s2,ires/2);
    return output/(1.+_omegaKstarwgt);
  }

private:

  /**
   *  Parameters for the \f$\rho\f$ in the axial-vector terms
   */
  //@{
  /**
   *  Weight for the different resonances
   */
  vector<double> _rho1wgts;

  /**
   *  Masses
   */
  vector<Energy> _rho1mass;

  /**
   *  Widths
   */
  vector<Energy> _rho1width;
  //@}

  /**
   *  Parameters for the \f$\rho\f$ in the vector terms
   */
  //@{
  /**
   *  Weight for the different resonances
   */
  vector<double> _rho2wgts;

  /**
   *  Masses
   */
  vector<Energy> _rho2mass;

  /**
   *  Widths
   */
  vector<Energy> _rho2width;
  //@}

  /**
   *  Parameters for the \f$K^*\f$ in the axial-vector terms
   */
  //@{
  /**
   *  Weight for the different resonances
   */
  vector<double> _kstar1wgts;

  /**
   *  Masses
   */
  vector<Energy> _kstar1mass;

  /**
   *  Widths
   */
  vector<Energy> _kstar1width;
  //@}

  /**
   *  Parameters for the \f$K^*\f$ in the vector terms
   */
  //@{
  /**
   *  Weight for the different resonances
   */
  vector<double> _kstar2wgts;

  /**
   *  Masses
   */
  vector<Energy> _kstar2mass;

  /**
   *  Widths
   */
  vector<Energy> _kstar2width;
  //@}

  /**
   *  Parameters for the three meson resonances
   */
  //@{
  /**
   * The mass of the \f$a_1\f$ resonances.
   */
  Energy _a1mass;
  
  /**
   * The width of the \f$a_1\f$ resonances.
   */
  Energy _a1width;

  /**
   *  The masses of the \f$aK1\f$ resonances.
   */
  vector<Energy> _k1mass;

  /**
   *  The widths of the \f$K_1\f$ resonances.
   */
  vector<Energy> _k1width;

  /**
   *  The weights for the different \f$K_1\f$ resonances for \f$K_1\to K^*\pi\f$
   */
  vector<double> _k1wgta;

  /**
   *  The weights for the different \f$K_1\f$ resonaces for \f$K_1\to\rho K\f$.
   */
  vector<double> _k1wgtb;

  /**
   * The \f$a_1\f$ width for the running \f$a_1\f$ width calculation.
   */
  vector<Energy>  _a1runwidth;

  /**
   * The \f$q^2\f$ for the running \f$a_1\f$  width calculation.
   */
  vector<Energy2> _a1runq2;

  /**
   * The interpolator for the running \f$a_1\f$ width calculation.
   */
  Interpolator<Energy,Energy2>::Ptr _a1runinter;
  //@}

  /**
   *  Parameters for the \f$T_\omega\f$ function
   */
  //@{
  /**
   *  Mixing parameter
   */
  double _epsomega;

  /**
   *  Mass of the \f$\omega\f$
   */
  Energy _omegamass;

  /**
   *  Width of the \f$\omega\f$
   */
  Energy _omegawidth;

  /**
   *  Mass of the \f$\phi\f$
   */
  Energy _phimass;

  /**
   *  Width of the \f$\phi\f$
   */
  Energy _phiwidth;
  //@}

  /**
   * The relative weight of the \f$\omega-\phi\f$ and \f$K^*\f$ where needed.
   */
  double _omegaKstarwgt;

  /**
   * The pion decay constant, \f$f_\pi\f$.
   */
  Energy _fpi;

  /**
   * The pion mass
   */
  Energy _mpi;

  /**
   * The kaon mass
   */
  Energy _mK;

  /**
   *  Initialization switches
   */
  //@{
  /**
   * Initialize the running \f$a_1\f$ width.
   */
  bool _initializea1;
  
  /**
   * use local values of the \f$\rho\f$ masses and widths
   */
  bool _rhoparameters;

  /**
   * use local values of the \f$K^*\f$ resonances masses and widths
   */
  bool _kstarparameters;

  /**
   * Use local values of the \f$a_1\f$ parameters
   */
  bool _a1parameters;

  /**
   * Use local values of the \f$K_1\f$ parameters
   */
  bool _k1parameters;

  /**
   * Option for the \f$a_1\f$ width
   */
  bool _a1opt;

  /**
   * Option for the \f$\omega-\phi\f$ parameters
   */
  bool _omegaopt;
  //@}

  /**
   *  The maximum mass of the hadronic system
   */
  Energy _maxmass;

  /**
   *  The maximum mass when the running width was calculated
   */
  Energy _maxcalc;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of KaonThreeMesonCurrent. */
template <>
struct BaseClassTrait<Herwig::KaonThreeMesonCurrent,1> {
  /** Typedef of the first base class of KaonThreeMesonCurrent. */
  typedef Herwig::ThreeMesonCurrentBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the KaonThreeMesonCurrent class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::KaonThreeMesonCurrent>
  : public ClassTraitsBase<Herwig::KaonThreeMesonCurrent> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::KaonThreeMesonCurrent"; }
  /**
   * The name of a file containing the dynamic library where the class
   * KaonThreeMesonCurrent is implemented. It may also include several, space-separated,
   * libraries if the class KaonThreeMesonCurrent depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwWeakCurrents.so"; }
};

/** @endcond */

}

#endif /* HERWIG_KaonThreeMesonCurrent_H */
