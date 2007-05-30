// -*- C++ -*-
#ifndef HERWIG_DtoKPiPiBaBar_H
#define HERWIG_DtoKPiPiBaBar_H
//
// This is the declaration of the DtoKPiPiBaBar class.
//

#include "Herwig++/Decay/DecayIntegrator.h"
#include "Herwig++/Decay/DecayPhaseSpaceMode.h"
#include "ThePEG/Utilities/Math.h"
#include "Herwig++/Utilities/Interpolator.h"
#include "Herwig++/Utilities/GaussianIntegrator.h"
#include "DtoKPiPiBaBar.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the DtoKPiPiBaBar class.
 *
 * @see \ref DtoKPiPiBaBarInterfaces "The interfaces"
 * defined for DtoKPiPiBaBar.
 */
class DtoKPiPiBaBar: public DecayIntegrator {

/**
 *  The struct for the integration is a friend
 */
friend struct DtoKPiPiBaBarInnerIntegrand;
  
public:

  /**
   * The default constructor.
   */
  DtoKPiPiBaBar();

  /**
   * Which of the possible decays is required
   * @param cc Is this mode the charge conjugate
   * @param dm The decay mode
   */
  virtual int modeNumber(bool & cc,const DecayMode & dm) const;

  /**
   * Return the matrix element squared for a given mode and phase-space channel.
   * @param vertex Output the information on the vertex for spin correlations
   * @param ichan The channel we are calculating the matrix element for. 
   * @param part The decaying Particle.
   * @param decay The particles produced in the decay.
   * @return The matrix element squared for the phase-space configuration.
   */
  double me2(bool vertex, const int ichan,const Particle & part,
	     const ParticleVector & decay) const;

  /**
   * Output the setup information for the particle database
   * @param os The stream to output the information to
   * @param header Whether or not to output the information for MySQL
   */
  virtual void dataBaseOutput(ofstream & os,bool header) const;

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

  /**
   *  Calculate the amplitude for a resonance
   * @param ispin The spin of the intermediate resonance
   * @param gs Whether to use the normal or Gounaris and Sakurai form for the 
   * Breit-Wigner
   * @param mD  The mass of the decaying particle
   * @param mA  The mass of the first  decay product
   * @param mB  The mass of the second decay product
   * @param mC  The mass of the third  decay product
   * @param mAB The mass of the pair AB 
   * @param mAC The mass of the pair AC
   * @param mBC The mass of the pair BC
   * @param mres The on-shell mass of the intermediate resonance
   * @param wres The width         of the intermediate resonance
   */
  inline Complex amplitude(int ispin, bool gs, Energy mD, 
			   Energy mA , Energy mB , Energy mC ,
			   Energy mAB, Energy mAC, Energy mBC,
			   Energy mres, Energy wres) const;

  /**
   *  The \f$F_1(s)\f$ contribution from the scalar \f$\pi\pi\f$ states
   * @param s The mass squared of the \f$\pi\pi\f$ status
   */
  Complex F1(Energy2 s) const;
 
  /**
   *  Integrand for the \f$\rho\f$ function for the \f$4\pi\f$ state
   */
  inline double rho3(Energy2 s,Energy2 s1, Energy2 s2) const;

protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr fullclone() const;
  //@}

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit() throw(InitException);
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<DtoKPiPiBaBar> initDtoKPiPiBaBar;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  DtoKPiPiBaBar & operator=(const DtoKPiPiBaBar &);

private:

  /**
   *  Switch to control which model is used
   */
  unsigned int _imodel;

  /**
   *  Parameters for the K-matrix based fit
   */
  //@{
  /**
   *  The \f$g^\alpha_{\pi^+\pi^-}\f$ coupling for the different resonances in
   *  the K-matrix fit
   */
  vector<Energy> _gpipi;

  /**
   *  The \f$g^\alpha_{K\bar{K}}\f$ coupling for the different resonances in
   *  the K-matrix fit
   */
  vector<Energy> _gKK;

  /**
   *  The \f$g^\alpha_{4\pi}\f$ coupling for the different resonances in
   *  the K-matrix fit
   */
  vector<Energy> _g4pi;

  /**
   *  The \f$g^\alpha_{\eta\eta}\f$ coupling for the different resonances in
   *  the K-matrix fit
   */
  vector<Energy> _getaeta;

  /**
   *  The \f$g^\alpha_{\eta\eta'}\f$ coupling for the different resonances in
   *  the K-matrix fit
   */
  vector<Energy> _getaetap;

  /**
   *  The masses of the resonances
   */
  vector<Energy> _malpha;

  /**
   *  The \f$f^{\rm scatt}_{1\alpha}\f$ parameters for the K-matrix
   */
  vector<double> _fscatt;

  /**
   *  \f$s^{\rm scatt}_0\f$
   */
  Energy2 _s0scatt;

  /**
   *  \f$s_{A_0}\f$
   */
  Energy2 _sA0;

  /**
   * \f$s_A\f$
   */
  Energy2 _sA;

  /**
   *  The \f$\rho_0\f$ parameter for the multi-body function
   */
  Energy2 _rho0;

  /**
   *  \f$g\f$ couplings for easy access
   */
  vector<vector<Energy> > _gialpha;

  /**
   *  Pion mass
   */
  Energy _mpi;

  /**
   *  Kaon mass
   */
  Energy _mK;

  /**
   * \f$\eta\f$ mass
   */
  Energy _meta;

  /**
   *  \f$\eta'\f$ mass
   */
  Energy _metap;
  //@}

  /**
   *  The \f$\beta\f$ production coupling
   */
  //@{
  /**
   *  The real part
   */
  vector<Energy> _betare;
  
  /**
   *  The imaginary part
   */
  vector<Energy> _betaim;
  
  /**
   *  The full coupling
   */
  vector<complex<Energy> > _beta;
  //@}

  /**
   *  The \f$f^{\rm prod}\f$ parameters
   */
  //@{
  /**
   *  The real part
   */
  double _fprodre;

  /**
   *  The imaginary part
   */
  double _fprodim;

  /**
   *  The full coupling
   */
  Complex _fprod;

  /**
   *  \f$s\f$ values for the interpolation table for \f$\rho_3\f$.
   */
  vector<double> _rho3scale;

  /**
   *  \f$\rho_3\f$ values for the interpolation table for \f$\rho_3\f$.
   */
  vector<double> _rho3;

  /**
   *  Initialise the interpolation table for \f$\rho_3\f$
   */
  bool _initrho3;

  /**
   *  The interpolator for \f$\rho_3\f$
   */
  InterpolatorPtr _rho3inter;
  //@}

  /**
   *  The amplitudes and phases for the K-matrix fit
   */
  //@{
  /**
   *  Real part of the amplitude for \f$K^*(892)^-\f$
   */
  double _k892mre;

  /**
   *  Imaginary part of the amplitude for \f$K^*(892)^-\f$
   */
  double _k892mim;

  /**
   *  Real part of the amplitude for \f$K^*_0(1430)^-\f$
   */
  Energy2 _k1430mre0;

  /**
   *  Imaginary part of the amplitude for \f$K^*_0(1430)^-\f$
   */
  Energy2 _k1430mim0;

  /**
   *  Real part of the amplitude for \f$K^*_2(1430)^-\f$
   */
  InvEnergy2 _k1430mre2;

  /**
   *  Imaginary part of the amplitude for \f$K^*_2(1430)^-\f$
   */
  InvEnergy2 _k1430mim2;

  /**
   *  Real part of the amplitude for \f$K^*(1410)^-\f$
   */
  double _k1410mre;

  /**
   *  Imaginary part of the amplitude for \f$K^*(1410)^-\f$
   */
  double _k1410mim;

  /**
   *  Real part of the amplitude for \f$K^*(1680)^-\f$
   */
  double _k1680mre;

  /**
   *  Imaginary part of the amplitude for \f$K^*(1680)^-\f$
   */
  double _k1680mim;

  /**
   *  Real part of the amplitude for \f$K^*(892)^+\f$
   */
  double _k892pre;

  /**
   *  Imaginary part of the amplitude for \f$K^*(892)^+\f$
   */
  double _k892pim;

  /**
   *  Real part of the amplitude for \f$K^*_0(1430)^+\f$
   */
  Energy2 _k1430pre0;

  /**
   *  Imaginary part of the amplitude for \f$K^*_0(1430)^+\f$
   */
  Energy2 _k1430pim0;

  /**
   *  Real part of the amplitude for \f$K^*_2(1430)^+\f$
   */
  InvEnergy2 _k1430pre2;

  /**
   *  Imaginary part of the amplitude for \f$K^*_2(1430)^+\f$
   */
  InvEnergy2 _k1430pim2;

  /**
   *  Real part of the amplitude for \f$\rho(770)\f$
   */
  double _rho770re;

  /**
   *  Imaginary part of the amplitude for \f$\rho(770)\f$
   */
  double _rho770im;

  /**
   *  Real part of the amplitude for \f$\omega(782)\f$
   */
  double _omegare;

  /**
   *  Imaginary part of the amplitude for \f$\omega(782)\f$
   */
  double _omegaim;

  /**
   *  Real part of the amplitude for \f$f_2(1270)\f$
   */
  InvEnergy2 _f2re;

  /**
   *  Imaginary part of the amplitude for \f$f_2(1270)\f$
   */
  InvEnergy2 _f2im;

  /**
   *  Real part of the amplitude for \f$\rho(1450)\f$
   */
  double _rho1450re;

  /**
   *  Imaginary part of the amplitude for \f$\rho(1450)\f$
   */
  double _rho1450im;
  //@}

  /**
   *  The amplitudes and phases for the normal fit
   */
  //@{
  /**
   *  Amplitude for \f$K^*(892)^-\f$
   */
  double _k892mamp;

  /**
   *  Phase for \f$K^*(892)^-\f$
   */
  double _k892mphase;

  /**
   *  Amplitude for \f$K^*_0(1430)^-\f$
   */
  Energy2 _k1430mamp0;

  /**
   *  Phase for \f$K^*_0(1430)^-\f$
   */
  double _k1430mphase0;

  /**
   *  Amplitude for \f$K^*_2(1430)^-\f$
   */
  InvEnergy2 _k1430mamp2;

  /**
   *  Phase for \f$K^*_2(1430)^-\f$
   */
  double _k1430mphase2;

  /**
   *  Amplitude for \f$K^*(1410)^-\f$
   */
  double _k1410mamp;

  /**
   *  Phase for \f$K^*(1410)^-\f$
   */
  double _k1410mphase;

  /**
   *  Amplitude for \f$K^*(1680)^-\f$
   */
  double _k1680mamp;

  /**
   *  Phase for \f$K^*(1680)^-\f$
   */
  double _k1680mphase;

  /**
   *  Amplitude for \f$K^*(892)^+\f$
   */
  double _k892pamp;

  /**
   *  Phase for \f$K^*(892)^+\f$
   */
  double _k892pphase;

  /**
   *  Amplitude for \f$K^*_0(1430)^+\f$
   */
  Energy2 _k1430pamp0;

  /**
   *  Phase for \f$K^*_0(1430)^+\f$
   */
  double _k1430pphase0;

  /**
   *  Amplitude for \f$K^*_2(1430)^+\f$
   */
  InvEnergy2 _k1430pamp2;

  /**
   *  Phase for \f$K^*_2(1430)^+\f$
   */
  double _k1430pphase2;

  /**
   *  Amplitude for \f$\rho(770)\f$
   */
  double _rho770amp;

  /**
   *  Phase for \f$\rho(770)\f$
   */
  double _rho770phase;

  /**
   *  Amplitude for \f$\omega(782)\f$
   */
  double _omegaamp;

  /**
   *  Phase for \f$\omega(782)\f$
   */
  double _omegaphase;

  /**
   *  Amplitude for \f$f_2(1270)\f$
   */
  InvEnergy2 _f2amp;

  /**
   *  Phase for \f$f_2(1270)\f$
   */
  double _f2phase;

  /**
   *  Amplitude for \f$\rho(1450)\f$
   */
  double _rho1450amp;

  /**
   *  Phase for \f$\rho(1450)\f$
   */
  double _rho1450phase;

  /**
   *  Amplitude for \f$f_0(980)\f$
   */
  Energy2 _f980amp;

  /**
   *  Phase for \f$f_0(980)\f$
   */
  double _f980phase;

  /**
   *  Amplitude for \f$f_0(1370)\f$
   */
  Energy2 _f1370amp;

  /**
   *  Phase for \f$f_0(1370)\f$
   */
  double _f1370phase;

  /**
   *  Amplitude for \f$\sigma\f$
   */
  Energy2 _sigmaamp;

  /**
   *  Phase for \f$\sigma\f$
   */
  double _sigmaphase;

  /**
   *  Amplitude for \f$\sigma'\f$
   */
  Energy2 _sigmapamp;

  /**
   *  Phase for \f$\sigma'\f$
   */
  double _sigmapphase;

  /**
   *  Amplitude for the non-resonant component
   */
  double _nonamp;

  /**
   *  Phase for the non-resonant component
   */
  double _nonphase;
  //@}

  /**
   *  Parameters for the Blatt-Weisskopf form-factors
   */
  //@{
  /**
   *  Radial size for the \f$D^0\f$
   */
  InvEnergy _rD0;

  /**
   *  Radial size for the light resonances
   */
  InvEnergy _rres;
  //@}

  /**
   *  Amplitudes as complex numbers
   */
  //@{
  /**
   *  Amplitude for \f$K^*(892)^-$
   */
  Complex _aKm892;

  /**
   *  Amplitude for \f$K^*_0(1430)^-$
   */
  complex<Energy2> _aKm14300;

  /**
   *  Amplitude for \f$K^*_2(1430)^-$
   */
  complex<InvEnergy2> _aKm14302;

  /**
   *  Amplitude for \f$K^*(1410)^-$
   */
  Complex _aKm1410;

  /**
   *  Amplitude for \f$K^*(1680)^-$
   */
  Complex _aKm1680;

  /**
   *  Amplitude for \f$K^*(892)^+$
   */
  Complex _aKp892;

  /**
   *  Amplitude for \f$K^*_0(1430)^+$
   */
  complex<Energy2> _aKp14300;

  /**
   *  Amplitude for \f$K^*_2(1430)^+$
   */
  complex<InvEnergy2> _aKp14302;

  /**
   *  Amplitude for \f$\rho(770)\f$
   */
  Complex _arho770;

  /**
   *  Amplitude for \f$\omega(782)\f$
   */
  Complex _aomega;

  /**
   *  Amplitude for \f$f_2(1270)\f$
   */
  complex<InvEnergy2> _af2;

  /**
   *  Amplitude for \f$\rho(1450)\f$
   */
  Complex _arho1450;

  /**
   *  Amplitude for \f$f_0(980)\f$
   */
  complex<Energy2> _af980;

  /**
   *  Amplitude for \f$f_0(1370)\f$
   */
  complex<Energy2> _af1370;

  /**
   *  Amplitude for \f$\sigma\f$
   */
  complex<Energy2> _asigma;

  /**
   *  Amplitude for \f$\sigma'\f$
   */
  complex<Energy2> _asigmap;

  /**
   *  Amplitude for the non-resonant component
   */
  Complex _aNR;
  //@}

  /**
   *  Masses and widths of the various resonances
   */
  //@{
  /**
   *  Mass of the \f$K^*(892)\f$
   */
  Energy _mK892;

  /**
   *  Mass of the \f$K^*_0(1430)\f$
   */
  Energy _mK14300;

  /**
   *  Mass of the \f$K^*_2(1430)\f$
   */
  Energy _mK14302;

  /**
   *  Mass of the \f$K^*(1410)\f$
   */
  Energy _mK1410;

  /**
   *  Mass of the \f$K^*(1680)\f$
   */
  Energy _mK1680;

  /**
   *  Mass of the \f$\rho(770)\f$
   */
  Energy _mrho770;

  /**
   *  Mass of the \f$\omega(782)\f$
   */
  Energy _momega;

  /**
   *  Mass of the \f$f_2(1270)\f$
   */
  Energy _mf2;

  /**
   *  Mass of the \f$\rho(1450)\f$
   */
  Energy _mrho1450;

  /**
   *  Mass of the \f$f_0(980)\f$
   */
  Energy _mf980;

  /**
   *  Mass of the \f$f_0(1370)\f$
   */
  Energy _mf1370;

  /**
   *  Mass of the \f$\sigma\f$
   */
  Energy _msigma;

  /**
   *  Mass of the \f$\sigma'\f$
   */
  Energy _msigmap;

  /**
   *  Width of the \f$K^*(892)\f$
   */
  Energy _wK892;

  /**
   *  Width of the \f$K^*_0(1430)\f$
   */
  Energy _wK14300;

  /**
   *  Width of the \f$K^*_2(1430)\f$
   */
  Energy _wK14302;

  /**
   *  Width of the \f$K^*(1410)\f$
   */
  Energy _wK1410;

  /**
   *  Width of the \f$K^*(1680)\f$
   */
  Energy _wK1680;

  /**
   *  Width of the \f$\rho(770)\f$
   */
  Energy _wrho770;

  /**
   *  Width of the \f$\omega(782)\f$
   */
  Energy _womega;

  /**
   *  Width of the \f$f_2(1270)\f$
   */
  Energy _wf2;

  /**
   *  Width of the \f$\rho(1450)\f$
   */
  Energy _wrho1450;

  /**
   *  Width of the \f$f_0(980)\f$
   */
  Energy _wf980;

  /**
   *  Width of the \f$f_0(1370)\f$
   */
  Energy _wf1370;

  /**
   *  Width of the \f$\sigma\f$
   */
  Energy _wsigma;

  /**
   *  Width of the \f$\sigma'\f$
   */
  Energy _wsigmap;
  //@}

  /**
   *  Integration parameters
   */
  //@{
  /**
   *  The maximum weight
   */
  double _maxwgt;

  /**
   *  The weights for the different channels
   */
  vector<double> _weights;
  //@}
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of DtoKPiPiBaBar. */
template <>
struct BaseClassTrait<Herwig::DtoKPiPiBaBar,1> {
  /** Typedef of the first base class of DtoKPiPiBaBar. */
  typedef Herwig::DecayIntegrator NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the DtoKPiPiBaBar class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::DtoKPiPiBaBar>
  : public ClassTraitsBase<Herwig::DtoKPiPiBaBar> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig++::DtoKPiPiBaBar"; }
  /**
   * The name of a file containing the dynamic library where the class
   * DtoKPiPiBaBar is implemented. It may also include several, space-separated,
   * libraries if the class DtoKPiPiBaBar depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwWeakCurrents.so HwSMDecay.so"; }
};

/** @endcond */

}

namespace Herwig {

/**
 * A struct for the integration of \f$\rho_3\f$, this is the outer integrand
 */
struct DtoKPiPiBaBarInnerIntegrand {

  /**
   *  The constructor
   */
  inline DtoKPiPiBaBarInnerIntegrand(DtoKPiPiBaBarPtr);

  /**
   * Get the function value
   */
  inline double operator ()(double argument) const;

  /**
   * Set the value of s
   */
  inline void s(Energy2) const;

  /**
   * Set the value of \f$s_1\f$
   */
  inline void s1(Energy2) const;

  /**
   * Mass squared of the system
   */
  mutable Energy _s;

  /**
   * Mass squared of the first rho
   */
  mutable Energy _s1;

  /**
   *  Pointer to the DtoKPiPiBaBar class
   */
  DtoKPiPiBaBarPtr _decayer;
};

/**
 * A struct for the integration of \f$\rho_3\f$, this is the outer integrand
 */
struct DtoKPiPiBaBarOuterIntegrand {

  /**
   *  The constructor
   */
  inline DtoKPiPiBaBarOuterIntegrand(DtoKPiPiBaBarPtr,Energy );

  /**
   * Get the function value
   */
  inline double operator ()(double argument) const;

  /**
   * Set the value of s
   */
  inline void s(Energy2) const;

  /**
   *  Pointer to the DtoKPiPiBaBar class
   */
  DtoKPiPiBaBarPtr _decayer;

  /**
   *  The inner integrand
   */
  DtoKPiPiBaBarInnerIntegrand integral;

  /**
   *  The value of \f$m_\pi\f$
   */
  Energy _mpi;

  /**
   * Mass squared of the system
   */
  mutable Energy _s;

  /**
   *  Integrator
   */
  GaussianIntegrator _integrator;
};

}


#include "DtoKPiPiBaBar.icc"

#endif /* HERWIG_DtoKPiPiBaBar_H */
