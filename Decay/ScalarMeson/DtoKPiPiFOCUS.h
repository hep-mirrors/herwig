// -*- C++ -*-
#ifndef HERWIG_DtoKPiPiFOCUS_H
#define HERWIG_DtoKPiPiFOCUS_H
//
// This is the declaration of the DtoKPiPiFOCUS class.
//

#include "Herwig++/Decay/DecayIntegrator.h"
#include "Herwig++/Decay/DecayPhaseSpaceMode.h"
#include "DtoKPiPiFOCUS.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the DtoKPiPiFOCUS class.
 *
 * @see \ref DtoKPiPiFOCUSInterfaces "The interfaces"
 * defined for DtoKPiPiFOCUS.
 */
class DtoKPiPiFOCUS: public DecayIntegrator {

public:

  /**
   * The default constructor.
   */
  DtoKPiPiFOCUS();

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
   *  Methods to calculate the amplitudes for a given channel
   */
  //@{
  /**
   * Calculate the decay angle for the amplitude, the angle is the 
   * angle between the 2 and 3 for the decay $D\to(12)3$ in the rest frame of
   * the resonance which decays to 1 and 2.
   * @param pparent The momentum of the parent
   * @param pres    The momentum of the resonance
   * @param p1      The momentum of the first decay product of the resonance
   * @param cost    The cosine of the angle between 1 and 3
   * @param ac      The product of the magnitudes of the three-momenta of 1 and 3.
   */
  inline void decayAngle(const Lorentz5Momentum & pparent,
			 const Lorentz5Momentum & pres,
			 const Lorentz5Momentum & p1,
			 double & cost, Energy2 & ac) const;

  /**
   * Calculate the amplitude for a resonance
   * @param ispin The spin of the intermediate resonance
   * @param mD  The mass of the decaying particle
   * @param mA  The mass of the first  decay product
   * @param mB  The mass of the second decay product
   * @param mC  The mass of the third  decay product
   * @param mAB The mass of the pair AB 
   * @param mres The on-shell mass of the intermediate resonance
   * @param wres The width         of the intermediate resonance
   * @param ac   The product \f$|a||c|\f$
   * @param cost The cosine of the angle between a and c.
   */
  inline Complex amplitude(int ispin, Energy mD, 
			   Energy mA , Energy mB , Energy mC ,
			   Energy mAB,Energy mres, Energy wres,
			   Energy2 ac, double cost) const;

  /**
   *  Return the \f$K-matrix\f$ result for the \f$s\f$-wave contribution
   * @param s The mass squared of the \f$K\pi\f$ system
   */
  Complex F(Energy2 s);
  //@}

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
  static ClassDescription<DtoKPiPiFOCUS> initDtoKPiPiFOCUS;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  DtoKPiPiFOCUS & operator=(const DtoKPiPiFOCUS &);

private:

  /**
   *  The choice of model 
   */
  unsigned int _imodel;

  /**
   *  Amplitudes and Phases for the isobar model
   */
  //@{
  /**
   *  Amplitude for the non-resonant component
   */
  double _a1NR;

  /**
   *  Phase for the non-resonant component
   */
  double _phi1NR;

  /**
   *  Amplitude for the \f$\kappa\f$ component
   */
  Energy2 _a1kappa;

  /**
   *  Phase for the \f$\kappa\f$ component
   */
  double _phi1kappa;

  /**
   *  Amplitude for the \f$K^*(892)\f$
   */
  double _a1K892;

  /**
   *  Phase for the \f$K^*(892)\f$
   */
  double _phi1K892;

  /**
   *  Amplitude for the \f$K^*(1410)\f$
   */
  double _a1K1410;

  /**
   *  Phase for the \f$K^*(1410)\f$
   */
  double _phi1K1410;

  /**
   *  Amplitude for the \f$K^*(1680)\f$
   */
  double _a1K1680;

  /**
   *  Phase for the \f$K^*(1680)\f$
   */
  double _phi1K1680;

  /**
   *  Amplitude for the \f$K^*_0(1430)\f$
   */
  Energy2 _a1K14300;

  /**
   *  Phase for the \f$K^*_0(1430)\f$
   */
  double _phi1K14300;

  /**
   *  Amplitude for the \f$K^*_2(1430)\f$
   */
  InvEnergy2 _a1K14302;

  /**
   *  Phase for the \f$K^*_2(1430)\f$
   */
  double _phi1K14302;
  //@}

  /**
   *  Amplitudes, Phases and other parameters for the K-matrix model
   */
  //@{
  /**
   *  Amplitude for the \f$K^*(892)\f$
   */
  double _a2K892;

  /**
   *  Phase for the \f$K^*(892)\f$
   */
  double _phi2K892;

  /**
   *  Amplitude for the \f$K^*(1410)\f$
   */
  double _a2K1410;

  /**
   *  Phase for the \f$K^*(1410)\f$
   */
  double _phi2K1410;

  /**
   *  Amplitude for the \f$K^*(1680)\f$
   */
  double _a2K1680;

  /**
   *  Phase for the \f$K^*(1680)\f$
   */
  double _phi2K1680;

  /**
   *  Amplitude for the \f$K^*_2(1430)\f$
   */
  InvEnergy2 _a2K14302;

  /**
   *  Phase for the \f$K^*_2(1430)\f$
   */
  double _phi2K14302;
  //@}

  /**
   *  Parameters for ther K-matrix
   */
  //@{
  /**
   *  The coupling \f$g_1\f$
   */
  Energy _g1;

  /**
   *  The coupling \f$g_2\f$
   */
  Energy _g2;

  /**
   *  The coefficent \f$C_{110}\f$
   */
  double _c110;

  /**
   *  The coefficent \f$C_{111}\f$
   */
  double _c111;

  /**
   *  The coefficent \f$C_{112}\f$
   */
  double _c112;

  /**
   *  The coefficent \f$C_{120}\f$
   */
  double _c120;

  /**
   *  The coefficent \f$C_{121}\f$
   */
  double _c121;

  /**
   *  The coefficent \f$C_{122}\f$
   */
  double _c122;

  /**
   *  The coefficent \f$C_{220}\f$
   */
  double _c220;

  /**
   *  The coefficent \f$C_{221}\f$
   */
  double _c221;

  /**
   *  The coefficent \f$C_{222}\f$
   */
  double _c222;

  /**
   *  The coefficent \f$D_{110}\f$
   */
  double _d110;

  /**
   *  The coefficent \f$D_{111}\f$
   */
  double _d111;

  /**
   *  The coefficent \f$D_{112}\f$
   */
  double _d112;

  /**
   *  The position of the Adler zero in the \f$I=\frac12\f$ channel
   */
  Energy _s0half;

  /**
   *  The position of the Adler zero in the \f$I=\frac32\f$ channel
   */
  Energy _s0threehalf;

  /**
   *  The position of the pole in the \f$I=\frac12\f$ channel
   */
  Energy _s1;

  /**
   *  Magnitude of the coupling to the pole for the 'initial' production
   */
  Energy _beta;

  /**
   *  Phase of the coupling to the pole for the 'initial' production
   */
  double _theta;

  /**
   *  Coefficient \f$c_{10}\f$ for the production
   */
  double _c10;

  /**
   *  Coefficient \f$c_{20}\f$ for the production
   */
  double _c20;

  /**
   *  Coefficient \f$c_{30}\f$ for the production
   */
  double _c30;

  /**
   *  Coefficient \f$c_{11}\f$ for the production
   */
  double _c11;

  /**
   *  Coefficient \f$c_{21}\f$ for the production
   */
  double _c21;

  /**
   *  Coefficient \f$c_{31}\f$ for the production
   */
  double _c31;

  /**
   *  Coefficient \f$c_{12}\f$ for the production
   */
  double _c12;

  /**
   *  Coefficient \f$c_{22}\f$ for the production
   */
  double _c22;

  /**
   *  Coefficient \f$c_{32}\f$ for the production
   */
  double _c32;

  /**
   *  The phase \f$\gamma_1\f$ for the production
   */
  double _gamma1;

  /**
   *  The phase \f$\gamma_2\f$ for the production
   */
  double _gamma2;

  /**
   *  The phase \f$\gamma_3\f$ for the production
   */
  double _gamma3;

  /**
   *  The mass of the kaon
   */
  Energy _mK;

  /**
   *  The mass of the pion
   */
  Energy _mpi;

  /**
   *  The mass of the \f$\eta'\f$
   */
  Energy _metap;
  //@}

  /**
   *  Complex Amplitudes
   */
  //@{
  /**
   *  Amplitude for the non-resonant component
   */
  Complex _cNR;

  /**
   *  Amplitude for the \f$\kappa\f$ component
   */
  complex<Energy2> _ckappa;

  /**
   *  Amplitude for the \f$K^*(892)\f$
   */
  Complex _cK892;

  /**
   *  Amplitude for the \f$K^*(1410)\f$
   */
  Complex _cK1410;

  /**
   *  Amplitude for the \f$K^*(1680)\f$
   */
  Complex _cK1680;

  /**
   *  Amplitude for the \f$K^*_0(1430)\f$
   */
  complex<Energy2> _cK14300;

  /**
   *  Amplitude for the \f$K^*_2(1430)\f$
   */
  complex<InvEnergy2> _cK14302;
  //@}

  /**
   *  Masses and Widths for the intermediate resonances
   */
  //@{
  /**
   *  Use local values of the parameters
   */
  bool _localparameters;

  /**
   *  Mass of the \f$\kappa\f$
   */
  Energy _mkappa;

  /**
   *  Width of the \f$\kappa\f$
   */
  Energy _wkappa;

  /**
   *  Mass of the \f$K^*(892)\f$
   */
  Energy _mK892;

  /**
   *  Width of the \f$K^*(892)\f$
   */
  Energy _wK892;

  /**
   *  Mass of the \f$K^*(1410)\f$
   */
  Energy _mK1410;

  /**
   *  Width of the \f$K^*(1410)\f$
   */
  Energy _wK1410;

  /**
   *  Mass of the \f$K^*(1680)\f$
   */
  Energy _mK1680;

  /**
   *  Width of the \f$K^*(1680)\f$
   */
  Energy _wK1680;

  /**
   *  Mass of the \f%K^*_0(1430)\f$
   */
  Energy _mK14300;

  /**
   *  Width of the \f%K^*_0(1430)\f$
   */
  Energy _wK14300;

  /**
   *  Mass of the \f%K^*_2(1430)\f$
   */
  Energy _mK14302;

  /**
   *  Width of the \f%K^*_2(1430)\f$
   */
  Energy _wK14302;
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
   *  Parameters for the integration
   */
  //@{
  /**
   *  Maximum weight for the integration
   */
  double _maxwgt;

  /**
   *  Weights for the integration channel
   */
  vector<double> _weights;
  //@}
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of DtoKPiPiFOCUS. */
template <>
struct BaseClassTrait<Herwig::DtoKPiPiFOCUS,1> {
  /** Typedef of the first base class of DtoKPiPiFOCUS. */
  typedef Herwig::DecayIntegrator NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the DtoKPiPiFOCUS class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::DtoKPiPiFOCUS>
  : public ClassTraitsBase<Herwig::DtoKPiPiFOCUS> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig++::DtoKPiPiFOCUS"; }
  /**
   * The name of a file containing the dynamic library where the class
   * DtoKPiPiFOCUS is implemented. It may also include several, space-separated,
   * libraries if the class DtoKPiPiFOCUS depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwWeakCurrents.so HwSMDecay.so"; }
};

/** @endcond */

}

#include "DtoKPiPiFOCUS.icc"

#endif /* HERWIG_DtoKPiPiFOCUS_H */
