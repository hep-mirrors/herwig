// -*- C++ -*-
#ifndef HERWIG_DtoKPiPiE791_H
#define HERWIG_DtoKPiPiE791_H
//
// This is the declaration of the DtoKPiPiE791 class.
//

#include "Herwig++/Decay/DecayIntegrator.h"
#include "Herwig++/Decay/DecayPhaseSpaceMode.h"
#include "DtoKPiPiE791.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the DtoKPiPiE791 class.
 *
 * @see \ref DtoKPiPiE791Interfaces "The interfaces"
 * defined for DtoKPiPiE791.
 */
class DtoKPiPiE791: public DecayIntegrator {

public:

  /**
   * The default constructor.
   */
  DtoKPiPiE791();

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
   * angle between the 2 and 3 for the decay \f$D\to(12)3\f$ in the rest frame of
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
			   Energy mAB, double cost, Energy2 ac,
			   Energy mres, Energy wres) const;
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
  static ClassDescription<DtoKPiPiE791> initDtoKPiPiE791;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  DtoKPiPiE791 & operator=(const DtoKPiPiE791 &);

private:

  /**
   *  Switch to control which model is used
   */
  unsigned int _imodel;

  /**
   *  Amplitudes and Phases for the different components
   */
  //@{
  /**
   *  Amplitude for the non-resonant component in model A
   */
  double _aANR;

  /**
   *  Phase for the non-resonant component in model A
   */
  double _phiANR;

  /**
   *  Amplitude for the \f$K^*(892)^-\f$ component in model A
   */
  double _aAK892;

  /**
   *  Phase for the \f$K^*(892)^-\f$ component in model A
   */
  double _phiAK892;

  /**
   *  Amplitude for the \f$K^*_0(1430)^-\f$ component in model A
   */
  Energy2 _aAK14300;

  /**
   *  Phase for the \f$K^*_0(1430)^-\f$ component in model A
   */
  double _phiAK14300;

  /**
   *  Amplitude for the \f$K^*_2(1430)^-\f$ component in model A
   */
  InvEnergy2 _aAK14302;

  /**
   *  Phase for the \f$K^*_2(1430)^-\f$ component in model A
   */
  double _phiAK14302;

  /**
   *  Amplitude for the \f$K^*(1680)^-\f$ component in model A
   */
  double _aAK1680;

  /**
   *  Phase for the \f$K^*(1680)^-\f$ component in model A
   */
  double _phiAK1680;

  /**
   *  Amplitude for the non-resonant component in model B
   */
  double _aBNR;

  /**
   *  Phase for the non-resonant component in model B
   */
  double _phiBNR;

  /**
   *  Amplitude for the \f$K^*(892)^-\f$ component in model B
   */
  double _aBK892;

  /**
   *  Phase for the \f$K^*(892)^-\f$ component in model B
   */
  double _phiBK892;

  /**
   *  Amplitude for the \f$K^*_0(1430)^-\f$ component in model B
   */
  Energy2 _aBK14300;

  /**
   *  Phase for the \f$K^*_0(1430)^-\f$ component in model B
   */
  double _phiBK14300;

  /**
   *  Amplitude for the \f$K^*_2(1430)^-\f$ component in model B
   */
  InvEnergy2 _aBK14302;

  /**
   *  Phase for the \f$K^*_2(1430)^-\f$ component in model B
   */
  double _phiBK14302;

  /**
   *  Amplitude for the \f$K^*(1680)^-\f$ component in model B
   */
  double _aBK1680;

  /**
   *  Phase for the \f$K^*(1680)^-\f$ component in model B
   */
  double _phiBK1680;

  /**
   *  Amplitude for the non-resonant component in model C
   */
  double _aCNR;

  /**
   *  Phase for the non-resonant component in model C
   */
  double _phiCNR;

  /**
   *  Amplitude for the \f$\kappa\f$ component in model C
   */
  Energy2 _aCkappa;

  /**
   *  Phase for the \f$\kappa\f$ component in model C
   */
  double _phiCkappa;

  /**
   *  Amplitude for the \f$K^*(892)^-\f$ component in model C
   */
  double _aCK892;

  /**
   *  Phase for the \f$K^*(892)^-\f$ component in model C
   */
  double _phiCK892;

  /**
   *  Amplitude for the \f$K^*_0(1430)^-\f$ component in model C
   */
  Energy2 _aCK14300;

  /**
   *  Phase for the \f$K^*_0(1430)^-\f$ component in model C
   */
  double _phiCK14300;

  /**
   *  Amplitude for the \f$K^*_2(1430)^-\f$ component in model C
   */
  InvEnergy2 _aCK14302;

  /**
   *  Phase for the \f$K^*_2(1430)^-\f$ component in model C
   */
  double _phiCK14302;

  /**
   *  Amplitude for the \f$K^*(1680)^-\f$ component in model C
   */
  double _aCK1680;

  /**
   *  Phase for the \f$K^*(1680)^-\f$ component in model C
   */
  double _phiCK1680;
  //@}

  /**
   *  Complex amplitudes for the various components
   */
  //@{
  /**
   *  Complex amplitude for the non-resonant term
   */
  Complex _cNR;

  /**
   *  Complex amplitude for the \f$\kappa\f$ component in 
   */
  complex<Energy2> _ckappa;

  /**
   *  Complex amplitude for the \f$K^*(892)^-\f$ component in 
   */
  Complex _cK892;

  /**
   *  Complex amplitude for the \f$K^*_0(1430)^-\f$ component in 
   */
  complex<Energy2> _cK14300;

  /**
   *  Complex amplitude for the \f$K^*_2(1430)^-\f$ component in 
   */
  complex<InvEnergy2> _cK14302;

  /**
   *  Complex amplitude for the \f$K^*(1680)^-\f$ component in 
   */
  Complex _cK1680;
  //@}

  /**
   *  Parameters for the Blatt-Weisskopf form-factors
   */
  //@{
  /**
   *  Radial size for the \f$D^0\f$ for model A
   */
  InvEnergy _rD0A;

  /**
   *  Radial size for the light resonances for model A
   */
  InvEnergy _rresA;

  /**
   *  Radial size for the \f$D^0\f$ for model B
   */
  InvEnergy _rD0B;

  /**
   *  Radial size for the light resonances for model B
   */
  InvEnergy _rresB;

  /**
   *  Radial size for the \f$D^0\f$ for model C
   */
  InvEnergy _rD0C;

  /**
   *  Radial size for the light resonances for model C
   */
  InvEnergy _rresC;
  //@}

  /**
   *  Masses and Widths of the resonances
   */
  //@{
  /**
   *  Use local values for the masses and widths
   */
  bool _localparameters;

  /**
   *  Mass of \f$\kappa\f$
   */
  Energy _mkappa;

  /**
   *  Width of \f$\kappa\f$
   */
  Energy _wkappa;

  /**
   *  Mass of \f$K^*(892)\f$
   */
  Energy _mK892;

  /**
   *  Width of \f$K^*(892)\f$
   */
  Energy _wK892;

  /**
   *  Mass of \f$K^*_0(1430)\f$ actually used 
   */
  Energy _mK1430;

  /**
   *  Width of \f$K^*_0(1430)\f$ actually used
   */
  Energy _wK1430;

  /**
   *  Mass of \f$K^*_0(1430)\f$ in model A
   */
  Energy _mK1430A;

  /**
   *  Width of \f$K^*_0(1430)\f$ in model A
   */
  Energy _wK1430A;

  /**
   *  Mass of \f$K^*_0(1430)\f$ in model B
   */
  Energy _mK1430B;

  /**
   *  Width of \f$K^*_0(1430)\f$ in model B
   */
  Energy _wK1430B;

  /**
   *  Mass of \f$K^*_0(1430)\f$ in model C
   */
  Energy _mK1430C;

  /**
   *  Width of \f$K^*_0(1430)\f$ in model C
   */
  Energy _wK1430C;

  /**
   *  Mass of \f$K^*_2(1430)\f$ 
   */
  Energy _mK14302;

  /**
   *  Width of \f$K^*_2(1430)\f$
   */
  Energy _wK14302;

  /**
   *  Mass of \f$K^*(1680)\f$
   */
  Energy _mK1680;

  /**
   *  Width of \f$K^*(1680)\f$
   */
  Energy _wK1680;
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
 *  base classes of DtoKPiPiE791. */
template <>
struct BaseClassTrait<Herwig::DtoKPiPiE791,1> {
  /** Typedef of the first base class of DtoKPiPiE791. */
  typedef Herwig::DecayIntegrator NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the DtoKPiPiE791 class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::DtoKPiPiE791>
  : public ClassTraitsBase<Herwig::DtoKPiPiE791> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::DtoKPiPiE791"; }
  /**
   * The name of a file containing the dynamic library where the class
   * DtoKPiPiE791 is implemented. It may also include several, space-separated,
   * libraries if the class DtoKPiPiE791 depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwWeakCurrents.so HwSMDecay.so"; }
};

/** @endcond */

}

#include "DtoKPiPiE791.icc"

#endif /* HERWIG_DtoKPiPiE791_H */
