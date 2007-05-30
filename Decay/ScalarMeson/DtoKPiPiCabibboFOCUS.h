// -*- C++ -*-
#ifndef HERWIG_DtoKPiPiCabibboFOCUS_H
#define HERWIG_DtoKPiPiCabibboFOCUS_H
//
// This is the declaration of the DtoKPiPiCabibboFOCUS class.
//

#include "Herwig++/Decay/DecayIntegrator.h"
#include "Herwig++/Decay/DecayPhaseSpaceMode.h"
#include "DtoKPiPiCabibboFOCUS.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the DtoKPiPiCabibboFOCUS class.
 *
 * @see \ref DtoKPiPiCabibboFOCUSInterfaces "The interfaces"
 * defined for DtoKPiPiCabibboFOCUS.
 */
class DtoKPiPiCabibboFOCUS: public DecayIntegrator {

public:

  /**
   * The default constructor.
   */
  DtoKPiPiCabibboFOCUS();

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
			   Energy mAB,Energy mres, Energy wres,
			   Energy2 ac, double cost) const;
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
  static ClassDescription<DtoKPiPiCabibboFOCUS> initDtoKPiPiCabibboFOCUS;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  DtoKPiPiCabibboFOCUS & operator=(const DtoKPiPiCabibboFOCUS &);

private:

  /**
   *  Masses and widths of the resonances
   */
  //@{
  /**
   *  Use local values for the masses and widths
   */
  bool _localparameters;

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
   *  Mass of the \f$K^*_0(1430)\f$
   */
  Energy _mK14300;

  /**
   *  Width of the \f$K^*_0(1430)\f$
   */
  Energy _wK14300;

  /**
   *  Mass of the \f$K^*_2(1430)\f$
   */
  Energy _mK14302;

  /**
   *  Width of the \f$K^*_2(1430)\f$
   */
  Energy _wK14302;

  /**
   *  Mass of the \f$\rho(770)\f$
   */
  Energy _mrho770;

  /**
   *  Width of the \f$\rho(770)\f$
   */
  Energy _wrho770;

  /**
   *  Mass of the \f$f_0(980)\f$
   */
  Energy _mf980;

  /**
   *  Width of the \f$f_0(980)\f$
   */
  Energy _wf980;

  /**
   *  Mass of the \f$\rho(1450)\f$
   */
  Energy _mrho1450;

  /**
   *  Width of the \f$\rho(1450)\f$
   */
  Energy _wrho1450;
  //@}

  /**
   *  Amplitudes and Phases for the \f$D^+\to K^+\pi^-\pi^+\f$
   */
  //@{
  /**
   *  Amplitude for \f$\rho(770)\f$
   */
  double _aDrho770;

  /**
   *  Phase for \f$\rho(770)\f$
   */
  double _phiDrho770;

  /**
   *  Amplitude for \f$K^*(892)\f$
   */
  double _aDK892;

  /**
   *  Phase for \f$K^*(892)\f$
   */
  double _phiDK892;

  /**
   *  Amplitude for \f$f_0(980)\f$
   */
  Energy2 _aDf980;

  /**
   *  Phase for \f$f_0(980)\f$
   */
  double _phiDf980;

  /**
   *  Amplitude for \f$K^*_2(1430)\f$
   */
  InvEnergy2 _aDK1430;

  /**
   *  Phase for \f$K^*_2(1430)\f$
   */
  double _phiDK1430;

  /**
   *  Complex Amplitude for \f$\rho(770)\f$
   */
  Complex _cDrho770;

  /**
   *  Complex Amplitude for \f$K^*(892)\f$
   */
  Complex _cDK892;

  /**
   *  Complex Amplitude for \f$f_0(980)\f$
   */
  complex<Energy2> _cDf980;

  /**
   *  Complex Amplitude for \f$K^*_2(1430)\f$
   */
  complex<InvEnergy2> _cDK1430;
  //@}

  /**
   *  Amplitudes and Phases for the \f$D^+_s\to K^+\pi^-\pi^+\f$
   */
  //@{
  /**
   *  Amplitude for the non-resonant component
   */
  double _aDsNR;

  /**
   *  Phase for the non-resonant component
   */
  double _phiDsNR;

  /**
   *  Amplitude for the \f$\rho(770)\f$
   */
  double _aDsrho770;

  /**
   *  Phase for the \f$\rho(770)\f$
   */
  double _phiDsrho770;

  /**
   *  Amplitude for \f$K^*(892)\f$
   */
  double _aDsK892;

  /**
   *  Phase for \f$K^*(892)\f$
   */
  double _phiDsK892;

  /**
   *  Amplitude for \f$K^*(1410)\f$
   */
  double _aDsK1410;

  /**
   *  Phase for \f$K^*(1410)\f$
   */
  double _phiDsK1410;

  /**
   *  Amplitude for \f$K^*_0(1430)\f$
   */
  Energy2 _aDsK1430;

  /**
   *  Phase for \f$K^*_0(1430)\f$
   */
  double _phiDsK1430;

  /**
   *  Amplitude for the \f$\rho(1450)\f$
   */
  double _aDsrho1450;

  /**
   *  Phase for the \f$\rho(1450)\f$
   */
  double _phiDsrho1450;

  /**
   *  Complex Amplitude for the non-resonant component
   */
  Complex _cDsNR;

  /**
   *  Complex Amplitude for the \f$\rho(770)\f$
   */
  Complex _cDsrho770;

  /**
   *  Complex Amplitude for \f$K^*(892)\f$
   */
  Complex _cDsK892;

  /**
   *  Complex Amplitude for \f$K^*(1410)\f$
   */
  Complex _cDsK1410;

  /**
   *  Complex Amplitude for \f$K^*_0(1430)\f$
   */
  complex<Energy2> _cDsK1430;

  /**
   *  Complex Amplitude for the \f$\rho(1450)\f$
   */
  Complex _cDsrho1450;
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
   *  Parameters for the phase-space integration
   */
  //@{
  /**
   *  Maximum weights for the different modes
   */
  vector<double> _maxweight;

  /**
   *  Weights for the different channels
   */
  vector<double> _weights;
  //@}
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of DtoKPiPiCabibboFOCUS. */
template <>
struct BaseClassTrait<Herwig::DtoKPiPiCabibboFOCUS,1> {
  /** Typedef of the first base class of DtoKPiPiCabibboFOCUS. */
  typedef Herwig::DecayIntegrator NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the DtoKPiPiCabibboFOCUS class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::DtoKPiPiCabibboFOCUS>
  : public ClassTraitsBase<Herwig::DtoKPiPiCabibboFOCUS> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig++::DtoKPiPiCabibboFOCUS"; }
  /**
   * The name of a file containing the dynamic library where the class
   * DtoKPiPiCabibboFOCUS is implemented. It may also include several, space-separated,
   * libraries if the class DtoKPiPiCabibboFOCUS depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwWeakCurrents.so HwSMDecay.so"; }
};

/** @endcond */

}

#include "DtoKPiPiCabibboFOCUS.icc"

#endif /* HERWIG_DtoKPiPiCabibboFOCUS_H */
