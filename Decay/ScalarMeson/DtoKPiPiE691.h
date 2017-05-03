// -*- C++ -*-
//
// DtoKPiPiE691.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_DtoKPiPiE691_H
#define HERWIG_DtoKPiPiE691_H
//
// This is the declaration of the DtoKPiPiE691 class.
//

#include "Herwig/Decay/DecayIntegrator.h"
#include "Herwig/Decay/DecayPhaseSpaceMode.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the DtoKPiPiE691 class.
 *
 * @see \ref DtoKPiPiE691Interfaces "The interfaces"
 * defined for DtoKPiPiE691.
 */
class DtoKPiPiE691: public DecayIntegrator {

public:
  /**
   * The default constructor.
   */
  DtoKPiPiE691();
  
  /**
   * Which of the possible decays is required
   * @param cc Is this mode the charge conjugate
   * @param parent The decaying particle
   * @param children The decay products
   */
  virtual int modeNumber(bool & cc, tcPDPtr parent, 
			 const tPDVector & children) const;

  /**
   * Return the matrix element squared for a given mode and phase-space channel.
   * @param ichan The channel we are calculating the matrix element for. 
   * @param part The decaying Particle.
   * @param decay The particles produced in the decay.
   * @param meopt Option for the matrix element
   * @return The matrix element squared for the phase-space configuration.
   */
  double me2( const int ichan,const Particle & part,
	     const ParticleVector & decay, MEOption meopt) const;

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
   */
  double decayAngle(const Lorentz5Momentum & pparent,
		    const Lorentz5Momentum & pres,
		    const Lorentz5Momentum & p1) const {
    Energy2 dot   = pparent*p1,   mREp  = pres*pparent;
    Energy2 mRE1  = pres*p1,      mp2   = pparent.mass2();
    Energy2 mres2 = pres.mass2(), m12   = p1.mass2();
    return (dot*mres2-mREp*mRE1)/
      sqrt((mREp*mREp-mres2*mp2)*(mRE1*mRE1-mres2*m12));
  }

  /**
   * Calculate the amplitude
   * @param ispin The spin of the resonance
   * @param costheta The decay angle
   * @param mAB The off-shell mass of the resonance
   * @param wres The width of the resonance
   * @param mres The on-shell mass of the resonance
   */
  Complex amplitude(int ispin, double costheta,Energy mAB,
		    Energy wres, Energy mres) const {
    double s = 0.;
    switch(ispin) {
    case 0: s = 1.;                    break;
    case 1: s = costheta;              break;
    case 2: s = 1.5*sqr(costheta)-0.5; break;
    default: assert(false);
    }
    Complex bw = sqrt(0.5*wres/GeV/Constants::pi)*GeV/
      (mAB-mres-complex<Energy>(ZERO,0.5*wres));
    return s*bw;
  }
  //@}

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
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();

  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  virtual void doinitrun();
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<DtoKPiPiE691> initDtoKPiPiE691;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  DtoKPiPiE691 & operator=(const DtoKPiPiE691 &);

private:

  /**
   *  Amplitudes and phases for the different components
   */
  //@{
  /**
   *  Amplitude of the non-resonant component for \f$D^+\to K^-\pi^+\pi^+\f$
   */
  double _a1NR;

  /**
   *  Phase of the non-resonant component for \f$D^+\to K^-\pi^+\pi^+\f$
   */
  double _phi1NR;

  /**
   *  Amplitude of the \f$\bar{K}^*(892)^0\f$ component for \f$D^+\to K^-\pi^+\pi^+\f$
   */
  double _a1K892;

  /**
   *  Phase of the \f$\bar{K}^*(892)^0\f$ component for \f$D^+\to K^-\pi^+\pi^+\f$
   */
  double _phi1K892;

  /**
   *  Amplitude of the \f$\bar{K}^*_0(1430)^0\f$ component for \f$D^+\to K^-\pi^+\pi^+\f$
   */
  double _a1K1430;

  /**
   *  Phase of the \f$\bar{K}^*_0(1430)^0\f$ component for \f$D^+\to K^-\pi^+\pi^+\f$
   */
  double _phi1K1430;

  /**
   *  Amplitude of the \f$\bar{K}^*_0(1680)^0\f$ component for \f$D^+\to K^-\pi^+\pi^+\f$
   */
  double _a1K1680;

  /**
   *  Phase of the \f$\bar{K}^*_0(1680)^0\f$ component for \f$D^+\to K^-\pi^+\pi^+\f$
   */
  double _phi1K1680;

  /**
   *  Amplitude of the non-resonant component for \f$D^0\to K^-\pi^+\pi^0\f$
   */
  double _a2NR;

  /**
   *  Phase of the non-resonant component for \f$D^0\to K^-\pi^+\pi^0\f$
   */
  double _phi2NR;

  /**
   *  Amplitude of the \f$\bar{K}^*(892)^0\f$  component for \f$D^0\to K^-\pi^+\pi^0\f$
   */
  double _a2K8920;

  /**
   *  Phase of the \f$\bar{K}^*(892)^0\f$  component for \f$D^0\to K^-\pi^+\pi^0\f$
   */
  double _phi2K8920;

  /**
   *  Amplitude of the \f$K^*(892)^-\f$ component  for \f$D^0\to K^-\pi^+\pi^0\f$
   */
  double _a2K892m;

  /**
   *  Phase of the \f$K^*(892)^-\f$ component  for \f$D^0\to K^-\pi^+\pi^0\f$
   */
  double _phi2K892m;

  /**
   *  Amplitude of the \f$\rho^+\f$ component  for \f$D^0\to K^-\pi^+\pi^0\f$
   */
  double _a2rho;

  /**
   *  Phase of the \f$\rho^+\f$  component for \f$D^0\to K^-\pi^+\pi^0\f$
   */
  double _phi2rho;

  /**
   *  Amplitude of the non-resonant component for \f$D^0\to \bar{K}^0\pi^+\pi^-\f$
   */
  double _a3NR;

  /**
   *  Phase of the non-resonant component for \f$D^0\to \bar{K}^0\pi^+\pi^-\f$
   */
  double _phi3NR;

  /**
   *  Amplitude of the  \f$K^*(892)^-\f$ component for \f$D^0\to \bar{K}^0\pi^+\pi^-\f$
   */
  double _a3K892;

  /**
   *  Phase of the  \f$K^*(892)^-\f$ component for \f$D^0\to \bar{K}^0\pi^+\pi^-\f$
   */
  double _phi3K892;

  /**
   *  Amplitude of the  \f$\rho^0\f$ component for \f$D^0\to \bar{K}^0\pi^+\pi^-\f$
   */
  double _a3rho;

  /**
   *  Phase of the  \f$\rho^0\f$ component for \f$D^0\to \bar{K}^0\pi^+\pi^-\f$
   */
  double _phi3rho;
  //@}

  /**
   *  Complex amplitudes for use in the matrix element
   */
  //@{
  /**
   *  Amplitude of the non-resonant component for \f$D^+\to K^-\pi^+\pi^+\f$
   */
  Complex _c1NR;

  /**
   *  Amplitude of the \f$\bar{K}^*(892)^0\f$ component for \f$D^+\to K^-\pi^+\pi^+\f$
   */
  Complex _c1K892;

  /**
   *  Amplitude of the \f$\bar{K}^*_0(1430)^0\f$ component for \f$D^+\to K^-\pi^+\pi^+\f$
   */
  Complex _c1K1430;

  /**
   *  Amplitude of the \f$\bar{K}^*_0(1680)^0\f$ component for \f$D^+\to K^-\pi^+\pi^+\f$
   */
  Complex _c1K1680;

  /**
   *  Amplitude of the non-resonant component for \f$D^0\to K^-\pi^+\pi^0\f$
   */
  Complex _c2NR;

  /**
   *  Amplitude of the \f$\bar{K}^*(892)^0\f$ for \f$D^0\to K^-\pi^+\pi^0\f$
   */
  Complex _c2K8920;

  /**
   *  Amplitude of the \f$K^*(892)^-\f$ for \f$D^0\to K^-\pi^+\pi^0\f$
   */
  Complex _c2K892m;

  /**
   *  Amplitude of the \f$\rho^+\f$ for \f$D^0\to K^-\pi^+\pi^0\f$
   */
  Complex _c2rho;

  /**
   *  Amplitude of the non-resonant component for \f$D^0\to \bar{K}^0\pi^+\pi^-\f$
   */
  Complex _c3NR;

  /**
   *  Amplitude of the  \f$K^*(892)^-\f$ component for \f$D^0\to \bar{K}^0\pi^+\pi^-\f$
   */
  Complex _c3K892;

  /**
   *  Amplitude of the  \f$\rho^0\f$ component for \f$D^0\to \bar{K}^0\pi^+\pi^-\f$
   */
  Complex _c3rho;
  //@}

  /**
   *  Masses and widths of the various resonances
   */
  //@{
  /**
   *  Use local values for the masses and widths
   */
  bool _localparameters;

  /**
   * Mass of the \f$K^*(892)^0\f$
   */
  Energy _mK8920;

  /**
   * Width of the \f$K^*(892)^0\f$
   */
  Energy _wK8920;

  /**
   * Mass of the \f$K^*(892)^-\f$
   */
  Energy _mK892m;

  /**
   * Width of the \f$K^*(892)^-\f$
   */
  Energy _wK892m;

  /**
   * Mass of the \f$K^*(1680)^0\f$
   */
  Energy _mK1680;

  /**
   * Width of the \f$K^*(1680)^0\f$
   */
  Energy _wK1680;

  /**
   * Mass of the \f$K^*_0(1430)^0\f$
   */
  Energy _mK1430;

  /**
   * Width of the \f$K^*_0(1430)^0\f$
   */
  Energy _wK1430;

  /**
   * Mass of the \f$\rho^0\f$
   */
  Energy _mrho0;

  /**
   * Width of the \f$\rho^0\f$
   */
  Energy _wrho0;

  /**
   * Mass of the \f$\rho^+\f$
   */
  Energy _mrhop;

  /**
   * Width of the \f$\rho^+\f$
   */
  Energy _wrhop;
  //@}

  /**
   *  Parameters for the phase-space integration
   */
  //@{
  /**
   *  Maximum weights for the various modes
   */
  vector<double> _maxwgt;

  /**
   *  Weights for the different integration channels
   */
  vector<double> _weights;
  //@}

  /**
   *  Spin density matrix
   */
  mutable RhoDMatrix _rho;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of DtoKPiPiE691. */
template <>
struct BaseClassTrait<Herwig::DtoKPiPiE691,1> {
  /** Typedef of the first base class of DtoKPiPiE691. */
  typedef Herwig::DecayIntegrator NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the DtoKPiPiE691 class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::DtoKPiPiE691>
  : public ClassTraitsBase<Herwig::DtoKPiPiE691> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::DtoKPiPiE691"; }
  /**
   * The name of a file containing the dynamic library where the class
   * DtoKPiPiE691 is implemented. It may also include several, space-separated,
   * libraries if the class DtoKPiPiE691 depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwSMDecay.so"; }
};

/** @endcond */

}

#endif /* HERWIG_DtoKPiPiE691_H */
