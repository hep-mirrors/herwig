// -*- C++ -*-
//
// DtoKPiPiMarkIII.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_DtoKPiPiMarkIII_H
#define HERWIG_DtoKPiPiMarkIII_H
//
// This is the declaration of the DtoKPiPiMarkIII class.
//

#include "Herwig/Decay/DecayIntegrator.h"
#include "Herwig/Decay/DecayPhaseSpaceMode.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the DtoKPiPiMarkIII class.
 *
 * @see \ref DtoKPiPiMarkIIIInterfaces "The interfaces"
 * defined for DtoKPiPiMarkIII.
 */
class DtoKPiPiMarkIII: public DecayIntegrator {

public:

  /**
   * The default constructor.
   */
  DtoKPiPiMarkIII();
  
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
   * Calculate the amplitude for a resonance
   * @param rho True for rho resonances and false for \f$K^*\f$
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
  Complex amplitude(bool rho, Energy mD, 
			   Energy mA , Energy mB , Energy mC ,
			   Energy mAB, Energy mAC, Energy mBC,
			   Energy mres, Energy wres) const;

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
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  DtoKPiPiMarkIII & operator=(const DtoKPiPiMarkIII &) = delete;

private:

  /**
   *  Amplitudes and phases for \f$D^0\to K^-\pi^+\pi^0\f$
   */
  //@{
  /**
   *  Magnitude of the \f$\rho\f$ component
   */
  double _a1rho;

  /**
   *  Phase of the \f$\rho\f$ component
   */
  double _phi1rho;

  /**
   *  Magnitude of the \f$K^{*-}\f$ component
   */
  double _a1Kstarm;

  /**
   *  Phase of the \f$K^{*-}\f$ component
   */
  double _phi1Kstarm;

  /**
   *  Magnitude of the \f$\bar{K}^{*0}\f$ component
   */
  double _a1Kstar0;

  /**
   *  Phase of the \f$\bar{K}^{*0}\f$ component
   */
  double _phi1Kstar0;

  /**
   *  Magnitude of the non-resonant component
   */
  double _a1NR;

  /**
   *  Phase of the non-resonant component
   */
  double _phi1NR;

  /**
   *  Magnitude of the \f$\rho\f$ component
   */
  Complex _c1rho;

  /**
   *  Magnitude of the \f$K^{*-}\f$ component
   */
  Complex _c1Kstarm;

  /**
   *  Magnitude of the \f$\bar{K}^{*0}\f$ component
   */
  Complex _c1Kstar0;

  /**
   *  Magnitude of the non-resonant component
   */
  Complex _c1NR;
  //@}

  /**
   *  Amplitudes and phases for \f$D^0\to \bar{K}^0\pi^+\pi^-\f$
   */
  //@{
  /**
   *  Magnitude of the \f$\rho^0\f$ component
   */
  double _a2rho;

  /**
   *  Phase of the \f$\rho^0\f$ component
   */
  double _phi2rho;

  /**
   *  Magnitude of the \f$K^{*-}\f$ component
   */
  double _a2Kstar;

  /**
   *  Phase of the \f$K^{*-}\f$ component
   */
  double _phi2Kstar;

  /**
   *  Magnitude of the non-resonant component
   */
  double _a2NR;

  /**
   *  Phase of the non-resonant component
   */
  double _phi2NR;

  /**
   *  Amplitude of the \f$\rho^0\f$ component
   */
  Complex _c2rho;

  /**
   *  Amplitude of the \f$K^{*-}\f$ component
   */
  Complex _c2Kstar;

  /**
   *  Amplitude of the non-resonant component
   */
  Complex _c2NR;
  //@}

  /**
   *  Amplitudes and phases for \f$D^+\to \bar{K}^0\pi^+\pi^0\f$
   */
  //@{
  /**
   *  Magnitude of the \f$\rho\f$ component
   */
  double _a3rho;

  /**
   *  Phase of the \f$\rho\f$ component
   */
  double _phi3rho;

  /**
   *  Magnitude of the \f$\bar{K}^{*0}\f$ component
   */
  double _a3Kstar;

  /**
   *  Phase of the \f$\bar{K}^{*0}\f$ component
   */
  double _phi3Kstar;

  /**
   *  Magnitude of the non-resonant component
   */
  double _a3NR;

  /**
   *  Phase of the non-resonant component
   */
  double _phi3NR;

  /**
   *  Amplitude of the \f$\rho\f$ component
   */
  Complex _c3rho;

  /**
   *  Amplitude of the \f$\bar{K}^{*0}\f$ component
   */
  Complex _c3Kstar;

  /**
   *  Amplitude of the non-resonant component
   */
  Complex _c3NR;
  //@}

  /**
   *  Amplitudes and phases for \f$D^+\to K^-\pi^+\pi^+\f$
   */
  //@{
  /**
   *  Magnitude of the \f$\bar{K}^{*0}\f$ component
   */
  double _a4Kstar;

  /**
   *  Phase of the \f$\bar{K}^{*0}\f$ component
   */
  double _phi4Kstar;

  /**
   *  Magnitude of the non-resonant component
   */
  double _a4NR;

  /**
   *  Phase of the non-resonant component
   */
  double _phi4NR;

  /**
   *  Amplitude of the \f$\bar{K}^{*0}\f$ component
   */
  Complex _c4Kstar;

  /**
   *  Amplitude of the non-resonant component
   */
  Complex _c4NR;
  //@}

  /**
   *  Masses and Widths of the resonances
   */
  //@{
  /**
   *  Use local values of the masses and widths
   */
  bool _localparameters;

  /**
   *  Mass of the \f$\rho^+\f$
   */
  Energy _mrhop;

  /**
   *  Width of the \f$\rho^+\f$
   */
  Energy _wrhop;

  /**
   *  Mass of the \f$\rho^0\f$
   */
  Energy _mrho0;

  /**
   *  Width of the \f$\rho^0\f$
   */
  Energy _wrho0;

  /**
   *  Mass of the \f$K^{*-}\f$
   */
  Energy _mKstarm;

  /**
   *  Width of the \f$K^{*-}\f$
   */
  Energy _wKstarm;

  /**
   *  Mass of the \f$\bar{K}^{*0}\f$
   */
  Energy _mKstar0;

  /**
   *  Width of the \f$\bar{K}^{*0}\f$
   */
  Energy _wKstar0;
  //@}

  /**
   *  The radii of the mesons for the form-factors
   */
  //@{
  /**
   *  \f$\rho\f$ radius
   */
  InvEnergy _rrho;

  /**
   *  \f$K^*\f$ radius
   */
  InvEnergy _rKstar;
  //@}

  /**
   *  Parameters for the phase-space integration
   */
  //@{
  /**
   *  Maximum weights for the different modes
   */
  vector<double> _maxwgt;

  /**
   *  Weights for the different phase-space channels
   */
  vector<double> _weights;
  //@}

  /**
   *  Spin density matrix
   */
  mutable RhoDMatrix _rho;
};

}

#endif /* HERWIG_DtoKPiPiMarkIII_H */
