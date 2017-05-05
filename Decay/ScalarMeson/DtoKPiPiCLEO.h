// -*- C++ -*-
//
// DtoKPiPiCLEO.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_DtoKPiPiCLEO_H
#define HERWIG_DtoKPiPiCLEO_H
//
// This is the declaration of the DtoKPiPiCLEO class.
//

#include "Herwig/Decay/DecayIntegrator.h"
#include "Herwig/Decay/DecayPhaseSpaceMode.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The documentation of the DtoKPiPiCLEO class implements the Dalitz plot fits
 * of the CLEO collaboration for \f$D^0\to\bar{K}^0\pi^+\pi^-\f$,
 *  Phys. Rev. Lett. 89 (2002) 251802,
 * and \f$D^0\to K^-\pi^+\pi^0\f$, Phys. Rev. D63 (2001) 092001.
 *
 * @see \ref DtoKPiPiCLEOInterfaces "The interfaces"
 * defined for DtoKPiPiCLEO.
 */
class DtoKPiPiCLEO: public DecayIntegrator {

public:

  /**
   * The default constructor.
   */
  DtoKPiPiCLEO();
  
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
   *  Calculate the amplitude for a resonance
   * @param ispin The spin of the intermediate resonance
   * @param f0 Whether to use the special form of the propagator for the \f$f_0(980)\f$.
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
  Complex amplitude(int ispin, bool f0, Energy mD, 
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
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<DtoKPiPiCLEO> initDtoKPiPiCLEO;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  DtoKPiPiCLEO & operator=(const DtoKPiPiCLEO &);

private:

  /**
   *  Mass, Widths and related parameters
   */
  //@{
  /**
   *  Whether to use local values for the masses and widths or
   *  those from the ParticleData objects
   */
  bool _localparameters;

  /**
   *  Mass of the \f$\omega\f$
   */
  Energy _momega;

  /**
   *  Width of the \f$\omega\f$
   */
  Energy _womega;

  /**
   *  Mass of the \f$f_0(980)\f$
   */
  Energy _mf980;

  /**
   *  Width of the \f$f_0(980)\f$
   */
  Energy _wf980;

  /**
   * \f$g_\pi\f$ coupling for the \f$f_0(980)\f$ width
   */
  double _gpi;

  /**
   *\f$g_K\f$ coupling for the \f$f_0(980)\f$ width
   */
  double _gK;

  /**
   *  Option for handling the width of the \f$f_0(980)\f$
   */
  bool _f0opt;

  /**
   *  Mass of the \f$f_2(1270)\f$
   */
  Energy _mf2;

  /**
   *  Width of the \f$f_2(1270)\f$
   */
  Energy _wf2;

  /**
   *  Mass of the \f$f_0(1370)\f$
   */
  Energy _mf1370;

  /**
   *  Width of the \f$f_0(1370)\f$
   */
  Energy _wf1370;

  /**
   *  Mass of the \f$K_0^*(1430)\f$
   */
  Energy _mK14300;

  /**
   *  Width of the \f$K_0^*(1430)\f$
   */
  Energy _wK14300;

  /**
   *  Mass of the \f$K_2^*(1430)\f$
   */
  Energy _mK14302;

  /**
   *  Width of the \f$K_2^*(1430)\f$
   */
  Energy _wK14302;

  /**
   *  Mass of the \f$K^*(1680)\f$
   */
  Energy _mK1680;

  /**
   *  Width of the \f$K^*(1680)\f$
   */
  Energy _wK1680;

  /**
   *  Mass of the \f$\rho(1700)\f$
   */
  Energy _mrho1700;

  /**
   *  Width of the \f$\rho(1700)\f$
   */
  Energy _wrho1700;

  /**
   *  Mass of the \f$K^{*0}(892)\f$
   */
  Energy _mK8920;

  /**
   *  Width of the \f$K^{*0}(892)\f$
   */
  Energy _wK8920;

  /**
   *  Mass of the \f$K^{*+}(892)\f$ for \f$D^0\to K^-\pi^+\pi^0\f$
   */
  Energy _mK892A;

  /**
   *  Width of the \f$K^{*+}(892)\f$ for \f$D^0\to K^-\pi^+\pi^0\f$
   */
  Energy _wK892A;

  /**
   *  Mass of the \f$K^{*+}(892)\f$ for \f$D^0\to \bar{K}^0\pi^+\pi^-\f$
   */
  Energy _mK892B;

  /**
   *  Width of the \f$K^{*+}(892)\f$ for \f$D^0\to \bar{K}^0\pi^+\pi^-\f$
   */
  Energy _wK892B;

  /**
   *  Mass of the \f$\rho(770)\f$ for \f$D^0\to K^-\pi^+\pi^0\f$
   */
  Energy _mrhoA;

  /**
   *  Width of the \f$\rho(770)\f$ for \f$D^0\to K^-\pi^+\pi^0\f$
   */
  Energy _wrhoA;

  /**
   *  Mass of the \f$\rho(770)\f$ for \f$D^0\to \bar{K}^0\pi^+\pi^-\f$
   */
  Energy _mrhoB;

  /**
   *  Width of the \f$\rho(770)\f$ for \f$D^0\to \bar{K}^0\pi^+\pi^-\f$
   */
  Energy _wrhoB;
  //@}

  /**
   *  Magnitudes and phases of the amplitudes for \f$D^0\to K^-\pi^+\pi^0\f$
   */
  //@{
  /**
   *  Amplitude of the non-resonant component
   */
  double _a1NR;

  /**
   *  Phase of the non=resonant component
   */
  double _phi1NR;

  /**
   *  Amplitude of the \f$\rho^+\f$ component
   */
  double _a1rho;

  /**
   *  Phase of the \f$\rho^+\f$ component
   */
  double _phi1rho;

  /**
   *  Amplitude of the \f$K^{*-}\f$ component
   */
  double _a1Kstarm;

  /**
   *  Phase of the \f$K^{*-}\f$ component
   */
  double _phi1Kstarm;

  /**
   *  Amplitude of the \f$\bar{K}^{*0}\f$ component
   */
  double _a1Kstar0;

  /**
   *  Phase of the \f$\bar{K}^{*0}\f$ component
   */
  double _phi1Kstar0;

  /**
   *  Amplitude for the \f$K_0(1430)^-\f$ component
   */
  Energy2 _a1K1430m;

  /**
   *  Phase for the \f$K_0(1430)^-\f$ component
   */
  double _phi1K1430m;

  /**
   *  Amplitude for the \f$\bar{K}_0(1430)^0\f$ component
   */
  Energy2 _a1K14300;

  /**
   *  Phase for the \f$\bar{K}_0(1430)^0\f$ component
   */
  double _phi1K14300;

  /**
   *  Amplitude for the \f$\rho(1700)^+\f$ component
   */
  double _a1rho1700;

  /**
   * Phase for the \f$\rho(1700)^+\f$ component
   */
  double _phi1rho1700;

  /**
   *  Amplitude of the \f$K^*(1680)^-\f$ component
   */
  double _a1K1680;

  /**
   *  Phase of the \f$K^*(1680)^-\f$ component
   */
  double _phi1K1680;

  /**
   *  Complex amplitude of the non-resonant component
   */
  Complex _c1NR;

  /**
   *  Complex amplitude of the \f$\rho^+\f$ component
   */
  Complex _c1rho;

  /**
   *  Complex amplitude of the \f$K^{*-}\f$ component
   */
  Complex _c1Kstarm;

  /**
   *  Complex amplitude of the \f$\bar{K}^{*0}\f$ component
   */
  Complex _c1Kstar0;

  /**
   *  Complex amplitude for the \f$K_0(1430)^-\f$ component
   */
  complex<Energy2> _c1K1430m;

  /**
   *  Complex amplitude for the \f$\bar{K}_0(1430)^0\f$ component
   */
  complex<Energy2> _c1K14300;

  /**
   *  Complex amplitude for the \f$\rho(1700)^+\f$ component
   */
  Complex _c1rho1700;

  /**
   *  Complex amplitude of the \f$K^*(1680)^-\f$ component
   */
  Complex _c1K1680;
  //@}

  /**
   *  Magnitudes and phases of the amplitudes for \f$D^0\to \bar{K}^0\pi^+\pi^-\f$ 
   */
  //@{
  /**
   *  Amplitude for the \f$K^{*+}\f$
   */
  double _a2Kstarp;

  /**
   *  Phase for the \f$K^{*+}\f$
   */
  double _phi2Kstarp;

  /**
   *  Amplitude for the \f$\rho^0(770)\f$
   */
  double _a2rho;

  /**
   *  Phase for the \f$\rho^0(770)\f$
   */
  double _phi2rho;

  /**
   *  Amplitude for the \f$\omega\f$
   */
  double _a2omega;

  /**
   *  Phase for the \f$\omega\f$
   */
  double _phi2omega;

  /**
   *  Amplitude for the \f$K^{*-}\f$
   */
  double _a2Kstarm;

  /**
   *  Phase for the \f$K^{*-}\f$
   */
  double _phi2Kstarm;

  /**
   *  Amplitude for the \f$f_0(980)\f$
   */
  Energy2 _a2f980;

  /**
   *  Phase for the \f$f_0(980)\f$
   */
  double _phi2f980;

  /**
   *  Amplitude for the \f$f_2(1270)\f$
   */
  InvEnergy2 _a2f2;

  /**
   *  Phase for the \f$f_2(1270)\f$
   */
  double _phi2f2;

  /**
   *  Amplitude for the \f$f_0(1370)\f$
   */
  Energy2 _a2f1370;

  /**
   *  Phase for the \f$f_0(1370)\f$
   */
  double _phi2f1370;

  /**
   *  Amplitude for the \f$K^*_0(1430)^-\f$
   */
  Energy2 _a2K14300;

  /**
   *  Phase for the \f$K^*_0(1430)^-\f$
   */
  double _phi2K14300;

  /**
   *  Amplitude for the \f$K^*_2(1430)^-\f$
   */
  InvEnergy2 _a2K14302;

  /**
   *  Phase for the \f$K^*_2(1430)^-\f$
   */
  double _phi2K14302;

  /**
   *  Amplitude for the \f$K^*(1680)^-\f$
   */
  double _a2K1680;

  /**
   *  Phase for the \f$K^*(1680)^-\f$
   */
  double _phi2K1680;

  /**
   *  Amplitude of the non-resonant component
   */
  double _a2NR;

  /**
   *  Phase of the non=resonant component
   */
  double _phi2NR;

  /**
   *  Complex amplitude for the \f$K^{*+}\f$
   */
  Complex _c2Kstarp;

  /**
   *  Complex amplitude for the \f$\rho^0(770)\f$
   */
  Complex _c2rho;

  /**
   *  Complex amplitude for the \f$\omega\f$
   */
  Complex _c2omega;

  /**
   *  Complex amplitude for the \f$K^{*-}\f$
   */
  Complex _c2Kstarm;

  /**
   *  Complex amplitude for the \f$f_0(980)\f$
   */
  complex<Energy2> _c2f980;

  /**
   *  Complex amplitude for the \f$f_2(1270)\f$
   */
  complex<InvEnergy2> _c2f2;

  /**
   *  Complex amplitude for the \f$f_0(1370)\f$
   */
  complex<Energy2> _c2f1370;

  /**
   *  Complex amplitude for the \f$K^*_0(1430)^-\f$
   */
  complex<Energy2> _c2K14300;

  /**
   *  Complex amplitude for the \f$K^*_2(1430)^-\f$
   */
  complex<InvEnergy2> _c2K14302;

  /**
   *  Complex amplitude for the \f$K^*(1680)^-\f$
   */
  Complex _c2K1680;

  /**
   *  Complex amplitude of the non-resonant component
   */
  Complex _c2NR;
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
   *  Maximum weights for the modes
   */
  vector<double> _maxwgt;

  /**
   *  Weights for the channels
   */
  vector<double> _weights;
  //@}

  /**
   *  Masses for the \f$f_0(980)\f$ Breit-Wigner
   */
  //@{
  /**
   *  The pion mass
   */
  Energy _mpi;

  /**
   *  The charged kaon mass
   */
  Energy _mkp;

  /**
   *  The neutral kaon mass
   */
  Energy _mk0;
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
 *  base classes of DtoKPiPiCLEO. */
template <>
struct BaseClassTrait<Herwig::DtoKPiPiCLEO,1> {
  /** Typedef of the first base class of DtoKPiPiCLEO. */
  typedef Herwig::DecayIntegrator NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the DtoKPiPiCLEO class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::DtoKPiPiCLEO>
  : public ClassTraitsBase<Herwig::DtoKPiPiCLEO> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::DtoKPiPiCLEO"; }
  /**
   * The name of a file containing the dynamic library where the class
   * DtoKPiPiCLEO is implemented. It may also include several, space-separated,
   * libraries if the class DtoKPiPiCLEO depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwSMDecay.so"; }
};

/** @endcond */

}

#endif /* HERWIG_DtoKPiPiCLEO_H */
